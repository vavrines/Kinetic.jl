using Kinetic, OrdinaryDiffEq, SciMLSensitivity, Plots
using KitML.Solaris
using KitML.Solaris.Flux: Adam, throttle

set = Setup(
    case = "shock",
    space = "1d1f1v",
    collision = "shakhov",
    interpOrder = 1,
    maxTime = 10,
)
ps = PSpace1D(-25, 25, 50, 0)
gas = Gas(Kn = 1, Ma = 2, Pr = 2 / 3, K = 0, γ = 3)
vs = VSpace1D(-10, 10, 72)
fw, ff, bc, p = Kinetic.KitBase.ib_rh(set, ps, vs, gas)
ib = IB1F(fw, ff, bc, p)
ks = SolverSet(set, ps, vs, gas, ib)
ctr, face = init_fvm(ks)

dt = timestep(ks, ctr, 0.0)
sumRes = zeros(3)
for iter = 1:123
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, sumRes)
end

X = Array{Float32}(undef, ks.vSpace.nu, ks.pSpace.nx)
for i = 1:ks.pSpace.nx
    X .= ctr[i].f
end

function shakhov!(df, f, p, t)
    M, tau = p
    df .= (M .- f) ./ tau
end

M = Array{Float32}(undef, vs.nu, size(X, 2))
SM = Array{Float32}(undef, vs.nu, size(X, 2))
τ = Array{Float32}(undef, 1, size(X, 2))
for i in axes(X, 2)
    M[:, i] .= maxwellian(ks.vSpace.u, ctr[i].prim)
    q = heat_flux(ctr[i].f, ctr[i].prim, ks.vSpace.u, ks.vSpace.weights)
    S = shakhov(ks.vSpace.u, M[:, i], q, ctr[i].prim, ks.gas.Pr)
    SM[:, i] .= M[:, i] .+ S
    τ[1, i] = vhs_collision_time(ctr[i].prim, ks.gas.μᵣ, ks.gas.ω)
end
P = [SM, τ]

tspan = (0, dt)
prob = ODEProblem(shakhov!, X, tspan, P)
Y = solve(prob, Euler(), dt = dt) |> Array

model_univ = FnChain(FnDense(vs.nu, vs.nu * 2, tanh), FnDense(vs.nu * 2, vs.nu))
p_model = init_params(model_univ)

function dfdt(f, p, t)
    df = (M .- f) ./ τ .+ model_univ(f .- M, p)
end

prob_ube = ODEProblem(dfdt, X, tspan, p_model)

function loss(p)
    sol_ube = solve(prob_ube, Euler(), u0 = X, p = p, dt = dt)
    loss = sum(abs2, Array(sol_ube) .- Y)
    return loss
end

cb = function (p, l)
    display(l)
    return false
end

res = sci_train(loss, p_model, Adam(), cb = throttle(cb, 1), maxiters = 200)

sol = solve(prob_ube, Euler(), u0 = X, p = res.u, dt = dt)

contour(ks.ps.x, ks.vs.u, sol.u[end])
