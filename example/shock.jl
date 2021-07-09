# This is a demo for solving shock structure problem with neural networks

using Kinetic
using KitML.DiffEqFlux
using OrdinaryDiffEq
using Plots

set = Setup(
    "gas",
    "shock",
    "1d1f1v",
    "kfvs",
    "shakhov",
    1,
    1,
    "vanleer",
    "fix",
    0.5,
    10.0,
)

ps = PSpace1D(-25, 25, 50, 0)
gas = Gas(1.0, 2.0, 2/3, 0.0, 3.0, 0.81, 1.0, 0.5)
vs = VSpace1D(-10, 10, 72)
wL, primL, fL, bcL, wR, primR, fR, bcR = ib_rh(gas.Ma, gas.γ, vs.u)
ib = IB1F(wL, primL, fL, bcL, wR, primR, fR, bcR)
ks = SolverSet(set, ps, vs, gas, ib, pwd())

ctr = Array{ControlVolume1D1F}(undef, ks.pSpace.nx)
face = Array{Interface1D1F}(undef, ks.pSpace.nx + 1)
for i in eachindex(ctr)
    if i <= ks.pSpace.nx ÷ 2
        ctr[i] = ControlVolume1D1F(
            ks.pSpace.x[i],
            ks.pSpace.dx[i],
            Float32.(ks.ib.wL),
            Float32.(ks.ib.primL),
            Float32.(ks.ib.fL),
        )
    else
        ctr[i] = ControlVolume1D1F(
            ks.pSpace.x[i],
            ks.pSpace.dx[i],
            Float32.(ks.ib.wR),
            Float32.(ks.ib.primR),
            Float32.(ks.ib.fR),
        )
    end
end
for i = 1:ks.pSpace.nx+1
    face[i] = Interface1D1F(ks.ib.wL, ks.ib.fL)
end

dt = timestep(ks, ctr, 0.0)
sumRes = zeros(3)
for iter = 1:123
    Kinetic.evolve!(ks, ctr, face, dt)
    Kinetic.update!(ks, ctr, face, dt, sumRes)
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
    S = shakhov(ks.vSpace.u, M[:,i], q, ctr[i].prim, ks.gas.Pr)
    SM[:,i] .= M[:,i] .+ S
    τ[1, i] = vhs_collision_time(ctr[i].prim, ks.gas.μᵣ, ks.gas.ω)
end
P = [SM, τ]

tspan = (0, dt)
prob = ODEProblem(shakhov!, X, tspan, P)
Y = solve(prob, Euler(), dt=dt) |> Array

model_univ = FastChain(FastDense(vs.nu, vs.nu * 2, tanh), FastDense(vs.nu * 2, vs.nu))

p_model = initial_params(model_univ)

function dfdt(f, p, t)
    df = (M .- f) ./ τ .+ model_univ(f .- M, p)
end

prob_ube = ODEProblem(dfdt, X, tspan, p_model)

function loss(p)
    sol_ube = solve(prob_ube, Euler(), u0=X, p=p, dt=dt)
    loss = sum(abs2, Array(sol_ube) .- Y)
    return loss
end

cb = function (p, l)
    display(l)
    return false
end

res = sci_train(loss, p_model, ADAM(), cb=Flux.throttle(cb, 1), maxiters=200)

sol = solve(prob_ube, Euler(), u0=X, p=res.u, dt=dt)

contour(ks.ps.x, ks.vs.u, sol.u[end])
