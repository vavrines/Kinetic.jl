using KitBase, Plots

set = Setup(;
    case="sod",
    space="1d0f0v",
    flux="godunov",
    limiter="minmod",
    cfl=0.3,
    maxTime=0.1,
)
# other avaiable flux solvers include Roe, HLL, HLLC, etc.
ps = PSpace1D(0.0, 1.0, 300, 1)
vs = nothing
gas = Gas(; Kn=1e-4, K=4.0, γ=7 / 5)
fw, bc, p = KB.config_ib(set, ps, vs, gas)

# arbitrary initial condition
primL = [1.0, 0.3, 1.2]
primR = [0.3, -0.1, 2.0]
wL = prim_conserve(primL, gas.γ)
wR = prim_conserve(primR, gas.γ)
p = (x0=ps.x0, x1=ps.x1, wL=wL, wR=wR, primL=primL, primR=primR, γ=gas.γ)
fw = function (x, p)
    if x <= (p.x0 + p.x1) / 2
        return p.wL
    else
        return p.wR
    end
end
bc = function (x, p)
    if x <= (p.x0 + p.x1) / 2
        return p.primL
    else
        return p.primR
    end
end

ib = IB(fw, bc, p)
ks = SolverSet(set, ps, vs, gas, ib)
ctr, face = init_fvm(ks)

t = 0.0
dt = KB.timestep(ks, ctr, t)
nt = ks.set.maxTime ÷ dt |> Int
res = zeros(3)

for iter in 1:nt
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, res)
end

# low-level godunov flux
for iter in 1:nt
    for i in 1:ps.nx+1
        flux_godunov!(face[i].fw, ctr[i-1].w, ctr[i].w, gas.γ, dt)
    end
    update!(ks, ctr, face, dt, res)
end

sol = extract_sol(ks, ctr)
_sole = KB.sample_riemann_solution(
    (primL[1], primL[2], 0.5 * primL[1] / primL[end]),
    (primR[1], primR[2], 0.5 * primR[1] / primR[end]),
    ps.x[1:end-1] .- 0.5,
    0.2,
    gas.γ,
)
sole = [_sole[i][1] for i in 1:length(_sole)]

plot(ps.x[1:end-1], sol[:, 1])
plot!(ps.x[1:end-1], sole[:, 1])
