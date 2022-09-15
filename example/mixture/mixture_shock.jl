using KitBase, Plots
using KitBase.ProgressMeter: @showprogress

cd(@__DIR__)

ks = SolverSet("mixture_shock.txt")
ctr, face = init_fvm(ks)

begin
    iter = 0
    res = zeros(3, 2)
    t = 0.0
    dt = timestep(ks, ctr, t)
    nt = Int(floor(ks.set.maxTime / dt))
end

@showprogress for iter = 1:nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, res; bc=:fix)
end

sol = zeros(ks.ps.nx, 6)
for i in axes(sol, 1)
    sol[i, 1:3] .= ctr[i].prim[:, 1]
    sol[i, 4:6] .= ctr[i].prim[:, 2]
end

plot(ks.ps.x[1:ks.ps.nx], sol[:, 1] ./ ks.gas.mi)
plot!(ks.ps.x[1:ks.ps.nx], sol[:, 4] ./ ks.gas.me)
