using Kinetic, Plots
using Kinetic.KitBase.ProgressMeter: @showprogress

cd(@__DIR__)

ks = SolverSet("mixture_shock.txt")
ctr, face = init_fvm(ks)

begin
    iter = 0
    res = zeros(3)
    simTime = 0.0
    dt = Kinetic.timestep(ks, ctr, simTime)
    nt = Int(floor(ks.set.maxTime / dt))
end

res = zeros(3, 2)
@showprogress for iter = 1:nt
    Kinetic.reconstruct!(ks, ctr)
    Kinetic.evolve!(ks, ctr, face, dt)
    Kinetic.update!(ks, ctr, face, dt, res; bc=:fix)
end

sol = zeros(ks.ps.nx, 6)
for i in axes(sol, 1)
    sol[i, 1:3] .= ctr[i].prim[:, 1]
    sol[i, 4:6] .= ctr[i].prim[:, 2]
end

plot(ks.ps.x[1:ks.ps.nx], sol[:, 1] ./ ks.gas.mi)
plot!(ks.ps.x[1:ks.ps.nx], sol[:, 4] ./ ks.gas.me)
