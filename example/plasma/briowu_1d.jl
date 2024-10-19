using KitBase
using KitBase.ProgressMeter: @showprogress

cd(@__DIR__)
ks, ctr, face, t = initialize("briowu_1d.txt")

dt = timestep(ks, ctr, t)
nt = Int(floor(ks.set.maxTime / dt)) + 1
res = zeros(5, 2)

@showprogress for iter in 1:nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt; mode=:kcu, isPlasma=true)
    update!(ks, ctr, face, dt, res; coll=:bgk, bc=:extra, isMHD=true)
end

sol = zeros(ks.ps.nx, 10, 2)
for i in 1:ks.ps.nx
    sol[i, 1, 1] = ctr[i].prim[1, 1]
    sol[i, 1, 2] = ctr[i].prim[1, 2] / ks.gas.me
    sol[i, 2:4, 1] .= ctr[i].prim[2:4, 1]
    sol[i, 2:4, 2] .= ctr[i].prim[2:4, 2]
    sol[i, 5, 1] = 1.0 / ctr[i].prim[5, 1]
    sol[i, 5, 2] = ks.gas.me / ctr[i].prim[5, 2]

    sol[i, 6, 1] = ctr[i].B[2]
    sol[i, 6, 2] = ctr[i].E[1]
end

using Plots
plot(ks.ps.x[1:ks.ps.nx], sol[:, 1, :])
plot(ks.ps.x[1:ks.ps.nx], sol[:, 6, :])
