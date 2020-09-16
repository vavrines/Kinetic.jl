using Revise, ProgressMeter, Kinetic

cd(@__DIR__)
ks, ctr, face, simTime = Kinetic.initialize("briowu_1d.txt")

dt = Kinetic.timestep(KS, ctr, simTime)
nt = Int(floor(ks.set.maxTime / dt))+1
res = zeros(5, 2)

@showprogress for iter in 1:nt
    #dt = timestep(KS, ctr, simTime)
    #Kinetic.reconstruct!(ks, ctr)
    Kinetic.evolve!(ks, ctr, face, dt; mode=:kfvs)
    Kinetic.update!(ks, ctr, face, dt, res)
    # it's equivalent to the following process
    #=
    sumRes = zeros(5, 2)
    sumAvg = zeros(5, 2)
    Threads.@threads for i = 1:KS.pSpace.nx
        Kinetic.step!(
            ks,
            face[i],
            ctr[i],
            face[i+1],
            dt,
            sumRes,
            sumAvg,
        )
    end
    =#
end

sol = zeros(ks.pSpace.nx, 10, 2)
for i in 1:ks.pSpace.nx
    sol[i, 1, 1] = ctr[i].prim[1,1]
    sol[i, 1, 2] = ctr[i].prim[1,2] / ks.gas.me
    sol[i, 2:4, 1] .= ctr[i].prim[2:4,1]
    sol[i, 2:4, 2] .= ctr[i].prim[2:4,2]
    sol[i, 5, 1] = 1. / ctr[i].prim[5,1]
    sol[i, 5, 2] = ks.gas.me / ctr[i].prim[5,2]

    sol[i, 6, 1] = ctr[i].B[2]
    sol[i, 6, 2] = ctr[i].E[1]
end
using Plots
plot(ks.pSpace.x[1:ks.pSpace.nx], sol[:,1,:])
plot(ks.pSpace.x[1:ks.pSpace.nx], sol[:,6,:])

using JLD2
@save "refsol.jld2" ks ctr