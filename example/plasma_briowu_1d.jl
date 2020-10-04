using Revise, ProgressMeter, Kinetic

cd(@__DIR__)
ks, ctr, face, simTime = Kinetic.initialize("briowu_1d.txt")

dt = Kinetic.timestep(ks, ctr, simTime)
nt = Int(floor(ks.set.maxTime / dt))+1
res = zeros(5, 2)

KS = ks
@showprogress for iter in 1:nt
    #dt = timestep(KS, ctr, simTime)
    Kinetic.reconstruct!(ks, ctr)
    Kinetic.evolve!(ks, ctr, face, dt; mode=:kcu, isPlasma=true)
    #=
    @inbounds Threads.@threads for i = 1:KS.pSpace.nx+1
        flux_kfvs!(
            face[i].fw,
            face[i].fh0,
            face[i].fh1,
            face[i].fh2,
            face[i].fh3,
            ctr[i-1].h0 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh0,
            ctr[i-1].h1 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh1,
            ctr[i-1].h2 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh2,
            ctr[i-1].h3 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh3,
            ctr[i].h0 .- 0.5 .* ctr[i].dx .* ctr[i].sh0,
            ctr[i].h1 .- 0.5 .* ctr[i].dx .* ctr[i].sh1,
            ctr[i].h2 .- 0.5 .* ctr[i].dx .* ctr[i].sh2,
            ctr[i].h3 .- 0.5 .* ctr[i].dx .* ctr[i].sh3,
            KS.vSpace.u,
            KS.vSpace.weights,
            dt,
            ctr[i-1].sh0,
            ctr[i-1].sh1,
            ctr[i-1].sh2,
            ctr[i-1].sh3,
            ctr[i].sh0,
            ctr[i].sh1,
            ctr[i].sh2,
            ctr[i].sh3,
        )
    end
    @inbounds Threads.@threads for i = 1:KS.pSpace.nx+1
        flux_em!(
            face[i].femL,
            face[i].femR,
            ctr[i-2].E,
            ctr[i-2].B,
            ctr[i-1].E,
            ctr[i-1].B,
            ctr[i].E,
            ctr[i].B,
            ctr[i+1].E,
            ctr[i+1].B,
            ctr[i-1].ϕ,
            ctr[i].ϕ,
            ctr[i-1].ψ,
            ctr[i].ψ,
            ctr[i-1].dx,
            ctr[i].dx,
            KS.gas.Ap,
            KS.gas.An,
            KS.gas.D,
            KS.gas.sol,
            KS.gas.χ,
            KS.gas.ν,
            dt,
        )
    end
    =#
    #Kinetic.update!(ks, ctr, face, dt, res; coll=:bgk, bc=:extra, isMHD=true)
    # it's equivalent to the following process
    sumRes = zeros(5, 2)
    sumAvg = zeros(5, 2)
    Threads.@threads for i = 1:ks.pSpace.nx
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
    
    for i = 1:2
        ctr[1-i].w .= ctr[1].w
        ctr[1-i].prim .= ctr[1].prim
        ctr[1-i].h0 .= ctr[1].h0
        ctr[1-i].h1 .= ctr[1].h1
        ctr[1-i].h2 .= ctr[1].h2
        ctr[1-i].E .= ctr[1].E
        ctr[1-i].B .= ctr[1].B
        ctr[1-i].ϕ = deepcopy(ctr[1].ϕ)
        ctr[1-i].ψ = deepcopy(ctr[1].ψ)
        ctr[1-i].lorenz .= ctr[1].lorenz

        ctr[KS.pSpace.nx+i].w .= ctr[KS.pSpace.nx].w
        ctr[KS.pSpace.nx+i].prim .= ctr[KS.pSpace.nx].prim
        ctr[KS.pSpace.nx+i].h0 .= ctr[KS.pSpace.nx].h0
        ctr[KS.pSpace.nx+i].h1 .= ctr[KS.pSpace.nx].h1
        ctr[KS.pSpace.nx+i].h2 .= ctr[KS.pSpace.nx].h2
        ctr[KS.pSpace.nx+i].E .= ctr[KS.pSpace.nx].E
        ctr[KS.pSpace.nx+i].B .= ctr[KS.pSpace.nx].B
        ctr[KS.pSpace.nx+i].ϕ = deepcopy(ctr[KS.pSpace.nx].ϕ)
        ctr[KS.pSpace.nx+i].ψ = deepcopy(ctr[KS.pSpace.nx].ψ)
        ctr[KS.pSpace.nx+i].lorenz .= ctr[KS.pSpace.nx].lorenz
    end
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