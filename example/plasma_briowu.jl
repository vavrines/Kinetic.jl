using Revise, ProgressMeter
#using Kinetic
begin
    using Dates
    using LinearAlgebra
    using FastGaussQuadrature
    using OffsetArrays
    using SpecialFunctions
    using FFTW
    using OrdinaryDiffEq
    using FileIO
    using JLD2
    using Plots
    using PyCall
    using ProgressMeter
    include("/home/vavrines/Coding/Kinetic.jl/src/Config/config.jl")
    include("/home/vavrines/Coding/Kinetic.jl/src/Data/data.jl")
    include("/home/vavrines/Coding/Kinetic.jl/src/Flux/flux.jl")
    include("/home/vavrines/Coding/Kinetic.jl/src/Geometry/geometry.jl")
    include("/home/vavrines/Coding/Kinetic.jl/src/IO/io.jl")
    include("/home/vavrines/Coding/Kinetic.jl/src/Math/math.jl")
    include("/home/vavrines/Coding/Kinetic.jl/src/Phase/phase.jl")
    include("/home/vavrines/Coding/Kinetic.jl/src/Reconstruction/reconstruction.jl")
    include("/home/vavrines/Coding/Kinetic.jl/src/Solver/solver.jl")
    include("/home/vavrines/Coding/Kinetic.jl/src/Theory/theory.jl")
end

cd(@__DIR__)
ks, ctr, face, simTime = initialize("config.txt")
KS = ks

dt = timestep(KS, ctr, simTime)
nt = Int(floor(ks.set.maxTime / dt))+1
res = zeros(5, 2)

@showprogress for iter in 1:10000
    #dt = timestep(KS, ctr, simTime)
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt; mode=:kfvs)
    update!(ks, ctr, face, dt, res)
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
end
import Plots
Plots.plot(ks.pSpace.x[1:ks.pSpace.nx], sol[:,1,1])
Plots.plot!(ks.pSpace.x[1:ks.pSpace.nx], sol[:,1,2])

Plots.plot(ks.pSpace.x[1:ks.pSpace.nx], sol[:,6,1])

@save "refsol.jld2" ks ctr