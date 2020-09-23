using Revise#, Kinetic

cd(@__DIR__)
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
    include("../src/Data/data.jl")
    include("../src/IO/io.jl")
    include("../src/Math/math.jl")
    include("../src/Geometry/geometry.jl")
    include("../src/Theory/theory.jl")
    include("../src/Phase/phase.jl")
    include("../src/Reconstruction/reconstruction.jl")
    include("../src/Flux/flux.jl")
    include("../src/Config/config.jl")
    include("../src/Solver/solver.jl")
end

function step_tst!(
    KS::SolverSet,
    faceL::Interface1D4F,
    cell::ControlVolume1D4F,
    faceR::Interface1D4F,
    dt::AbstractFloat,
)
    #--- update conservative flow variables: step 1 ---#
    # w^n
    w_old = deepcopy(cell.w)
    prim_old = deepcopy(cell.prim)

    # flux -> w^{n+1}
    @. cell.w += (faceL.fw - faceR.fw) / cell.dx
    cell.prim .= mixture_conserve_prim(cell.w, KS.gas.γ)

    #--- update electromagnetic variables ---#
    # flux -> E^{n+1} & B^{n+1}
    cell.E[1] -= dt * (faceL.femR[1] + faceR.femL[1]) / cell.dx
    cell.E[2] -= dt * (faceL.femR[2] + faceR.femL[2]) / cell.dx
    cell.E[3] -= dt * (faceL.femR[3] + faceR.femL[3]) / cell.dx
    cell.B[1] -= dt * (faceL.femR[4] + faceR.femL[4]) / cell.dx
    cell.B[2] -= dt * (faceL.femR[5] + faceR.femL[5]) / cell.dx
    cell.B[3] -= dt * (faceL.femR[6] + faceR.femL[6]) / cell.dx
    cell.ϕ -= dt * (faceL.femR[7] + faceR.femL[7]) / cell.dx
    cell.ψ -= dt * (faceL.femR[8] + faceR.femL[8]) / cell.dx

    for i = 1:3
        if 1 ∈ vcat(isnan.(cell.E), isnan.(cell.B))
            @warn "NaN electromagnetic update"
        end
    end

    # source -> U^{n+1}, E^{n+1} and B^{n+1}
    mr = KS.gas.mi / KS.gas.me
    A, b =
        em_coefficients(cell.prim, cell.E, cell.B, mr, KS.gas.lD, KS.gas.rL, dt)
    x = A \ b

    #--- calculate lorenz force ---#
    cell.lorenz[1, 1] =
        0.5 * (
            x[1] + cell.E[1] + (cell.prim[3, 1] + x[5]) * cell.B[3] -
            (cell.prim[4, 1] + x[6]) * cell.B[2]
        ) / KS.gas.rL
    cell.lorenz[2, 1] =
        0.5 * (
            x[2] + cell.E[2] + (cell.prim[4, 1] + x[6]) * cell.B[1] -
            (cell.prim[2, 1] + x[4]) * cell.B[3]
        ) / KS.gas.rL
    cell.lorenz[3, 1] =
        0.5 * (
            x[3] + cell.E[3] + (cell.prim[2, 1] + x[4]) * cell.B[2] -
            (cell.prim[3, 1] + x[5]) * cell.B[1]
        ) / KS.gas.rL
    cell.lorenz[1, 2] =
        -0.5 *
        (
            x[1] + cell.E[1] + (cell.prim[3, 2] + x[8]) * cell.B[3] -
            (cell.prim[4, 2] + x[9]) * cell.B[2]
        ) *
        mr / KS.gas.rL
    cell.lorenz[2, 2] =
        -0.5 *
        (
            x[2] + cell.E[2] + (cell.prim[4, 2] + x[9]) * cell.B[1] -
            (cell.prim[2, 2] + x[7]) * cell.B[3]
        ) *
        mr / KS.gas.rL
    cell.lorenz[3, 2] =
        -0.5 *
        (
            x[3] + cell.E[3] + (cell.prim[2, 2] + x[7]) * cell.B[2] -
            (cell.prim[3, 2] + x[8]) * cell.B[1]
        ) *
        mr / KS.gas.rL

    cell.E[1] = x[1]
    cell.E[2] = x[2]
    cell.E[3] = x[3]

    #--- update conservative flow variables: step 2 ---#
    cell.prim[2, 1] = x[4]
    cell.prim[3, 1] = x[5]
    cell.prim[4, 1] = x[6]
    cell.prim[2, 2] = x[7]
    cell.prim[3, 2] = x[8]
    cell.prim[4, 2] = x[9]

    cell.w .= mixture_prim_conserve(cell.prim, KS.gas.γ)

    #--- update particle distribution function ---#
    # flux -> f^{n+1}
    @. cell.h0 += (faceL.fh0 - faceR.fh0) / cell.dx
    @. cell.h1 += (faceL.fh1 - faceR.fh1) / cell.dx
    @. cell.h2 += (faceL.fh2 - faceR.fh2) / cell.dx
    @. cell.h3 += (faceL.fh3 - faceR.fh3) / cell.dx

    # force -> f^{n+1} : step 1
    for j in axes(cell.h0, 2)
        _h0 = @view cell.h0[:, j]
        _h1 = @view cell.h1[:, j]
        _h2 = @view cell.h2[:, j]
        _h3 = @view cell.h3[:, j]

        shift_pdf!(_h0, cell.lorenz[1, j], KS.vSpace.du[1, j], dt)
        shift_pdf!(_h1, cell.lorenz[1, j], KS.vSpace.du[1, j], dt)
        shift_pdf!(_h2, cell.lorenz[1, j], KS.vSpace.du[1, j], dt)
        shift_pdf!(_h3, cell.lorenz[1, j], KS.vSpace.du[1, j], dt)
    end

    # force -> f^{n+1} : step 2
    for k in axes(cell.h1, 3)
        @. cell.h3[:, k] +=
            2.0 * dt * cell.lorenz[2, k] * cell.h1[:, k] +
            (dt * cell.lorenz[2, k])^2 * cell.h0[:, k] +
            2.0 * dt * cell.lorenz[3, k] * cell.h2[:, k] +
            (dt * cell.lorenz[3, k])^2 * cell.h0[:, k]
        @. cell.h2[:, k] += dt * cell.lorenz[3, k] * cell.h0[:, k]
        @. cell.h1[:, k] += dt * cell.lorenz[2, k] * cell.h0[:, k]
    end

    # source -> f^{n+1}
    tau = aap_hs_collision_time(
        cell.prim,
        KS.gas.mi,
        KS.gas.ni,
        KS.gas.me,
        KS.gas.ne,
        KS.gas.Kn[1],
    )

    # interspecies interaction
    prim = copy(cell.prim)
    prim = aap_hs_prim(cell.prim, tau, KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.Kn[1])
    g = mixture_maxwellian(KS.vSpace.u, prim)

    # BGK term
    Mu, Mv, Mw, MuL, MuR = mixture_gauss_moments(prim, KS.gas.K)
    for k in axes(cell.h0, 2)
        @. cell.h0[:, k] =
            (cell.h0[:, k] + dt / tau[k] * g[:, k]) / (1.0 + dt / tau[k])
        @. cell.h1[:, k] =
            (cell.h1[:, k] + dt / tau[k] * Mv[1, k] * g[:, k]) /
            (1.0 + dt / tau[k])
        @. cell.h2[:, k] =
            (cell.h2[:, k] + dt / tau[k] * Mw[1, k] * g[:, k]) /
            (1.0 + dt / tau[k])
        @. cell.h3[:, k] =
            (cell.h3[:, k] + dt / tau[k] * (Mv[2, k] + Mw[2, k]) * g[:, k]) /
            (1.0 + dt / tau[k])
    end
end

cd(@__DIR__)
ks, ctr, face, simTime = initialize("briowu.txt")
KS = ks

dt = timestep(ks, ctr, simTime)
nt = Int(floor(ks.set.maxTime / dt))+1
res = zeros(5, 2)

@showprogress for iter in 1:nt
    #dt = timestep(KS, ctr, simTime)
    reconstruct!(ks, ctr)

    evolve!(ks, ctr, face, dt; mode=:kfvs, isPlasma=true)
    
    #update!(ks, ctr, face, dt, res; coll=:bgk, bc=:extra, isMHD=true)
    @inbounds Threads.@threads for i in 1:KS.pSpace.nx
        step_tst!(
            KS,
            face[i],
            ctr[i],
            face[i+1],
            dt,
        )
    end
    for i in 1:2
        ctr[1-i].w .= ctr[1].w
        ctr[1-i].prim .= ctr[1].prim
        ctr[KS.pSpace.nx+i].w .= ctr[KS.pSpace.nx].w
        ctr[KS.pSpace.nx+i].prim .= ctr[KS.pSpace.nx].prim

        ctr[1-i].h0 .= ctr[1].h0
        ctr[1-i].h1 .= ctr[1].h1
        ctr[1-i].h2 .= ctr[1].h2
        ctr[1-i].h3 .= ctr[1].h3
        ctr[1-i].E .= ctr[1].E
        ctr[1-i].B .= ctr[1].B
        ctr[1-i].ϕ = deepcopy(ctr[1].ϕ)
        ctr[1-i].ψ = deepcopy(ctr[1].ψ)
        ctr[1-i].lorenz .= ctr[1].lorenz

        ctr[KS.pSpace.nx+i].h0 .= ctr[KS.pSpace.nx].h0
        ctr[KS.pSpace.nx+i].h1 .= ctr[KS.pSpace.nx].h1
        ctr[KS.pSpace.nx+i].h2 .= ctr[KS.pSpace.nx].h2
        ctr[KS.pSpace.nx+i].h3 .= ctr[KS.pSpace.nx].h3
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
plot(ks.pSpace.x[1:ks.pSpace.nx], sol[:,2,:])

@save "kfvs.jld2" ks ctr