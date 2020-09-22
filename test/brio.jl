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

function flux_tst!(
    fw::AbstractArray{<:AbstractFloat,2},
    fh0::AbstractArray{<:AbstractFloat,2},
    fh1::AbstractArray{<:AbstractFloat,2},
    fh2::AbstractArray{<:AbstractFloat,2},
    fh3::AbstractArray{<:AbstractFloat,2},
    wL::AbstractArray{<:Real,2},
    h0L::AbstractArray{<:AbstractFloat,2},
    h1L::AbstractArray{<:AbstractFloat,2},
    h2L::AbstractArray{<:AbstractFloat,2},
    h3L::AbstractArray{<:AbstractFloat,2},
    wR::AbstractArray{<:Real,2},
    h0R::AbstractArray{<:AbstractFloat,2},
    h1R::AbstractArray{<:AbstractFloat,2},
    h2R::AbstractArray{<:AbstractFloat,2},
    h3R::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
    inK::Real,
    γ::Real,
    mi::Real,
    ni::Real,
    me::Real,
    ne::Real,
    Kn::Real,
    dt::Real,
)

    #--- upwind reconstruction ---#
    δ = heaviside.(u)

    h0 = @. h0L * δ + h0R * (1.0 - δ)
    h1 = @. h1L * δ + h1R * (1.0 - δ)
    h2 = @. h2L * δ + h2R * (1.0 - δ)
    h3 = @. h3L * δ + h3R * (1.0 - δ)

    primL = mixture_conserve_prim(wL, γ)
    primR = mixture_conserve_prim(wR, γ)

    #--- construct interface distribution ---#
    Mu1, Mv1, Mw1, MuL1, MuR1 = mixture_gauss_moments(primL, inK)
    Muv1 = mixture_moments_conserve(MuL1, Mv1, Mw1, 0, 0, 0)
    Mu2, Mv2, Mw2, MuL2, MuR2 = mixture_gauss_moments(primR, inK)
    Muv2 = mixture_moments_conserve(MuR2, Mv2, Mw2, 0, 0, 0)

    w = similar(wL)
    for j in axes(w, 2)
        @. w[:, j] = primL[1, j] * Muv1[:, j] + primR[1, j] * Muv2[:, j]
    end
    prim = mixture_conserve_prim(w, γ)

    tau = aap_hs_collision_time(prim, mi, ni, me, ne, Kn)
    @. tau +=
        abs(primL[1, :] / primL[end, :] - primR[1, :] / primR[end, :]) /
        (primL[1, :] / primL[end, :] + primR[1, :] / primR[end, :]) *
        dt *
        5.0
    #prim = aap_hs_prim(prim, tau, mi, ni, me, ne, Kn)

    Mt = zeros(2, 2)
    @. Mt[2, :] = tau * (1.0 - exp(-dt / tau)) # f0
    @. Mt[1, :] = dt - Mt[2, :] # M0

    #--- calculate fluxes ---#
    Mu, Mv, Mw, MuL, MuR = mixture_gauss_moments(prim, inK)
    Muv = mixture_moments_conserve(Mu, Mv, Mw, 1, 0, 0)

    # flux from M0
    for j in axes(fw, 2)
        @. fw[:, j] = Mt[1, j] * prim[1, j] * Muv[:, j]
    end

    # flux from f0
    g0 = mixture_maxwellian(u, prim)

    g1 = similar(h0)
    g2 = similar(h0)
    g3 = similar(h0)
    for j in axes(g0, 2)
        g1[:, j] .= Mv[1, j] .* g0[:, j]
        g2[:, j] .= Mw[1, j] .* g0[:, j]
        g3[:, j] .= (Mv[2, j] + Mw[2, j]) .* g0[:, j]
    end

    for j in axes(fw, 2)
        fw[1, j] += Mt[2, j] * sum(ω[:, j] .* u[:, j] .* h0[:, j])
        fw[2, j] += Mt[2, j] * sum(ω[:, j] .* u[:, j] .^ 2 .* h0[:, j])
        fw[3, j] += Mt[2, j] * sum(ω[:, j] .* u[:, j] .* h1[:, j])
        fw[4, j] += Mt[2, j] * sum(ω[:, j] .* u[:, j] .* h2[:, j])
        fw[5, j] +=
            Mt[2, j] *
            0.5 *
            (
                sum(ω[:, j] .* u[:, j] .^ 3 .* h0[:, j]) +
                sum(ω[:, j] .* u[:, j] .* h3[:, j])
            )

        @. fh0[:, j] = Mt[1, j] * u[:, j] * g0[:, j] + Mt[2, j] * u[:, j] * h0[:, j]
        @. fh1[:, j] = Mt[1, j] * u[:, j] * g1[:, j] + Mt[2, j] * u[:, j] * h1[:, j]
        @. fh2[:, j] = Mt[1, j] * u[:, j] * g2[:, j] + Mt[2, j] * u[:, j] * h2[:, j]
        @. fh3[:, j] = Mt[1, j] * u[:, j] * g3[:, j] + Mt[2, j] * u[:, j] * h3[:, j]
    end

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
#=
    # temperature protection
    if cell.prim[5, 1] < 0
        @warn ("ion temperature update is negative")
        cell.w .= w_old
        cell.prim .= prim_old
    elseif cell.prim[5, 2] < 0
        @warn ("electron temperature update is negative")
        cell.w .= w_old
        cell.prim .= prim_old
    end

    # source -> w^{n+1}
    if isMHD == false
        #=
        # DifferentialEquations.jl
        tau = get_tau(cell.prim, KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.Kn[1])
        for j in axes(wRan, 2)
        prob = ODEProblem( mixture_source,
            vcat(cell.w[1:5,j,1], cell.w[1:5,j,2]),
            dt,
            (tau[1], tau[2], KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.Kn[1], KS.gas.γ) )
        sol = solve(prob, Rosenbrock23())

        cell.w[1:5,j,1] .= sol[end][1:5]
        cell.w[1:5,j,2] .= sol[end][6:10]
        for k=1:2
        cell.prim[:,j,k] .= Kinetic.conserve_prim(cell.w[:,j,k], KS.gas.γ)
        end
        end
        =#
        
        # explicit
        tau = aap_hs_collision_time(cell.prim, KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.Kn[1])
        mprim = aap_hs_prim(cell.prim, tau, KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.Kn[1])
        mw = mixture_prim_conserve(mprim, KS.gas.γ)
        for k=1:2
            @. cell.w[:,k] += (mw[:,k] - w_old[:,k]) * dt / tau[k]
        end
        cell.prim .= mixture_conserve_prim(cell.w, KS.gas.γ)
    end
=#
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
    prim = deepcopy(cell.prim)
    prim = aap_hs_prim(prim, tau, KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.Kn[1])
    g = mixture_maxwellian(KS.vSpace.u, prim)

    # BGK term
    Mu, Mv, Mw, MuL, MuR = mixture_gauss_moments(prim, KS.gas.K)
    for k in axes(cell.h0, 2)
        @. cell.h0[:, k] =
            (cell.h0[:, k] + dt / tau[k] * g[:, k]) / (1.0 + dt / tau[k])
        @. cell.h1[:, k] =
            (cell.h1[:, k] + dt / tau[k] * Mv[1, k] * g[:, k]) / (1.0 + dt / tau[k])
        @. cell.h2[:, k] =
            (cell.h2[:, k] + dt / tau[k] * Mw[1, k] * g[:, k]) / (1.0 + dt / tau[k])
        @. cell.h3[:, k] =
            (cell.h3[:, k] + dt / tau[k] * (Mv[2, k] + Mw[2, k]) * g[:, k]) / (1.0 + dt / tau[k])
    end

end

cd(@__DIR__)
ks, ctr, face, simTime = initialize("briowu_1d.txt")
KS = ks

dt = timestep(ks, ctr, simTime)
nt = Int(floor(ks.set.maxTime / dt))+1
res = zeros(5, 2)

@showprogress for iter in 1:nt
    #dt = timestep(KS, ctr, simTime)
    reconstruct!(ks, ctr)

    evolve!(ks, ctr, face, dt; mode=:kfvs, isPlasma=true)
    #=
    @inbounds Threads.@threads for i = 1:KS.pSpace.nx+1
        #flux_kfvs!(
        flux_tst!(
            face[i].fw,
            face[i].fh0,
            face[i].fh1,
            face[i].fh2,
            face[i].fh3,
            ctr[i-1].w .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sw,
            ctr[i-1].h0 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh0,
            ctr[i-1].h1 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh1,
            ctr[i-1].h2 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh2,
            ctr[i-1].h3 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh3,
            ctr[i].w .- 0.5 .* ctr[i].dx .* ctr[i].sw,
            ctr[i].h0 .- 0.5 .* ctr[i].dx .* ctr[i].sh0,
            ctr[i].h1 .- 0.5 .* ctr[i].dx .* ctr[i].sh1,
            ctr[i].h2 .- 0.5 .* ctr[i].dx .* ctr[i].sh2,
            ctr[i].h3 .- 0.5 .* ctr[i].dx .* ctr[i].sh3,
            KS.vSpace.u,
            KS.vSpace.weights,
            ks.gas.K,
            ks.gas.γ,
            ks.gas.mi,
            ks.gas.ni,
            ks.gas.me,
            ks.gas.ne,
            ks.gas.Kn[1],
            dt,
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
    #update!(ks, ctr, face, dt, res; coll=:bgk, bc=:extra, isMHD=true)
    
    @inbounds Threads.@threads for i in 2:KS.pSpace.nx-1
        step_tst!(
            KS,
            face[i],
            ctr[i],
            face[i+1],
            dt,
        )
    end
    update_boundary!(
        KS, 
        ctr, 
        face, 
        dt, 
        res; 
        coll=:bgk, 
        bc=:extra, 
        isMHD=true,
    )
    
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
plot(ks.pSpace.x[1:ks.pSpace.nx], sol[:,2,:])
plot(ks.pSpace.x[1:ks.pSpace.nx], sol[:,6,:])