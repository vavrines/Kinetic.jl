using Revise, Kinetic
#=
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
=#
function flux_mhd!(
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

cd(@__DIR__)
ks, ctr, face, simTime = Kinetic.initialize("briowu.txt")
KS = ks

dt = Kinetic.timestep(ks, ctr, simTime)
nt = Int(floor(ks.set.maxTime / dt))+1
res = zeros(5, 2)

for iter in 1:nt
    #dt = timestep(KS, ctr, simTime)
    Kinetic.reconstruct!(ks, ctr)
    Kinetic.evolve!(ks, ctr, face, dt; mode=:kfvs, isPlasma=true)
    #=    
    @inbounds Threads.@threads for i = 1:KS.pSpace.nx+1
        #=
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
        =#
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
    Kinetic.update!(ks, ctr, face, dt, res; coll=:bgk, bc=:extra, isMHD=true)
end
#=
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
=#