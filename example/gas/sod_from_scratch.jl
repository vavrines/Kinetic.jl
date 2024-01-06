"""
This file is an illustrative script built from scratch.
The target is to simulate the Sod shock tube problem using the BGK equation.
The numerical flux function is the kinetic flux vector splitting method.
"""

using Plots
using Base: @kwdef

#--- Data Structure ---#
"""
Collection of all structures
"""
struct SolverSet{TS,TP,TV,TG}
    # setup
    set::TS
    # physical space
    ps::TP
    # velocity space
    vs::TV
    # gas property
    gas::TG
end

"""
Structure of setup
"""
struct Setup{T1,T2}
    cfl::T1
    maxTime::T2
end

"""
Structure of physical space
"""
struct PSpace1D{TR,TI,TA}
    x0::TR
    x1::TR
    nx::TI
    x::TA
    dx::TA
end

"""
Generate uniform mesh
"""
function PSpace1D(X0, X1, NX)
    δ = (X1 - X0) / NX
    x = Array{Float64}(undef, NX)
    dx = similar(x)
    for i in eachindex(x)
        x[i] = X0 + (i - 0.5) * δ
        dx[i] = δ
    end

    return PSpace1D(X0, X1, NX, x, dx)
end

"""
Structure of velocity space
"""
struct VSpace1D{TR,TI,TA,TB}
    u0::TR
    u1::TR
    nu::TI
    u::TA
    du::TA
    weights::TB
end

"""
Generate uniform collocation points and quadrature rules
"""
function VSpace1D(U0, U1, NU)
    δ = (U1 - U0) / NU
    u = Array{Float64}(undef, NU)
    du = similar(u)
    weights = similar(u)

    for i in eachindex(u)
        u[i] = U0 + (i - 0.5) * δ
        du[i] = δ
        weights[i] = δ
    end

    return VSpace1D(U0, U1, NU, u, du, weights)
end

"""
Structure of gas properties
"""
@kwdef mutable struct Gas{T1,T2,T3,T4,T5,T6,T7}
    Kn::T1
    K::T2 = 2.0
    γ::T3 = 5 / 3
    ω::T4 = 0.81
    αᵣ::T5 = 1.0
    ωᵣ::T6 = 0.5
    μᵣ::T7 = ref_vhs_vis(Kn, αᵣ, ωᵣ)
end

"""
Structure of cell-averaged values
"""
struct ControlVolume2F{T}
    w::T
    prim::T
    h::T
    b::T
end

"""
Structure of face fluxes
"""
struct Interface2F{T}
    fw::T
    fh::T
    fb::T
end

#--- Physics ---#
"""
Transform conservative -> primitive variables
"""
function conserve_prim(W, γ)
    prim = zero(W)

    prim[1] = W[1]
    prim[2] = W[2] / W[1]
    prim[3] = 0.5 * W[1] / (γ - 1.0) / (W[3] - 0.5 * W[2]^2 / W[1])

    return prim
end

"""
Transform primitive -> conservative variables
"""
function prim_conserve(prim, γ)
    W = zero(prim)

    W[1] = prim[1]
    W[2] = prim[1] * prim[2]
    W[3] = 0.5 * prim[1] / prim[3] / (γ - 1.0) + 0.5 * prim[1] * prim[2]^2

    return W
end

"""
Calculate reduced Maxwellian distributions
"""
function maxwellian(u, prim, K)
    ρ, U, λ = prim
    h = @. ρ * sqrt(λ / π) * exp(-λ * (u - U)^2)
    b = @. h * K / (2.0 * λ)
    return h, b
end

"""
Calculate viscosity coefficient in the reference state using hard-sphere model
"""
ref_vhs_vis(Kn, alpha, omega) =
    5.0 * (alpha + 1.0) * (alpha + 2.0) * √π /
    (4.0 * alpha * (5.0 - 2.0 * omega) * (7.0 - 2.0 * omega)) * Kn

"""
Calculate viscosity coefficient using variable hard-sphere model
"""
function vhs_collision_time(prim, μᵣ, ω)
    ρ, λ = prim[1], prim[end]
    return μᵣ * 2.0 * λ^(1.0 - ω) / ρ
end

"""
Calculate speed of sound
"""
sound_speed(prim, γ) = (0.5 * γ / prim[end])^0.5

#--- Flux ---#
"""
Calculate fluxes
"""
function evolve!(KS, ctr, face, dt;)
    @inbounds for i = 2:KS.ps.nx
        flux_kfvs!(
            face[i].fw,
            face[i].fh,
            face[i].fb,
            ctr[i-1].h,
            ctr[i-1].b,
            ctr[i].h,
            ctr[i].b,
            KS.vs.u,
            KS.vs.weights,
            dt,
        )
    end

    return nothing
end

"""
Heaviside step function
"""
heaviside(x) = ifelse(x >= 0, one(x), zero(x))

"""
Kinetic flux vector splitting method
"""
function flux_kfvs!(fw, fh, fb, hL, bL, hR, bR, u, ω, dt)
    # upwind reconstruction
    δ = heaviside.(u)
    h = @. hL * δ + hR * (1.0 - δ)
    b = @. bL * δ + bR * (1.0 - δ)

    # calculate fluxes
    @. fh = dt * u * h
    @. fb = dt * u * b

    fw[1] = dt * sum(ω .* u .* h)
    fw[2] = dt * sum(ω .* u .^ 2 .* h)
    fw[3] = dt * 0.5 * (sum(ω .* u .^ 3 .* h) + sum(ω .* u .* b))

    return nothing
end

#--- Update ---#
"""
Calculate time step
"""
function timestep(KS, ctr)
    tmax = 0.0

    @inbounds for i = 1:KS.ps.nx
        prim = ctr[i].prim
        sos = sound_speed(prim, KS.gas.γ)
        vmax = max(KS.vs.u1, abs(prim[2])) + sos
        tmax = max(tmax, vmax / KS.ps.dx[i])
    end

    dt = KS.set.cfl / tmax

    return dt
end

"""
Update macroscopic and mesoscopic variables
"""
function update!(KS, ctr, face, dt)
    @inbounds for i = 2:KS.ps.nx-1
        step!(
            ctr[i].w,
            ctr[i].prim,
            ctr[i].h,
            ctr[i].b,
            face[i].fw,
            face[i].fh,
            face[i].fb,
            face[i+1].fw,
            face[i+1].fh,
            face[i+1].fb,
            KS.vs.u,
            KS.gas.γ,
            KS.gas.μᵣ,
            KS.gas.ω,
            KS.ps.dx[i],
            dt,
        )
    end

    return nothing
end

"""
Implicit-explicit step function
"""
function step!(w, prim, h, b, fwL, fhL, fbL, fwR, fhR, fbR, u, γ, μᵣ, ω, dx, dt)
    # update W^{n+1}
    @. w += (fwL - fwR) / dx
    prim .= conserve_prim(w, γ)
    
    # calculate M^{n+1} and tau^{n+1}
    MH, MB = maxwellian(u, prim, 2)
    τ = vhs_collision_time(prim, μᵣ, ω)

    # update distribution function
    for i in eachindex(u)
        h[i] = (h[i] + (fhL[i] - fhR[i]) / dx + dt / τ * MH[i]) / (1.0 + dt / τ)
        b[i] = (b[i] + (fbL[i] - fbR[i]) / dx + dt / τ * MB[i]) / (1.0 + dt / τ)
    end

    return nothing
end

#--- IO ---#
"""
Calculate initial condition of Sod shock tube
"""
function ic_sod(x)
    if x < 0.5
        return [1.0, 0.0, 0.5]
    else
        return [0.125, 0.0, 0.625]
    end
end

"""
Initialize cell and face variables
"""
function init_fvm(KS)
    γ = KS.gas.γ
    K = KS.gas.K
    ps = KS.ps
    vs = KS.vs

    ctr = Array{ControlVolume2F}(undef, KS.ps.nx)
    face = Array{Interface2F}(undef, KS.ps.nx + 1)

    for i in eachindex(ctr)
        prim = ic_sod(ps.x[i])
        w = prim_conserve(prim, γ)
        h, b = maxwellian(vs.u, prim, K)
        ctr[i] = ControlVolume2F(w, prim, h, b)
    end

    for i in eachindex(face)
        fw = zeros(3)
        fh = zeros(vs.nu)
        fb = zeros(vs.nu)
        face[i] = Interface2F(fw, fh, fb)
    end

    return ctr, face
end

"""
Extract solution for visualization
"""
function extract_sol(ps, ctr)
    sol = zeros(ps.nx, axes(ctr[1].prim)...)
    for i in axes(sol, 1)
        sol[i, :] .= ctr[i].prim
    end

    return sol
end

#--- Solution ---#
set = Setup(0.5, 0.2)
ps = PSpace1D(0, 1, 100)
vs = VSpace1D(-5.0, 5.0, 72)
gas = Gas(Kn = 1e-4)
ks = SolverSet(set, ps, vs, gas)
ctr, face = init_fvm(ks)
dt = timestep(ks, ctr)
nt = ks.set.maxTime ÷ dt |> Int
# main loop
for iter = 1:nt
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt)
end

#--- Visualization ---#
sol = extract_sol(ps, ctr)
plot(ks.ps.x, sol)
