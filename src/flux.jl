# ============================================================
# Flux functions for finite volume method
# ============================================================


export flux_gks!
export flux_kfvs!
export flux_kcu!
export flux_ugks!
export flux_boundary_maxwell!
export flux_em!


"""
Gas kinetic Navier-Stokes flux

> @param[in] : conservative variables and their slopes at left/right sides of interface
> @param[in] : thermodynamic and molecular parameters
> @param[in] : time step and cell size

< @return : flux of conservative variables

"""

function flux_gks!(
    fw::AbstractArray{<:AbstractFloat,1},
    wL::AbstractArray{<:AbstractFloat,1},
    wR::AbstractArray{<:AbstractFloat,1},
    γ::Real,
    inK::Real,
    μᵣ::Real,
    ω::Real,
    dt::Real,
    dx::Real,
    swL = zeros(eltype(fw), axes(wL))::AbstractArray{<:AbstractFloat,1},
    swR = zeros(eltype(fw), axes(wR))::AbstractArray{<:AbstractFloat,1},
)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    Mu1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Mu2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)

    w =
        primL[1] .* moments_conserve(MuL1, Mxi1, 0, 0) .+
        primR[1] .* moments_conserve(MuR2, Mxi2, 0, 0)
    prim = conserve_prim(w, γ)
    tau =
        vhs_collision_time(prim, μᵣ, ω) +
        2.0 * dt * abs(primL[1] / primL[end] - primR[1] / primR[end]) /
        (primL[1] / primL[end] + primR[1] / primR[end])

    faL = pdf_slope(primL, swL, inK)
    sw = -primL[1] .* moments_conserve_slope(faL, Mu1, Mxi1, 1)
    faTL = pdf_slope(primL, sw, inK)

    faR = pdf_slope(primR, swR, inK)
    sw = -primR[1] .* moments_conserve_slope(faR, Mu2, Mxi2, 1)
    faTR = pdf_slope(primR, sw, inK)

    Mu, Mxi, MuL, MuR = gauss_moments(prim, inK)
    swL = (w .- wL) ./ (0.5 * dx)
    swR = (wR .- w) ./ (0.5 * dx)
    gaL = pdf_slope(prim, swL, inK)
    gaR = pdf_slope(prim, swR, inK)
    sw =
        -prim[1] .* (
            moments_conserve_slope(gaL, MuL, Mxi, 1) .+
            moments_conserve_slope(gaR, MuR, Mxi, 1)
        )
    # ga = pdf_slope(prim, sw, inK)
    # sw = -prim[1] .* moments_conserve_slope(ga, Mu, Mxi, 1)
    gaT = pdf_slope(prim, sw, inK)

    # time-integration constants
    Mt = zeros(5)
    Mt[4] = tau * (1.0 - exp(-dt / tau))
    Mt[5] = -tau * dt * exp(-dt / tau) + tau * Mt[4]
    Mt[1] = dt - Mt[4]
    Mt[2] = -tau * Mt[1] + Mt[5]
    Mt[3] = 0.5 * dt^2 - tau * Mt[1]

    # flux related to central distribution
    Muv = moments_conserve(Mu, Mxi, 1, 0)
    MauL = moments_conserve_slope(gaL, MuL, Mxi, 2)
    MauR = moments_conserve_slope(gaR, MuR, Mxi, 2)
    # Mau = moments_conserve_slope(ga, MuR, Mxi, 2)
    MauT = moments_conserve_slope(gaT, Mu, Mxi, 1)

    fw .=
        Mt[1] .* prim[1] .* Muv .+ Mt[2] .* prim[1] .* (MauL .+ MauR) .+
        Mt[3] .* prim[1] .* MauT
    # fw .= Mt[1] .* prim[1] .* Muv .+ Mt[2] .* prim[1] .* Mau .+ Mt[3] .* prim[1] .* MauT

    # flux related to upwind distribution
    MuvL = moments_conserve(MuL1, Mxi1, 1, 0)
    MauL = moments_conserve_slope(faL, MuL1, Mxi1, 2)
    MauLT = moments_conserve_slope(faTL, MuL1, Mxi1, 1)

    MuvR = moments_conserve(MuR2, Mxi2, 1, 0)
    MauR = moments_conserve_slope(faR, MuR2, Mxi2, 2)
    MauRT = moments_conserve_slope(faTR, MuR2, Mxi2, 1)

    @. fw +=
        Mt[4] * primL[1] * MuvL - (Mt[5] + tau * Mt[4]) * primL[1] * MauL -
        tau * Mt[4] * primL[1] * MauLT + Mt[4] * primR[1] * MuvR -
        (Mt[5] + tau * Mt[4]) * primR[1] * MauR -
        tau * Mt[4] * primR[1] * MauRT
    # @. fw += Mt[4] * primL[1] * MuvL + Mt[4] * primR[1] * MuvR

    return fw

end


function flux_gks!(
    fw::AbstractArray{<:AbstractFloat,1},
    wL::AbstractArray{<:AbstractFloat,1},
    wR::AbstractArray{<:AbstractFloat,1},
    γ::Real,
    inK::Real,
    μᵣ::Real,
    ω::Real,
    dt::Real,
    dx::Real,
    dy::Real,
    swL = zeros(eltype(fw), axes(wL))::AbstractArray{<:AbstractFloat,1},
    swR = zeros(eltype(fw), axes(wR))::AbstractArray{<:AbstractFloat,1},
)

    Mu1, Mv1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Mu2, Mv2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)

    w =
        primL[1] .* moments_conserve(MuL1, Mv1, Mxi1, 0, 0, 0) .+
        primR[1] .* moments_conserve(MuR2, Mv2, Mxi2, 0, 0, 0)
    prim = conserve_prim(w, γ)
    tau =
        vhs_collision_time(prim, μᵣ, ω) +
        2.0 * dt * abs(primL[1] / primL[end] - primR[1] / primR[end]) /
        (primL[1] / primL[end] + primR[1] / primR[end])

    faL = pdf_slope(primL, swL, inK)
    sw = -primL[1] .* moments_conserve_slope(faL, Mu1, Mv1, Mxi1, 1, 0)
    faTL = pdf_slope(primL, sw, inK)

    faR = pdf_slope(primR, swR, inK)
    sw = -primR[1] .* moments_conserve_slope(faR, Mu2, Mv1, Mxi2, 1, 0)
    faTR = pdf_slope(primR, sw, inK)

    Mu, Mv, Mxi, MuL, MuR = gauss_moments(prim, inK)
    swL = (w .- wL) ./ (0.5 * dx)
    swR = (wR .- w) ./ (0.5 * dx)
    gaL = pdf_slope(prim, swL, inK)
    gaR = pdf_slope(prim, swR, inK)
    sw =
        -prim[1] .* (
            moments_conserve_slope(gaL, MuL, Mv, Mxi, 1, 0) .+
            moments_conserve_slope(gaR, MuR, Mv, Mxi, 1, 0)
        )
    # ga = pdf_slope(prim, sw, inK)
    # sw = -prim[1] .* moments_conserve_slope(ga, Mu, Mv, Mxi, 1, 0)
    gaT = pdf_slope(prim, sw, inK)

    # time-integration constants
    Mt = zeros(5)
    Mt[4] = tau * (1.0 - exp(-dt / tau))
    Mt[5] = -tau * dt * exp(-dt / tau) + tau * Mt[4]
    Mt[1] = dt - Mt[4]
    Mt[2] = -tau * Mt[1] + Mt[5]
    Mt[3] = 0.5 * dt^2 - tau * Mt[1]

    # flux related to central distribution
    Muv = moments_conserve(Mu, Mv, Mxi, 1, 0, 0)
    MauL = moments_conserve_slope(gaL, MuL, Mv, Mxi, 2, 0)
    MauR = moments_conserve_slope(gaR, MuR, Mv, Mxi, 2, 0)
    # Mau = moments_conserve_slope(ga, Mu, Mv, Mxi, 2, 0)
    MauT = moments_conserve_slope(gaT, Mu, Mv, Mxi, 1, 0)

    fw .=
        Mt[1] .* prim[1] .* Muv .+ Mt[2] .* prim[1] .* (MauL .+ MauR) .+
        Mt[3] .* prim[1] .* MauT
    # fw .= Mt[1] .* prim[1] .* Muv .+ Mt[2] .* prim[1] .* Mau .+ Mt[3] .* prim[1] .* MauT

    # flux related to upwind distribution
    MuvL = moments_conserve(MuL1, Mv1, Mxi1, 1, 0, 0)
    MauL = moments_conserve_slope(faL, MuL1, Mv1, Mxi1, 2, 0)
    MauLT = moments_conserve_slope(faTL, MuL1, Mv1, Mxi1, 1, 0)

    MuvR = moments_conserve(MuR2, Mv2, Mxi2, 1, 0, 0)
    MauR = moments_conserve_slope(faR, MuR2, Mv2, Mxi2, 2, 0)
    MauRT = moments_conserve_slope(faTR, MuR2, Mv2, Mxi2, 1, 0)

    @. fw +=
        Mt[4] * primL[1] * MuvL - (Mt[5] + tau * Mt[4]) * primL[1] * MauL -
        tau * Mt[4] * primL[1] * MauLT + Mt[4] * primR[1] * MuvR -
        (Mt[5] + tau * Mt[4]) * primR[1] * MauR -
        tau * Mt[4] * primR[1] * MauRT
    # @. fw += Mt[4] * primL[1] * MuvL + Mt[4] * primR[1] * MuvR

    return fw .* dy

end


"""
Kinetic flux vector splitting (KFVS) flux

> @param[in] : particle distribution functions and their slopes at left/right sides of interface
> @param[in] : particle velocity quadrature points and weights
> @param[in] : time step and cell size

< @return : flux of particle distribution function and its velocity moments on conservative variables

"""

# ------------------------------------------------------------
# 1D1F1V flux
# ------------------------------------------------------------
function flux_kfvs!(
    fw::AbstractArray{<:AbstractFloat,1},
    ff::AbstractArray{<:AbstractFloat,1},
    fL::AbstractArray{<:AbstractFloat,1},
    fR::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    ω::AbstractArray{<:AbstractFloat,1},
    dt::AbstractFloat,
    sfL = zeros(eltype(fL), axes(fL))::AbstractArray{<:AbstractFloat,1},
    sfR = zeros(eltype(fR), axes(fR))::AbstractArray{<:AbstractFloat,1},
)

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    f = @. fL * δ + fR * (1.0 - δ)
    sf = @. sfL * δ + sfR * (1.0 - δ)

    # --- calculate fluxes ---#
    fw[1] = dt * sum(ω .* u .* f) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sf)
    fw[2] = dt * sum(ω .* u .^ 2 .* f) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sf)
    fw[3] =
        dt * 0.5 * sum(ω .* u .^ 3 .* f) -
        0.5 * dt^2 * 0.5 * sum(ω .* u .^ 4 .* sf)

    @. ff = dt * u * f - 0.5 * dt^2 * u^2 * sf

end


# ------------------------------------------------------------
# 1D1F3V flux
# ------------------------------------------------------------
function flux_kfvs!(
    fw::AbstractArray{<:AbstractFloat,1},
    ff::AbstractArray{<:AbstractFloat,3},
    fL::AbstractArray{<:AbstractFloat,3},
    fR::AbstractArray{<:AbstractFloat,3},
    u::AbstractArray{<:AbstractFloat,3},
    v::AbstractArray{<:AbstractFloat,3},
    w::AbstractArray{<:AbstractFloat,3},
    ω::AbstractArray{<:AbstractFloat,3},
    dt::AbstractFloat,
    sfL = zeros(eltype(fL), axes(fL))::AbstractArray{<:AbstractFloat,3},
    sfR = zeros(eltype(fR), axes(fR))::AbstractArray{<:AbstractFloat,3},
)

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    f = @. fL * δ + fR * (1.0 - δ)
    sf = @. sfL * δ + sfR * (1.0 - δ)

    # --- calculate fluxes ---#
    fw[1] = dt * sum(ω .* u .* f) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sf)
    fw[2] = dt * sum(ω .* u .^ 2 .* f) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sf)
    fw[3] =
        dt * sum(ω .* u .* v .* f) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* v .* sf)
    fw[4] =
        dt * sum(ω .* u .* w .* f) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* w .* sf)
    fw[5] =
        dt * 0.5 * sum(ω .* u .* (u .^ 2 .+ v .^ 2 .+ w .^ 2) .* f) -
        0.5 *
        dt^2 *
        0.5 *
        sum(ω .* u .^ 2 .* (u .^ 2 .+ v .^ 2 .+ w .^ 2) .* sf)

    @. ff = dt * u * f - 0.5 * dt^2 * u^2 * sf

end


# ------------------------------------------------------------
# 1D2F flux
# ------------------------------------------------------------
function flux_kfvs!(
    fw::AbstractArray{<:AbstractFloat,1},
    fh::AbstractArray{<:AbstractFloat,1},
    fb::AbstractArray{<:AbstractFloat,1},
    hL::AbstractArray{<:AbstractFloat,1},
    bL::AbstractArray{<:AbstractFloat,1},
    hR::AbstractArray{<:AbstractFloat,1},
    bR::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    ω::AbstractArray{<:AbstractFloat,1},
    dt::AbstractFloat,
    shL = zeros(eltype(hL), axes(hL))::AbstractArray{<:AbstractFloat,1},
    sbL = zeros(eltype(bL), axes(bL))::AbstractArray{<:AbstractFloat,1},
    shR = zeros(eltype(hR), axes(hR))::AbstractArray{<:AbstractFloat,1},
    sbR = zeros(eltype(bR), axes(bR))::AbstractArray{<:AbstractFloat,1},
)

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    h = @. hL * δ + hR * (1.0 - δ)
    b = @. bL * δ + bR * (1.0 - δ)

    sh = @. shL * δ + shR * (1.0 - δ)
    sb = @. sbL * δ + sbR * (1.0 - δ)

    # --- calculate fluxes ---#
    fw[1] = dt * sum(ω .* u .* h) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh)
    fw[2] = dt * sum(ω .* u .^ 2 .* h) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sh)
    fw[3] =
        dt * 0.5 * (sum(ω .* u .^ 3 .* h) + sum(ω .* u .* b)) -
        0.5 * dt^2 * 0.5 * (sum(ω .* u .^ 4 .* sh) + sum(ω .* u .^ 2 .* sb))

    @. fh = dt * u * h - 0.5 * dt^2 * u^2 * sh
    @. fb = dt * u * b - 0.5 * dt^2 * u^2 * sb

end


# ------------------------------------------------------------
# 1D4F flux
# ------------------------------------------------------------
function flux_kfvs!(
    fw::AbstractArray{<:AbstractFloat,1},
    fh0::AbstractArray{<:AbstractFloat,1},
    fh1::AbstractArray{<:AbstractFloat,1},
    fh2::AbstractArray{<:AbstractFloat,1},
    fh3::AbstractArray{<:AbstractFloat,1},
    h0L::AbstractArray{<:AbstractFloat,1},
    h1L::AbstractArray{<:AbstractFloat,1},
    h2L::AbstractArray{<:AbstractFloat,1},
    h3L::AbstractArray{<:AbstractFloat,1},
    h0R::AbstractArray{<:AbstractFloat,1},
    h1R::AbstractArray{<:AbstractFloat,1},
    h2R::AbstractArray{<:AbstractFloat,1},
    h3R::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    ω::AbstractArray{<:AbstractFloat,1},
    dt::AbstractFloat,
    sh0L = zeros(eltype(h0L), axes(h0L))::AbstractArray{<:AbstractFloat,1},
    sh1L = zeros(eltype(h1L), axes(h1L))::AbstractArray{<:AbstractFloat,1},
    sh2L = zeros(eltype(h2L), axes(h2L))::AbstractArray{<:AbstractFloat,1},
    sh3L = zeros(eltype(h3L), axes(h3L))::AbstractArray{<:AbstractFloat,1},
    sh0R = zeros(eltype(h0R), axes(h0R))::AbstractArray{<:AbstractFloat,1},
    sh1R = zeros(eltype(h1R), axes(h1R))::AbstractArray{<:AbstractFloat,1},
    sh2R = zeros(eltype(h2R), axes(h2R))::AbstractArray{<:AbstractFloat,1},
    sh3R = zeros(eltype(h3R), axes(h3R))::AbstractArray{<:AbstractFloat,1},
)

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    h0 = @. h0L * δ + h0R * (1.0 - δ)
    h1 = @. h1L * δ + h1R * (1.0 - δ)
    h2 = @. h2L * δ + h2R * (1.0 - δ)
    h3 = @. h3L * δ + h3R * (1.0 - δ)

    sh0 = @. sh0L * δ + sh0R * (1.0 - δ)
    sh1 = @. sh1L * δ + sh1R * (1.0 - δ)
    sh2 = @. sh2L * δ + sh2R * (1.0 - δ)
    sh3 = @. sh3L * δ + sh3R * (1.0 - δ)

    # --- calculate fluxes ---#
    fw[1] = dt * sum(ω .* u .* h0) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh0)
    fw[2] = dt * sum(ω .* u .^ 2 .* h0) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sh0)
    fw[3] = dt * sum(ω .* u .* h1) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh1)
    fw[4] = dt * sum(ω .* u .* h2) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh2)
    fw[5] =
        dt * 0.5 * (sum(ω .* u .^ 3 .* h0) + sum(ω .* u .* h3)) -
        0.5 * dt^2 * 0.5 * (sum(ω .* u .^ 4 .* sh0) + sum(ω .* u .^ 2 .* sh3))

    @. fh0 = dt * u * h0 - 0.5 * dt^2 * u^2 * sh0
    @. fh1 = dt * u * h1 - 0.5 * dt^2 * u^2 * sh1
    @. fh2 = dt * u * h2 - 0.5 * dt^2 * u^2 * sh2
    @. fh3 = dt * u * h3 - 0.5 * dt^2 * u^2 * sh3

end


# ------------------------------------------------------------
# 2D1F flux
# ------------------------------------------------------------
function flux_kfvs!(
    fw::AbstractArray{<:AbstractFloat,1},
    ff::AbstractArray{<:AbstractFloat,2},
    fL::AbstractArray{<:AbstractFloat,2},
    fR::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
    dt::Real,
    len::Real,
    sfL = zeros(eltype(fL), axes(fL))::AbstractArray{<:AbstractFloat,2},
    sfR = zeros(eltype(fR), axes(fR))::AbstractArray{<:AbstractFloat,2},
)

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    f = @. fL * δ + fR * (1.0 - δ)
    sf = @. sfL * δ + sfR * (1.0 - δ)

    # --- calculate fluxes ---#
    fw[1] = dt * sum(ω .* u .* f) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sf)
    fw[2] = dt * sum(ω .* u .^ 2 .* f) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sf)
    fw[3] =
        dt * sum(ω .* v .* u .* f) - 0.5 * dt^2 * sum(ω .* v .* u .^ 2 .* sf)
    fw[4] =
        dt * 0.5 * sum(ω .* u .* (u .^ 2 .+ v .^ 2) .* f) -
        0.5 * dt^2 * 0.5 * sum(ω .* u .^ 2 .* (u .^ 2 .+ v .^ 2) .* sf)
    fw .*= len

    @. ff = (dt * u * f - 0.5 * dt^2 * u^2 * sf) * len

end


# ------------------------------------------------------------
# 2D2F flux
# ------------------------------------------------------------
function flux_kfvs!(
    fw::AbstractArray{<:AbstractFloat,1},
    fh::AbstractArray{<:AbstractFloat,2},
    fb::AbstractArray{<:AbstractFloat,2},
    hL::AbstractArray{<:AbstractFloat,2},
    bL::AbstractArray{<:AbstractFloat,2},
    hR::AbstractArray{<:AbstractFloat,2},
    bR::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
    dt::Real,
    len::Real,
    shL = zeros(eltype(hL), axes(hL))::AbstractArray{<:AbstractFloat,2},
    sbL = zeros(eltype(bL), axes(bL))::AbstractArray{<:AbstractFloat,2},
    shR = zeros(eltype(hR), axes(hR))::AbstractArray{<:AbstractFloat,2},
    sbR = zeros(eltype(bR), axes(bR))::AbstractArray{<:AbstractFloat,2},
)

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    h = @. hL * δ + hR * (1.0 - δ)
    b = @. bL * δ + bR * (1.0 - δ)
    sh = @. shL * δ + shR * (1.0 - δ)
    sb = @. sbL * δ + sbR * (1.0 - δ)

    # --- calculate fluxes ---#
    fw[1] = dt * sum(ω .* u .* h) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh)
    fw[2] = dt * sum(ω .* u .^ 2 .* h) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sh)
    fw[3] =
        dt * sum(ω .* v .* u .* h) - 0.5 * dt^2 * sum(ω .* v .* u .^ 2 .* sh)
    fw[4] =
        dt * 0.5 * (sum(ω .* u .* (u .^ 2 .+ v .^ 2) .* h) + sum(ω .* u .* b)) -
        0.5 *
        dt^2 *
        0.5 *
        (sum(ω .* u .^ 2 .* (u .^ 2 .+ v .^ 2) .* sh) + sum(ω .* u .^ 2 .* sb))
    fw .*= len

    @. fh = (dt * u * h - 0.5 * dt^2 * u^2 * sh) * len
    @. fb = (dt * u * b - 0.5 * dt^2 * u^2 * sb) * len

end


"""
Kinetic central-upwind (KCU) method

> @param[in] : particle distribution functions and their slopes at left/right sides of interface
> @param[in] : particle velocity quadrature points and weights
> @param[in] : time step and cell size

< @return : flux of particle distribution function and its velocity moments on conservative variables

"""

# ------------------------------------------------------------
# 1D1F flux
# ------------------------------------------------------------
function flux_kcu!(
    fw::AbstractArray{<:Real,1},
    ff::AbstractArray{<:AbstractFloat,1},
    wL::AbstractArray{<:Real,1},
    fL::AbstractArray{<:AbstractFloat,1},
    wR::AbstractArray{<:Real,1},
    fR::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    ω::AbstractArray{<:AbstractFloat,1},
    inK::Real,
    γ::Real,
    visRef::Real,
    visIdx::Real,
    pr::Real,
    dt::Real,
)

    # --- upwind reconstruction ---#
    δ = heaviside.(u)
    f = @. fL * δ + fR * (1.0 - δ)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    # --- construct interface distribution ---#
    Mu1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Muv1 = moments_conserve(MuL1, Mxi1, 0, 0)
    Mu2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)
    Muv2 = moments_conserve(MuR2, Mxi2, 0, 0)

    w = similar(wL, 3)
    @. w = primL[1] * Muv1 + primR[1] * Muv2

    prim = conserve_prim(w, γ)
    tau = vhs_collision_time(prim, visRef, visIdx)
    tau +=
        abs(primL[1] / primL[end] - primR[1] / primR[end]) /
        (primL[1] / primL[end] + primR[1] / primR[end]) *
        dt *
        2.0

    Mt = zeros(2)
    Mt[2] = tau * (1.0 - exp(-dt / tau)) # f0
    Mt[1] = dt - Mt[2] # M0

    # --- calculate fluxes ---#
    Mu, Mxi, MuL, MuR = gauss_moments(prim, inK)

    # flux from M0
    Muv = moments_conserve(Mu, Mxi, 1, 0)
    @. fw = Mt[1] * prim[1] * Muv

    # flux from f0
    g = maxwellian(u, prim)

    fw[1] += Mt[2] * sum(ω .* u .* f)
    fw[2] += Mt[2] * sum(ω .* u .^ 2 .* f)
    fw[3] += Mt[2] * 0.5 * (sum(ω .* u .^ 3 .* f))

    @. ff = Mt[1] * u * g + Mt[2] * u * f

end


# ------------------------------------------------------------
# 1D2F flux
# ------------------------------------------------------------
function flux_kcu!(
    fw::AbstractArray{<:Real,1},
    fh::AbstractArray{<:AbstractFloat,1},
    fb::AbstractArray{<:AbstractFloat,1},
    wL::AbstractArray{<:Real,1},
    hL::AbstractArray{<:AbstractFloat,1},
    bL::AbstractArray{<:AbstractFloat,1},
    wR::AbstractArray{<:Real,1},
    hR::AbstractArray{<:AbstractFloat,1},
    bR::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    ω::AbstractArray{<:AbstractFloat,1},
    inK::Real,
    γ::Real,
    visRef::Real,
    visIdx::Real,
    pr::Real,
    dt::Real,
)

    # --- upwind reconstruction ---#
    δ = heaviside.(u)
    h = @. hL * δ + hR * (1.0 - δ)
    b = @. bL * δ + bR * (1.0 - δ)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    # --- construct interface distribution ---#
    Mu1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Muv1 = moments_conserve(MuL1, Mxi1, 0, 0)
    Mu2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)
    Muv2 = moments_conserve(MuR2, Mxi2, 0, 0)

    w = @. primL[1] * Muv1 + primR[1] * Muv2

    prim = conserve_prim(w, γ)
    tau = vhs_collision_time(prim, visRef, visIdx)
    tau +=
        abs(primL[1] / primL[end] - primR[1] / primR[end]) /
        (primL[1] / primL[end] + primR[1] / primR[end]) *
        dt *
        2.0

    Mt = zeros(2)
    Mt[2] = tau * (1.0 - exp(-dt / tau)) # f0
    Mt[1] = dt - Mt[2] # M0

    # --- calculate fluxes ---#
    Mu, Mxi, MuL, MuR = gauss_moments(prim, inK)

    # flux from M0
    Muv = moments_conserve(Mu, Mxi, 1, 0)
    @. fw = Mt[1] * prim[1] * Muv

    # flux from f0
    Mh = maxwellian(u, prim)
    Mb = Mh .* inK ./ (2.0 * prim[end])

    fw[1] += Mt[2] * sum(ω .* u .* h)
    fw[2] += Mt[2] * sum(ω .* u .^ 2 .* h)
    fw[3] += Mt[2] * 0.5 * (sum(ω .* u .^ 3 .* h) + sum(ω .* u .* b))

    @. fh = Mt[1] * u * Mh + Mt[2] * u * h
    @. fb = Mt[1] * u * Mb + Mt[2] * u * b

end


# ------------------------------------------------------------
# 2D1F flux
# ------------------------------------------------------------
function flux_kcu!(
    fw::AbstractArray{<:Real,1},
    ff::AbstractArray{<:AbstractFloat,2},
    wL::AbstractArray{<:Real,1},
    fL::AbstractArray{<:AbstractFloat,2},
    wR::AbstractArray{<:Real,1},
    fR::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
    inK::Real,
    γ::Real,
    visRef::Real,
    visIdx::Real,
    pr::Real,
    dt::Real,
    len::Real,
)

    # --- prepare ---#
    delta = heaviside.(u)

    # --- reconstruct initial distribution ---#
    δ = heaviside.(u)
    f = @. fL * δ + fR * (1.0 - δ)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    # --- construct interface distribution ---#
    Mu1, Mv1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Muv1 = moments_conserve(MuL1, Mv1, Mxi1, 0, 0, 0)
    Mu2, Mv2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)
    Muv2 = moments_conserve(MuR2, Mv2, Mxi2, 0, 0, 0)

    w = @. primL[1] * Muv1 + primR[1] * Muv2
    prim = conserve_prim(w, γ)
    tau = vhs_collision_time(prim, visRef, visIdx)
    tau +=
        abs(primL[1] / primL[end] - primR[1] / primR[end]) /
        (primL[1] / primL[end] + primR[1] / primR[end]) *
        dt *
        2.0

    Mt = zeros(2)
    Mt[2] = tau * (1.0 - exp(-dt / tau)) # f0
    Mt[1] = dt - Mt[2] # M0

    # --- calculate interface flux ---#
    Mu, Mv, Mxi, MuL, MuR = gauss_moments(prim, inK)

    # flux from M0
    Muv = moments_conserve(Mu, Mv, Mxi, 1, 0, 0)
    @. fw = Mt[1] * prim[1] * Muv * len

    # flux from f0
    g = maxwellian(u, v, prim)

    fw[1] += Mt[2] * sum(ω .* u .* f) * len
    fw[2] += Mt[2] * sum(ω .* u .^ 2 .* f) * len
    fw[3] += Mt[2] * sum(ω .* v .* u .* f) * len
    fw[4] += Mt[2] * 0.5 * (sum(ω .* u .* (u .^ 2 .+ v .^ 2) .* f)) * len

    @. ff = (Mt[1] * u * g + Mt[2] * u * f) * len

end


# ------------------------------------------------------------
# 2D2F flux
# ------------------------------------------------------------
function flux_kcu!(
    fw::AbstractArray{<:AbstractFloat,1},
    fh::AbstractArray{<:AbstractFloat,2},
    fb::AbstractArray{<:AbstractFloat,2},
    wL::AbstractArray{<:Real,1},
    hL::AbstractArray{<:AbstractFloat,2},
    bL::AbstractArray{<:AbstractFloat,2},
    wR::AbstractArray{<:Real,1},
    hR::AbstractArray{<:AbstractFloat,2},
    bR::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
    inK::Real,
    γ::Real,
    visRef::Real,
    visIdx::Real,
    pr::Real,
    dt::Real,
    len::Real,
)

    # --- prepare ---#
    delta = heaviside.(u)

    # --- reconstruct initial distribution ---#
    δ = heaviside.(u)
    h = @. hL * δ + hR * (1.0 - δ)
    b = @. bL * δ + bR * (1.0 - δ)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    # --- construct interface distribution ---#
    Mu1, Mv1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Muv1 = moments_conserve(MuL1, Mv1, Mxi1, 0, 0, 0)
    Mu2, Mv2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)
    Muv2 = moments_conserve(MuR2, Mv2, Mxi2, 0, 0, 0)

    w = @. primL[1] * Muv1 + primR[1] * Muv2
    prim = conserve_prim(w, γ)
    tau = vhs_collision_time(prim, visRef, visIdx)
    tau +=
        abs(primL[1] / primL[end] - primR[1] / primR[end]) /
        (primL[1] / primL[end] + primR[1] / primR[end]) *
        dt *
        2.0

    Mt = zeros(2)
    Mt[2] = tau * (1.0 - exp(-dt / tau)) # f0
    Mt[1] = dt - Mt[2] # M0

    # --- calculate interface flux ---#
    Mu, Mv, Mxi, MuL, MuR = gauss_moments(prim, inK)

    # flux from M0
    Muv = moments_conserve(Mu, Mv, Mxi, 1, 0, 0)
    @. fw = Mt[1] * prim[1] * Muv * len

    # flux from f0
    H = maxwellian(u, v, prim)
    B = H .* inK ./ (2.0 * prim[end])

    fw[1] += Mt[2] * sum(ω .* u .* h) * len
    fw[2] += Mt[2] * sum(ω .* u .^ 2 .* h) * len
    fw[3] += Mt[2] * sum(ω .* v .* u .* h) * len
    fw[4] +=
        Mt[2] *
        0.5 *
        (sum(ω .* u .* (u .^ 2 .+ v .^ 2) .* h) + sum(ω .* u .* b)) *
        len

    @. fh = (Mt[1] * u * H + Mt[2] * u * h) * len
    @. fb = (Mt[1] * u * B + Mt[2] * u * b) * len

end


"""
Kinetic central-upwind (KCU) method for multi-component gas

> @param[in] : particle distribution functions and their slopes at left/right sides of interface
> @param[in] : particle velocity quadrature points and weights
> @param[in] : time step

< @return : flux of particle distribution function and its velocity moments on conservative variables

"""

# ------------------------------------------------------------
# 1D1F flux with AAP model
# ------------------------------------------------------------
function flux_kcu!(
    fw::AbstractArray{<:AbstractFloat,2},
    ff::AbstractArray{<:AbstractFloat,2},
    wL::AbstractArray{<:Real,2},
    fL::AbstractArray{<:AbstractFloat,2},
    wR::AbstractArray{<:Real,2},
    fR::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
    inK::Real,
    γ::Real,
    mi::Real,
    ni::Real,
    me::Real,
    ne::Real,
    kn::Real,
    dt::Real,
)

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    f = @. fL * δ + fR * (1.0 - δ)

    primL = zeros(axes(wL))
    primR = similar(primL)
    for j = 1:2
        primL[:, j] .= conserve_prim(wL[:, j], γ)
        primR[:, j] .= conserve_prim(wR[:, j], γ)
    end

    # --- construct interface distribution ---#
    Mu1 = OffsetArray{Float64}(undef, 0:6, 1:2)
    Mxi1 = similar(Mu1)
    MuL1 = similar(Mu1)
    MuR1 = similar(Mu1)
    Mu2 = similar(Mu1)
    Mxi2 = similar(Mu1)
    MuL2 = similar(Mu1)
    MuR2 = similar(Mu1)
    Muv1 = similar(wL)
    Muv2 = similar(wL)
    for j = 1:2
        Mu1[:, j], Mxi1[:, j], MuL1[:, j], MuR1[:, j] =
            gauss_moments(primL[:, j], inK)
        Muv1[:, j] = moments_conserve(MuL1[:, j], Mxi1[:, j], 0, 0)
        Mu2[:, j], Mxi2[:, j], MuL2[:, j], MuR2[:, j] =
            gauss_moments(primR[:, j], inK)
        Muv2[:, j] = moments_conserve(MuR2[:, j], Mxi2[:, j], 0, 0)
    end

    w = similar(wL)
    prim = similar(wL)
    for j = 1:2
        @. w[:, j] = primL[1, j] * Muv1[:, j] + primR[1, j] * Muv2[:, j]
        prim[:, j] .= conserve_prim(w[:, j], γ)
    end

    tau = aap_hs_collision_time(prim, mi, ni, me, ne, kn)
    @. tau +=
        abs(
            cellL.prim[1, :] / cellL.prim[end, :] -
            cellR.prim[1, :] / cellR.prim[end, :],
        ) / (
            cellL.prim[1, :] / cellL.prim[end, :] +
            cellR.prim[1, :] / cellR.prim[end, :]
        ) *
        dt *
        2.0
    prim = aap_hs_prim(prim, tau, mi, ni, me, ne, kn)

    Mt = zeros(2, 2)
    @. Mt[2, :] = tau * (1.0 - exp(-dt / tau)) # f0
    @. Mt[1, :] = dt - Mt[2, :] # M0

    # --- calculate fluxes ---#
    Mu = similar(Mu1)
    Mxi = similar(Mu1)
    MuL = similar(Mu1)
    MuR = similar(Mu1)
    Muv = similar(wL)
    for j in axes(Mu1, 2)
        Mu[:, j], Mxi[:, j], MuL[:, j], MuR[:, j] =
            gauss_moments(prim[:, j], inK)
        Muv[:, j] .= moments_conserve(Mu[:, j], Mxi[:, j], 1, 0)
    end

    # flux from M0
    for j = 1:2
        @. fw[:, j] = Mt[1, j] * prim[1, j] * Muv[:, j]
    end

    # flux from f0
    g = similar(f)
    for j = 1:2
        g[:, j] .= maxwellian(u[:, j], prim[:, j])
    end

    for j = 1:2
        fw[1, j] += Mt[2, j] * sum(ω[:, j] .* u[:, j] .* f[:, j])
        fw[2, j] += Mt[2, j] * sum(ω[:, j] .* u[:, j] .^ 2 .* f[:, j])
        fw[3, j] += Mt[2, j] * 0.5 * sum(ω[:, j] .* u[:, j] .^ 3 .* f[:, j])

        @. ff[:, j] =
            Mt[1, j] * u[:, j] * g[:, j] + Mt[2, j] * u[:, j] * f[:, j]
    end

end


# ------------------------------------------------------------
# 1D2F flux with AAP model
# ------------------------------------------------------------
function flux_kcu!(
    fw::AbstractArray{<:AbstractFloat,2},
    fh::AbstractArray{<:AbstractFloat,2},
    fb::AbstractArray{<:AbstractFloat,2},
    wL::AbstractArray{<:Real,2},
    hL::AbstractArray{<:AbstractFloat,2},
    bL::AbstractArray{<:AbstractFloat,2},
    wR::AbstractArray{<:Real,2},
    hR::AbstractArray{<:AbstractFloat,2},
    bR::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
    inK::Real,
    γ::Real,
    mi::Real,
    ni::Real,
    me::Real,
    ne::Real,
    kn::Real,
    dt::Real,
)

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    h = @. hL * δ + hR * (1.0 - δ)
    b = @. bL * δ + bR * (1.0 - δ)

    primL = similar(wL)
    primR = similar(primL)
    for j = 1:2
        primL[:, j] .= conserve_prim(wL[:, j], γ)
        primR[:, j] .= conserve_prim(wR[:, j], γ)
    end

    # --- construct interface distribution ---#
    Mu1 = OffsetArray{Float64}(undef, 0:6, 1:2)
    Mxi1 = similar(Mu1)
    MuL1 = similar(Mu1)
    MuR1 = similar(Mu1)
    Mu2 = similar(Mu1)
    Mxi2 = similar(Mu1)
    MuL2 = similar(Mu1)
    MuR2 = similar(Mu1)
    Muv1 = similar(wL)
    Muv2 = similar(wL)
    for j = 1:2
        Mu1[:, j], Mxi1[:, j], MuL1[:, j], MuR1[:, j] =
            gauss_moments(primL[:, j], inK)
        Muv1[:, j] = moments_conserve(MuL1[:, j], Mxi1[:, j], 0, 0)
        Mu2[:, j], Mxi2[:, j], MuL2[:, j], MuR2[:, j] =
            gauss_moments(primR[:, j], inK)
        Muv2[:, j] = moments_conserve(MuR2[:, j], Mxi2[:, j], 0, 0)
    end

    w = similar(wL)
    prim = similar(wL)
    for j = 1:2
        @. w[:, j] = primL[1, j] * Muv1[:, j] + primR[1, j] * Muv2[:, j]
        prim[:, j] .= conserve_prim(w[:, j], γ)
    end

    tau = aap_hs_collision_time(prim, mi, ni, me, ne, kn)
    @. tau +=
        abs(
            cellL.prim[1, :] / cellL.prim[end, :] -
            cellR.prim[1, :] / cellR.prim[end, :],
        ) / (
            cellL.prim[1, :] / cellL.prim[end, :] +
            cellR.prim[1, :] / cellR.prim[end, :]
        ) *
        dt *
        2.0
    prim = aap_hs_prim(prim, tau, mi, ni, me, ne, kn)

    Mt = zeros(2, 2)
    @. Mt[2, :] = tau * (1.0 - exp(-dt / tau)) # f0
    @. Mt[1, :] = dt - Mt[2, :] # M0

    # --- calculate fluxes ---#
    Mu = similar(Mu1)
    Mxi = similar(Mu1)
    MuL = similar(Mu1)
    MuR = similar(Mu1)
    Muv = similar(wL)
    for j in axes(Mu1, 2)
        Mu[:, j], Mxi[:, j], MuL[:, j], MuR[:, j] =
            gauss_moments(prim[:, j], inK)
        Muv[:, j] .= moments_conserve(Mu[:, j], Mxi[:, j], 1, 0)
    end

    # flux from M0
    for j = 1:2
        @. fw[:, j] = Mt[1, j] * prim[1, j] * Muv[:, j]
    end

    # flux from f0
    g0 = similar(h)
    g1 = similar(b)
    for j = 1:2
        g0[:, j] .= maxwellian(u[:, j], prim[:, j])
        g1[:, j] .= g0[:, j] .* inK ./ (2.0 * prim[end, j])
    end

    for j = 1:2
        fw[1, j] += Mt[2, j] * sum(ω[:, j] .* u[:, j] .* h[:, j])
        fw[2, j] += Mt[2, j] * sum(ω[:, j] .* u[:, j] .^ 2 .* h[:, j])
        fw[3, j] +=
            Mt[2, j] *
            0.5 *
            (
                sum(ω[:, j] .* u[:, j] .^ 3 .* h[:, j]) +
                sum(ω[:, j] .* u[:, j] .* b[:, j])
            )

        @. fh[:, j] =
            Mt[1, j] * u[:, j] * g0[:, j] + Mt[2, j] * u[:, j] * h[:, j]
        @. fb[:, j] =
            Mt[1, j] * u[:, j] * g1[:, j] + Mt[2, j] * u[:, j] * b[:, j]
    end

end


# ------------------------------------------------------------
# 1D4F flux with AAP model
# ------------------------------------------------------------
function flux_kcu!(
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
    kn::Real,
    dt::Real,
)

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    h0 = @. h0L * δ + h0R * (1.0 - δ)
    h1 = @. h1L * δ + h1R * (1.0 - δ)
    h2 = @. h2L * δ + h2R * (1.0 - δ)
    h3 = @. h3L * δ + h3R * (1.0 - δ)

    primL = mixture_conserve_prim(wL, γ)
    primR = mixture_conserve_prim(wR, γ)

    # --- construct interface distribution ---#
    Mu1, Mv1, Mw1, MuL1, MuR1 = mixture_gauss_moments(primL, inK)
    Muv1 = mixture_moments_conserve(MuL1, Mv1, Mw1, 0, 0, 0)
    Mu2, Mv2, Mw2, MuL2, MuR2 = mixture_gauss_moments(primR, inK)
    Muv2 = mixture_moments_conserve(MuR2, Mv2, Mw2, 0, 0, 0)

    w = similar(wL)
    for j = 1:2
        @. w[:, j] = primL[1, j] * Muv1[:, j] + primR[1, j] * Muv2[:, j]
    end
    prim = mixture_conserve_prim(w, γ)

    tau = aap_hs_collision_time(prim, mi, ni, me, ne, kn)
    @. tau +=
        abs(primL[1, :] / primL[end, :] - primR[1, :] / primR[end, :]) /
        (primL[1, :] / primL[end, :] + primR[1, :] / primR[end, :]) *
        dt *
        5.0
    # prim = aap_hs_prim(prim, tau, mi, ni, me, ne, kn)

    Mt = zeros(2, 2)
    @. Mt[2, :] = tau * (1.0 - exp(-dt / tau)) # f0
    @. Mt[1, :] = dt - Mt[2, :] # M0

    # --- calculate fluxes ---#
    Mu, Mv, Mw, MuL, MuR = mixture_gauss_moments(prim, inK)
    Muv = mixture_moments_conserve(Mu, Mv, Mw, 1, 0, 0)

    # flux from M0
    for j = 1:2
        @. fw[:, j] = Mt[1, j] * prim[1, j] * Muv[:, j]
    end

    # flux from f0
    g0 = mixture_maxwellian(u, prim)

    g1 = similar(h0)
    g2 = similar(h0)
    g3 = similar(h0)
    for j = 1:2
        g1[:, j] .= Mv[1, j] .* g0[:, j]
        g2[:, j] .= Mw[1, j] .* g0[:, j]
        g3[:, j] .= (Mv[2, j] + Mw[2, j]) .* g0[:, j]
    end

    for j = 1:2
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

        @. fh0[:, j] =
            Mt[1, j] * u[:, j] * g0[:, j] + Mt[2, j] * u[:, j] * h0[:, j]
        @. fh1[:, j] =
            Mt[1, j] * u[:, j] * g1[:, j] + Mt[2, j] * u[:, j] * h1[:, j]
        @. fh2[:, j] =
            Mt[1, j] * u[:, j] * g2[:, j] + Mt[2, j] * u[:, j] * h2[:, j]
        @. fh3[:, j] =
            Mt[1, j] * u[:, j] * g3[:, j] + Mt[2, j] * u[:, j] * h3[:, j]
    end

end


"""
Unified gas kinetic scheme (UGKS) flux

> @param[in] : particle distribution functions and their slopes at left/right sides of interface
> @param[in] : particle velocity quadrature points and weights
> @param[in] : time step

< @return : flux of particle distribution function and its velocity moments on conservative variables

"""

# ------------------------------------------------------------
# 2D2F flux
# ------------------------------------------------------------
function flux_ugks!(
    fw::AbstractArray{<:AbstractFloat,1},
    fh::AbstractArray{<:AbstractFloat,2},
    fb::AbstractArray{<:AbstractFloat,2},
    wL::AbstractArray{<:Real,1},
    hL::AbstractArray{<:AbstractFloat,2},
    bL::AbstractArray{<:AbstractFloat,2},
    wR::AbstractArray{<:Real,1},
    hR::AbstractArray{<:AbstractFloat,2},
    bR::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
    inK::Real,
    γ::Real,
    visRef::Real,
    visIdx::Real,
    pr::Real,
    dt::Real,
    dxL::Real,
    dxR::Real,
    len::Real,
    shL = zeros(eltype(hL), axes(hL))::AbstractArray{<:AbstractFloat,2},
    sbL = zeros(eltype(bL), axes(bL))::AbstractArray{<:AbstractFloat,2},
    shR = zeros(eltype(hR), axes(hR))::AbstractArray{<:AbstractFloat,2},
    sbR = zeros(eltype(bR), axes(bR))::AbstractArray{<:AbstractFloat,2},
)

    #--- reconstruct initial distribution ---#
    δ = heaviside.(u)

    h = @. hL * δ + hR * (1.0 - δ)
    b = @. bL * δ + bR * (1.0 - δ)
    sh = @. shL * δ + shR * (1.0 - δ)
    sb = @. sbL * δ + sbR * (1.0 - δ)

    #--- construct interface variables ---#
    #w = moments_conserve(h, b, u, v, ω)
    #prim = conserve_prim(w, γ)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)
    Mu1, Mv1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Muv1 = moments_conserve(MuL1, Mv1, Mxi1, 0, 0, 0)
    Mu2, Mv2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)
    Muv2 = moments_conserve(MuR2, Mv2, Mxi2, 0, 0, 0)
    w = @. primL[1] * Muv1 + primR[1] * Muv2
    prim = conserve_prim(w, γ)

    aL = pdf_slope(prim, (w .- wL) ./ dxL, inK)
    aR = pdf_slope(prim, (wR .- w) ./ dxR, inK)

    Mu, Mv, Mxi, MuL, MuR = gauss_moments(prim, inK)
    MauL = moments_conserve_slope(aL, MuL, Mv, Mxi, 1, 0)
    MauR = moments_conserve_slope(aR, MuR, Mv, Mxi, 1, 0)
    aT = pdf_slope(prim, -prim[1] .* (MauL .+ MauR), inK)

    #--- calculate integral time constants ---#
    tau = vhs_collision_time(prim, visRef, visIdx)

    Mt = zeros(5)
    Mt[4] = tau * (1.0 - exp(-dt / tau)) # f0
    Mt[5] = -tau * dt * exp(-dt / tau) + tau * Mt[4]
    Mt[1] = dt - Mt[4] # M0
    Mt[2] = -tau * Mt[1] + Mt[5]
    Mt[3] = dt^2 / 2.0 - tau * Mt[1]

    # --- calculate flux from M0 ---#
    Muv = moments_conserve(Mu, Mv, Mxi, 1, 0, 0)
    MauL = moments_conserve_slope(aL, MuL, Mv, Mxi, 2, 0)
    MauR = moments_conserve_slope(aR, MuR, Mv, Mxi, 2, 0)
    MauT = moments_conserve_slope(aT, Mu, Mv, Mxi, 1, 0)

    @. fw =
        Mt[1] * prim[1] * Muv +
        Mt[2] * prim[1] * (MauL + MauR) +
        Mt[3] * prim[1] * MauT

    # --- calculate flux from f0 ---#
    H = maxwellian(u, v, prim)
    B = H .* inK ./ (2.0 * prim[end])

    fw[1] += Mt[4] * sum(ω .* u .* h) - Mt[5] * sum(ω .* u .^ 2 .* sh)
    fw[2] += Mt[4] * sum(ω .* u .^ 2 .* h) - Mt[5] * sum(ω .* u .^ 2 .* sh)
    fw[3] += Mt[4] * sum(ω .* v .* u .* h) - Mt[5] * sum(ω .* v .* u .^ 2 .* sh)
    fw[4] +=
        Mt[4] *
        0.5 *
        (sum(ω .* u .* (u .^ 2 .+ v .^ 2) .* h) + sum(ω .* u .* b)) -
        Mt[5] *
        0.5 *
        (sum(ω .* u .^ 2 .* (u .^ 2 .+ v .^ 2) .* sh) + sum(ω .* u .^ 2 .* sb))

    @. fh =
        Mt[1] * u * H +
        Mt[2] *
        u^2 *
        (
            aL[1] * H +
            aL[2] * u * H +
            aL[3] * v * H +
            0.5 * aL[4] * ((u^2 + v^2) * H + B)
        ) *
        δ +
        Mt[2] *
        u^2 *
        (
            aR[1] * H +
            aR[2] * u * H +
            aR[3] * v * H +
            0.5 * aR[4] * ((u^2 + v^2) * H + B)
        ) *
        (1.0 - δ) +
        Mt[3] *
        u *
        (
            aT[1] * H +
            aT[2] * u * H +
            aT[3] * v * H +
            0.5 * aT[4] * ((u^2 + v^2) * H + B)
        ) +
        Mt[4] * u * h - Mt[5] * u^2 * sh
    @. fb =
        Mt[1] * u * B +
        Mt[2] *
        u^2 *
        (
            aL[1] * B +
            aL[2] * u * B +
            aL[3] * v * B +
            0.5 * aL[4] * ((u^2 + v^2) * B + Mxi[2] * H)
        ) *
        δ +
        Mt[2] *
        u^2 *
        (
            aR[1] * B +
            aR[2] * u * B +
            aR[3] * v * B +
            0.5 * aR[4] * ((u^2 + v^2) * B + Mxi[2] * H)
        ) *
        (1.0 - δ) +
        Mt[3] *
        u *
        (
            aT[1] * B +
            aT[2] * u * B +
            aT[3] * v * B +
            0.5 * aT[4] * ((u^2 + v^2) * B + Mxi[2] * H)
        ) +
        Mt[4] * u * b - Mt[5] * u^2 * sb

    # multiply interface length
    fw .*= len
    fh .*= len
    fb .*= len

end


"""
Maxwell's diffusive boundary flux

> @param[in] : particle distribution functions and their slopes at left/right sides of interface
> @param[in] : particle velocity quadrature points and weights
> @param[in] : time step

< @return : flux of particle distribution function and its velocity moments on conservative variables

"""

function flux_boundary_maxwell!(
    fw::AbstractArray{<:AbstractFloat,1},
    fh::AbstractArray{<:AbstractFloat,2},
    fb::AbstractArray{<:AbstractFloat,2},
    bc::Array{<:Real,1},
    h::AbstractArray{<:AbstractFloat,2},
    b::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
    inK::Real,
    dt::Real,
    len::Real,
    rot = 1::Real,
)

    @assert length(bc) == 4

    δ = heaviside.(u .* rot)
    SF = sum(ω .* u .* h .* (1.0 .- δ))
    SG =
        (bc[end] / π) * sum(
            ω .* u .*
            exp.(-bc[end] .* ((u .- bc[2]) .^ 2 .+ (v .- bc[3]) .^ 2)) .* δ,
        )
    prim = [-SF / SG; bc[2:end]]

    H = maxwellian(u, v, prim)
    B = H .* inK ./ (2.0 * prim[end])

    hWall = H .* δ .+ h .* (1.0 .- δ)
    bWall = B .* δ .+ b .* (1.0 .- δ)

    fw[1] = discrete_moments(hWall, u, ω, 1) * len * dt
    fw[2] = discrete_moments(hWall, u, ω, 2) * len * dt
    fw[3] = discrete_moments(hWall .* u, v, ω, 1) * len * dt
    fw[4] =
        (
            0.5 * discrete_moments(hWall .* (u .^ 2 .+ v .^ 2), u, ω, 1) +
            0.5 * discrete_moments(bWall, u, ω, 1)
        ) *
        len *
        dt

    @. fh = u * hWall * len * dt
    @. fb = u * bWall * len * dt

end


"""
Wave propagation method for the Maxwell's equations

> @param[in] : variables in left-left, left, right, and right-right cells
> @param[in] : eigenmatrix (A), eigenvalue (D)
> @param[in] : speed of light (sol)
> @param[in] : auxiliary parameters (χₑ, νᵦ)

< @return : flux of electromagnetic fields

"""

function flux_em!(
    femL::AbstractArray{<:AbstractFloat,1},
    femR::AbstractArray{<:AbstractFloat,1},
    ELL::AbstractArray{<:AbstractFloat,1},
    BLL::AbstractArray{<:AbstractFloat,1},
    EL::AbstractArray{<:AbstractFloat,1},
    BL::AbstractArray{<:AbstractFloat,1},
    ER::AbstractArray{<:AbstractFloat,1},
    BR::AbstractArray{<:AbstractFloat,1},
    ERR::AbstractArray{<:AbstractFloat,1},
    BRR::AbstractArray{<:AbstractFloat,1},
    ϕL::AbstractFloat,
    ϕR::AbstractFloat,
    ψL::AbstractFloat,
    ψR::AbstractFloat,
    dxL::AbstractFloat,
    dxR::AbstractFloat,
    A1p::AbstractArray{<:AbstractFloat,2},
    A1n::AbstractArray{<:AbstractFloat,2},
    D1::AbstractArray{<:AbstractFloat,1},
    sol::AbstractFloat,
    χ::AbstractFloat,
    ν::AbstractFloat,
    dt::AbstractFloat,
)

    @assert length(femL) == length(femR) == 8

    slop = zeros(8, 8)
    slop[3, 1] = -0.5 * sol^2 * (BR[2] - BL[2]) + 0.5 * sol * (ER[3] - EL[3])
    slop[5, 1] = 0.5 * sol * (BR[2] - BL[2]) - 0.5 * (ER[3] - EL[3])
    slop[2, 2] = 0.5 * sol^2 * (BR[3] - BL[3]) + 0.5 * sol * (ER[2] - EL[2])
    slop[6, 2] = 0.5 * sol * (BR[3] - BL[3]) + 0.5 * (ER[2] - EL[2])
    slop[1, 3] = 0.5 * sol * χ * (ER[1] - EL[1])
    slop[7, 3] = 0.5 * χ * (ER[1] - EL[1])
    slop[4, 4] = 0.5 * sol * ν * (BR[1] - BR[1])
    slop[8, 4] = 0.5 * sol^2 * χ * (BR[1] - BR[1])
    slop[3, 5] = -0.5 * sol^2 * (BR[2] - BL[2]) - 0.5 * sol * (ER[3] - EL[3])
    slop[5, 5] = -0.5 * sol * (BR[2] - BL[2]) - 0.5 * (ER[3] - EL[3])
    slop[2, 6] = 0.5 * sol^2 * (BR[3] - BL[3]) - 0.5 * sol * (ER[2] - EL[2])
    slop[6, 6] = -0.5 * sol * (BR[3] - BL[3]) + 0.5 * (ER[2] - EL[2])
    slop[1, 7] = -0.5 * sol * χ * (ER[1] - EL[1])
    slop[7, 7] = 0.5 * χ * (ER[1] - EL[1])
    slop[4, 8] = -0.5 * sol * ν * (BR[1] - BR[1])
    slop[8, 8] = 0.5 * sol^2 * χ * (BR[1] - BR[1])

    limiter = zeros(8, 8)
    limiter[3, 1] =
        -0.5 * sol^2 * (BL[2] - BLL[2]) + 0.5 * sol * (EL[3] - ELL[3])
    limiter[5, 1] = 0.5 * sol * (BL[2] - BLL[2]) - 0.5 * (EL[3] - ELL[3])
    limiter[2, 2] =
        0.5 * sol^2 * (BL[3] - BLL[3]) + 0.5 * sol * (EL[2] - ELL[2])
    limiter[6, 2] = 0.5 * sol * (BL[3] - BLL[3]) + 0.5 * (EL[2] - ELL[2])
    limiter[1, 3] = 0.5 * sol * χ * (EL[1] - ELL[1])
    limiter[7, 3] = 0.5 * χ * (EL[1] - ELL[1])
    limiter[4, 4] = 0.5 * sol * ν * (BL[1] - BL[1])
    limiter[8, 4] = 0.5 * sol^2 * χ * (BL[1] - BL[1])
    limiter[3, 5] =
        -0.5 * sol^2 * (BRR[2] - BR[2]) - 0.5 * sol * (ERR[3] - ER[3])
    limiter[5, 5] = -0.5 * sol * (BRR[2] - BR[2]) - 0.5 * (ERR[3] - ER[3])
    limiter[2, 6] =
        0.5 * sol^2 * (BRR[3] - BR[3]) - 0.5 * sol * (ERR[2] - ER[2])
    limiter[6, 6] = -0.5 * sol * (BRR[3] - BR[3]) + 0.5 * (ERR[2] - ER[2])
    limiter[1, 7] = -0.5 * sol * χ * (ERR[1] - ER[1])
    limiter[7, 7] = 0.5 * χ * (ERR[1] - ER[1])
    limiter[4, 8] = -0.5 * sol * ν * (BRR[1] - BRR[1])
    limiter[8, 8] = 0.5 * sol^2 * χ * (BRR[1] - BRR[1])

    for i = 1:8
        limiter_theta =
            sum(slop[:, i] .* limiter[:, i]) / (sum(slop[:, i] .^ 2) + 1.e-7)
        slop[:, i] .*= max(
            0.0,
            min(min((1.0 + limiter_theta) / 2.0, 2.0), 2.0 * limiter_theta),
        )
    end

    for i = 1:8
        femL[i] =
            sum(A1n[i, 1:3] .* (ER .- EL)) +
            sum(A1n[i, 4:6] .* (BR .- BL)) +
            A1n[i, 7] * (ϕR - ϕL) +
            A1n[i, 8] * (ψR - ψL) +
            0.5 * sum(
                fortsign.(1.0, D1) .*
                (1.0 .- dt ./ (0.5 * (dxL + dxR)) .* abs.(D1)) .* slop[i, :],
            )
        femR[i] =
            sum(A1p[i, 1:3] .* (ER .- EL)) +
            sum(A1p[i, 4:6] .* (BR .- BL)) +
            A1p[i, 7] * (ϕR - ϕL) +
            A1p[i, 8] * (ψR - ψL) -
            0.5 * sum(
                fortsign.(1.0, D1) .*
                (1.0 .- dt ./ (0.5 * (dxL + dxR)) .* abs.(D1)) .* slop[i, :],
            )
    end

end
