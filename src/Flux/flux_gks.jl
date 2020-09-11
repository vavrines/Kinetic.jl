"""
Gas kinetic Navier-Stokes flux

`flux_gks(uL, uR, μ, dt, dxL, dxR, suL, suR, a)`

* @arg: conservative scalars and their left/right slopes
* @arg: molecular and thermodynamic parameters
* @arg: time step and cell size
* @return: scalar flux

"""
function flux_gks(u::Real, μ::Real, dt::Real, su = 0.0::Real, a = 0::Real)

    prim = ifelse(a == 0, conserve_prim(u), conserve_prim(u, a))

    Mu, MuL, MuR = gauss_moments(prim)

    tau = 2.0 * μ

    fa = pdf_slope(u, su)
    Δ = -prim[1] * moments_conserve_slope(fa, Mu, 1)
    faT = pdf_slope(u, Δ)

    # time-integration constants
    Mt = zeros(5)
    Mt[4] = dt
    Mt[5] = -tau * dt * exp(-dt / tau) + tau * Mt[4]


    # flux related to upwind distribution
    Muv = moments_conserve(Mu, 1)
    Mau = moments_conserve_slope(fa, Mu, 2)
    MauT = moments_conserve_slope(faT, Mu, 1)

    fw =
        Mt[4] * prim[1] * Muv - (Mt[5] + tau * Mt[4]) * prim[1] * Mau -
        tau * Mt[4] * prim[1] * MauT

    return fw / dt

end


function flux_gks(
    uL::Real,
    uR::Real,
    μ::Real,
    dt::Real,
    dxL::Real,
    dxR::Real,
    suL = 0.0::Real,
    suR = 0.0::Real,
    a = 0::Real,
)

    primL = ifelse(a == 0, conserve_prim(uL), conserve_prim(uL, a))
    primR = ifelse(a == 0, conserve_prim(uR), conserve_prim(uR, a))

    Mu1, MuL1, MuR1 = gauss_moments(primL)
    Mu2, MuL2, MuR2 = gauss_moments(primR)

    u = primL[1] * moments_conserve(MuL1, 0) + primR[1] * moments_conserve(MuR2, 0)
    prim = ifelse(a == 0, conserve_prim(u), conserve_prim(u, a))
    tau = 2.0 * abs(uL - uR) / (abs(uL) + abs(uR)) * dt + 2.0 * μ

    faL = pdf_slope(uL, suL)
    Δ = -primL[1] * moments_conserve_slope(faL, Mu1, 1)
    faTL = pdf_slope(uL, Δ)

    faR = pdf_slope(uR, suR)
    Δ = -primR[1] * moments_conserve_slope(faR, Mu2, 1)
    faTR = pdf_slope(uR, Δ)

    Mu, MuL, MuR = gauss_moments(prim)
    sw0L = (u - uL) / dxL
    sw0R = (uR - u) / dxR
    gaL = pdf_slope(u, sw0L)
    gaR = pdf_slope(u, sw0R)
    Δ =
        -prim[1] *
        (moments_conserve_slope(gaL, MuL, 1) + moments_conserve_slope(gaR, MuR, 1))
    # sw = (uR - uL) / (dxL + dxR)
    # ga = pdf_slope(u, sw)
    # Δ = -prim[1] .* moments_conserve_slope(ga, Mu, 1)
    gaT = pdf_slope(u, Δ)

    # time-integration constants
    Mt = zeros(5)
    Mt[4] = tau * (1.0 - exp(-dt / tau))
    Mt[5] = -tau * dt * exp(-dt / tau) + tau * Mt[4]
    Mt[1] = dt - Mt[4]
    Mt[2] = -tau * Mt[1] + Mt[5]
    Mt[3] = 0.5 * dt^2 - tau * Mt[1]

    # flux related to central distribution
    Muv = moments_conserve(Mu, 1)
    MauL = moments_conserve_slope(gaL, MuL, 2)
    MauR = moments_conserve_slope(gaR, MuR, 2)
    # Mau = moments_conserve_slope(ga, MuR, 2)
    MauT = moments_conserve_slope(gaT, Mu, 1)

    fw = Mt[1] * prim[1] * Muv + Mt[2] * prim[1] * (MauL + MauR) + Mt[3] * prim[1] * MauT
    # fw = Mt[1] * prim[1] * Muv + Mt[2] * prim[1] * Mau + Mt[3] * prim[1] * MauT

    # flux related to upwind distribution
    MuvL = moments_conserve(MuL1, 1)
    MauL = moments_conserve_slope(faL, MuL1, 2)
    MauLT = moments_conserve_slope(faTL, MuL1, 1)

    MuvR = moments_conserve(MuR2, 1)
    MauR = moments_conserve_slope(faR, MuR2, 2)
    MauRT = moments_conserve_slope(faTR, MuR2, 1)

    fw +=
        Mt[4] * primL[1] * MuvL - (Mt[5] + tau * Mt[4]) * primL[1] * MauL -
        tau * Mt[4] * primL[1] * MauLT + Mt[4] * primR[1] * MuvR -
        (Mt[5] + tau * Mt[4]) * primR[1] * MauR - tau * Mt[4] * primR[1] * MauRT
    # fw += Mt[4] * primL[1] * MuvL + Mt[4] * primR[1] * MuvR

    return fw

end


"""
Gas kinetic Navier-Stokes flux

* 1D: flux_gks!(fw, wL, wR, γ, K, μᵣ, ω, dt, dx, swL, swR)
* 2D: flux_gks!(fw, wL, wR, γ, K, μᵣ, ω, dt, dx, dy, swL, swR)

* @param[in]: conservative variables and their left/right slopes
* @param[in]: molecular and thermodynamic parameters
* @param[in]: time step and cell size

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
    sw0L = (w .- wL) ./ (0.5 * dx)
    sw0R = (wR .- w) ./ (0.5 * dx)
    gaL = pdf_slope(prim, sw0L, inK)
    gaR = pdf_slope(prim, sw0R, inK)
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
        (Mt[5] + tau * Mt[4]) * primR[1] * MauR - tau * Mt[4] * primR[1] * MauRT
    # @. fw += Mt[4] * primL[1] * MuvL + Mt[4] * primR[1] * MuvR

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
    sw0L = (w .- wL) ./ (0.5 * dx)
    sw0R = (wR .- w) ./ (0.5 * dx)
    gaL = pdf_slope(prim, sw0L, inK)
    gaR = pdf_slope(prim, sw0R, inK)
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
        (Mt[5] + tau * Mt[4]) * primR[1] * MauR - tau * Mt[4] * primR[1] * MauRT
    # @. fw += Mt[4] * primL[1] * MuvL + Mt[4] * primR[1] * MuvR

    fw .*= dy

end


function flux_gks!(
    fw::AbstractArray{<:AbstractFloat,1},
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
    dt::Real,
    dxL::Real,
    dxR::Real,
    swL = zeros(eltype(fw), axes(wL))::AbstractArray{<:AbstractFloat,1},
    swR = zeros(eltype(fw), axes(wR))::AbstractArray{<:AbstractFloat,1},
    shL = zeros(eltype(hL), axes(hL))::AbstractArray{<:AbstractFloat,1},
    sbL = zeros(eltype(bL), axes(bL))::AbstractArray{<:AbstractFloat,1},
    shR = zeros(eltype(hR), axes(hR))::AbstractArray{<:AbstractFloat,1},
    sbR = zeros(eltype(bR), axes(bR))::AbstractArray{<:AbstractFloat,1},
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
        vhs_collision_time(prim, visRef, visIdx) +
        2.0 * dt * abs(primL[1] / primL[end] - primR[1] / primR[end]) /
        (primL[1] / primL[end] + primR[1] / primR[end])

    faL = pdf_slope(primL, swL, inK)
    sw = -primL[1] .* moments_conserve_slope(faL, Mu1, Mxi1, 1)
    faTL = pdf_slope(primL, sw, inK)

    faR = pdf_slope(primR, swR, inK)
    sw = -primR[1] .* moments_conserve_slope(faR, Mu2, Mxi2, 1)
    faTR = pdf_slope(primR, sw, inK)

    Mu, Mxi, MuL, MuR = gauss_moments(prim, inK)
    sw0L = (w .- wL) ./ dxL
    sw0R = (wR .- w) ./ dxR
    gaL = pdf_slope(prim, sw0L, inK)
    gaR = pdf_slope(prim, sw0R, inK)
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
        (Mt[5] + tau * Mt[4]) * primR[1] * MauR - tau * Mt[4] * primR[1] * MauRT
    # @. fw += Mt[4] * primL[1] * MuvL + Mt[4] * primR[1] * MuvR

    #--- fluxes of distribution functions ---#
    δ = heaviside.(u)

    HL = maxwellian(u, primL)
    BL = HL .* inK ./ (2.0 * primL[end])
    HR = maxwellian(u, primR)
    BR = HR .* inK ./ (2.0 * primR[end])

    H = maxwellian(u, prim)
    B = H .* inK ./ (2.0 * prim[end])

    @. fh =
        Mt[1] * u * H +
        Mt[2] * u^2 * (gaL[1] * H + gaL[2] * u * H + 0.5 * gaL[3] * (u^2 * H + B)) * δ +
        Mt[2] *
        u^2 *
        (gaR[1] * H + gaR[2] * u * H + 0.5 * gaR[3] * (u^2 * H + B)) *
        (1.0 - δ) +
        Mt[3] * u * (gaT[1] * H + gaT[2] * u * H + 0.5 * gaT[3] * (u^2 * H + B)) +
        Mt[4] * u * HL * δ -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (faL[1] * HL + faL[2] * u * HL + 0.5 * faL[3] * (u^2 * HL + BL)) *
        δ + Mt[4] * u * HR * (1.0 - δ) -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (faR[1] * HR + faR[2] * u * HR + 0.5 * faR[3] * (u^2 * HR + BR)) *
        (1.0 - δ) -
        tau *
        Mt[4] *
        u *
        (faTL[1] * HL + faTL[2] * u * HL + 0.5 * faTL[3] * (u^2 * HL + BL)) *
        δ -
        tau *
        Mt[4] *
        u *
        (faTR[1] * HR + faTR[2] * u * HR + 0.5 * faTR[3] * (u^2 * HR + BR)) *
        (1.0 - δ)

    @. fb =
        Mt[1] * u * B +
        Mt[2] *
        u^2 *
        (gaL[1] * B + gaL[2] * u * B + 0.5 * gaL[3] * (u^2 * B + Mxi[2] * H)) *
        δ +
        Mt[2] *
        u^2 *
        (gaR[1] * B + gaR[2] * u * B + 0.5 * gaR[3] * (u^2 * B + Mxi[2] * H)) *
        (1.0 - δ) +
        Mt[3] * u * (gaT[1] * H + gaT[2] * u * H + 0.5 * gaT[3] * (u^2 * B + Mxi[2] * H)) +
        Mt[4] * u * BL * δ -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (faL[1] * BL + faL[2] * u * BL + 0.5 * faL[3] * (u^2 * BL + Mxi[2] * HL)) *
        δ + Mt[4] * u * BR * (1.0 - δ) -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (faR[1] * BR + faR[2] * u * BR + 0.5 * faR[3] * (u^2 * BR + Mxi[2] * HR)) *
        (1.0 - δ) -
        tau *
        Mt[4] *
        u *
        (faTL[1] * BL + faTL[2] * u * BL + 0.5 * faTL[3] * (u^2 * BL + Mxi[2] * HL)) *
        δ -
        tau *
        Mt[4] *
        u *
        (faTR[1] * BR + faTR[2] * u * BR + 0.5 * faTR[3] * (u^2 * BR + Mxi[2] * HR)) *
        (1.0 - δ)

end


"""
Unified gas kinetic scheme (UGKS) flux

* @arg: particle distribution functions and their slopes at left/right sides of interface
* @arg: particle velocity quadrature points and weights
* @arg: time step

"""
function flux_ugks!(
    fw::AbstractArray{<:AbstractFloat,1},
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
    dxL::Real,
    dxR::Real,
    shL = zeros(eltype(hL), axes(hL))::AbstractArray{<:AbstractFloat,1},
    sbL = zeros(eltype(bL), axes(bL))::AbstractArray{<:AbstractFloat,1},
    shR = zeros(eltype(hR), axes(hR))::AbstractArray{<:AbstractFloat,1},
    sbR = zeros(eltype(bR), axes(bR))::AbstractArray{<:AbstractFloat,1},
) # 1D2F flux

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
    Mu1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Muv1 = moments_conserve(MuL1, Mxi1, 0, 0)
    Mu2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)
    Muv2 = moments_conserve(MuR2, Mxi2, 0, 0)

    w = @. primL[1] * Muv1 + primR[1] * Muv2
    prim = conserve_prim(w, γ)
    aL = pdf_slope(prim, (w .- wL) ./ dxL, inK)
    aR = pdf_slope(prim, (wR .- w) ./ dxR, inK)

    Mu, Mxi, MuL, MuR = gauss_moments(prim, inK)
    MauL = moments_conserve_slope(aL, MuL, Mxi, 1)
    MauR = moments_conserve_slope(aR, MuR, Mxi, 1)
    aT = pdf_slope(prim, -prim[1] .* (MauL .+ MauR), inK)

    #--- calculate integral time constants ---#
    tau = vhs_collision_time(prim, visRef, visIdx)

    Mt = zeros(5)
    Mt[4] = tau * (1.0 - exp(-dt / tau)) # f0
    Mt[5] = -tau * dt * exp(-dt / tau) + tau * Mt[4]
    Mt[1] = dt - Mt[4] # M0
    Mt[2] = -tau * Mt[1] + Mt[5]
    Mt[3] = dt^2 / 2.0 - tau * Mt[1]

    #--- calculate flux from M0 ---#
    Muv = moments_conserve(Mu, Mv, Mxi, 1, 0, 0)
    MauL = moments_conserve_slope(aL, MuL, Mv, Mxi, 2, 0)
    MauR = moments_conserve_slope(aR, MuR, Mv, Mxi, 2, 0)
    MauT = moments_conserve_slope(aT, Mu, Mv, Mxi, 1, 0)

    @. fw = Mt[1] * prim[1] * Muv + Mt[2] * prim[1] * (MauL + MauR) + Mt[3] * prim[1] * MauT

    #--- calculate flux from f0 ---#
    H = maxwellian(u, v, prim)
    B = H .* inK ./ (2.0 * prim[end])

    fw[1] += Mt[4] * sum(ω .* u .* h) - Mt[5] * sum(ω .* u .^ 2 .* sh)
    fw[2] += Mt[4] * sum(ω .* u .^ 2 .* h) - Mt[5] * sum(ω .* u .^ 3 .* sh)
    fw[3] +=
        Mt[4] * 0.5 * (sum(ω .* u .^ 3 .* h) + sum(ω .* u .* b)) -
        Mt[5] * 0.5 * (sum(ω .* u .^ 4 .* sh) + sum(ω .* u .^ 2 .* sb))

    @. fh =
        Mt[1] * u * H +
        Mt[2] * u^2 * (aL[1] * H + aL[2] * u * H + 0.5 * aL[3] * (u^2 * H + B)) * δ +
        Mt[2] *
        u^2 *
        (aR[1] * H + aR[2] * u * H + 0.5 * aR[3] * (u^2 * H + B)) *
        (1.0 - δ) +
        Mt[3] * u * (aT[1] * H + aT[2] * u * H + 0.5 * aT[3] * (u^2 * H + B)) +
        Mt[4] * u * h - Mt[5] * u^2 * sh
    @. fb =
        Mt[1] * u * B +
        Mt[2] *
        u^2 *
        (aL[1] * B + aL[2] * u * B + 0.5 * aL[3] * (u^2 * B + Mxi[2] * H)) *
        δ +
        Mt[2] *
        u^2 *
        (aR[1] * B + aR[2] * u * B + 0.5 * aR[3] * (u^2 * B + Mxi[2] * H)) *
        (1.0 - δ) +
        Mt[3] * u * (aT[1] * B + aT[2] * u * B + 0.5 * aT[3] * (u^2 * B + Mxi[2] * H)) +
        Mt[4] * u * b - Mt[5] * u^2 * sb

end


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
) # 2D2F flux

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

    @. fw = Mt[1] * prim[1] * Muv + Mt[2] * prim[1] * (MauL + MauR) + Mt[3] * prim[1] * MauT

    # --- calculate flux from f0 ---#
    H = maxwellian(u, v, prim)
    B = H .* inK ./ (2.0 * prim[end])

    fw[1] += Mt[4] * sum(ω .* u .* h) - Mt[5] * sum(ω .* u .^ 2 .* sh)
    fw[2] += Mt[4] * sum(ω .* u .^ 2 .* h) - Mt[5] * sum(ω .* u .^ 3 .* sh)
    fw[3] += Mt[4] * sum(ω .* v .* u .* h) - Mt[5] * sum(ω .* v .* u .^ 2 .* sh)
    fw[4] +=
        Mt[4] * 0.5 * (sum(ω .* u .* (u .^ 2 .+ v .^ 2) .* h) + sum(ω .* u .* b)) -
        Mt[5] *
        0.5 *
        (sum(ω .* u .^ 2 .* (u .^ 2 .+ v .^ 2) .* sh) + sum(ω .* u .^ 2 .* sb))

    @. fh =
        Mt[1] * u * H +
        Mt[2] *
        u^2 *
        (aL[1] * H + aL[2] * u * H + aL[3] * v * H + 0.5 * aL[4] * ((u^2 + v^2) * H + B)) *
        δ +
        Mt[2] *
        u^2 *
        (aR[1] * H + aR[2] * u * H + aR[3] * v * H + 0.5 * aR[4] * ((u^2 + v^2) * H + B)) *
        (1.0 - δ) +
        Mt[3] *
        u *
        (aT[1] * H + aT[2] * u * H + aT[3] * v * H + 0.5 * aT[4] * ((u^2 + v^2) * H + B)) +
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

* @param[in]: particle distribution functions and their slopes at left/right sides of interface
* @param[in]: particle velocity quadrature points and weights
* @param[in]: time step

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
        (bc[end] / π) *
        sum(ω .* u .* exp.(-bc[end] .* ((u .- bc[2]) .^ 2 .+ (v .- bc[3]) .^ 2)) .* δ)
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
