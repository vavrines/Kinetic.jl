"""
Gas kinetic Navier-Stokes flux

    flux_gks(u::Real, μ::Real, dt::Real, su = 0.0::Real, a = 0::Real)
    flux_gks(uL::Real, uR::Real, μ::Real, dt::Real, dxL::Real,
        dxR::Real, suL = 0.0::Real, suR = 0.0::Real, a = 0::Real)

* @args: conservative scalars and their slopes
* @args: viscosity
* @args: time step and cell size
* @return: scalar flux

"""
function flux_gks(u::T, μ, dt, su = 0.0, a = 0) where {T<:Real}

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
    uL::T,
    uR::T,
    μ,
    dt,
    dxL,
    dxR,
    suL = 0.0,
    suR = 0.0,
    a = 0,
) where {T<:Real}

    primL = ifelse(a == 0, conserve_prim(uL), conserve_prim(uL, a))
    primR = ifelse(a == 0, conserve_prim(uR), conserve_prim(uR, a))

    Mu1, MuL1, MuR1 = gauss_moments(primL)
    Mu2, MuL2, MuR2 = gauss_moments(primR)

    u =
        primL[1] * moments_conserve(MuL1, 0) +
        primR[1] * moments_conserve(MuR2, 0)
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
        -prim[1] * (
            moments_conserve_slope(gaL, MuL, 1) +
            moments_conserve_slope(gaR, MuR, 1)
        )
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

    fw =
        Mt[1] * prim[1] * Muv +
        Mt[2] * prim[1] * (MauL + MauR) +
        Mt[3] * prim[1] * MauT
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

    1D: flux_gks!(fw, wL, wR, γ, K, μᵣ, ω, dt, dx, swL, swR)
    2D: flux_gks!(fw, wL, wR, γ, K, μᵣ, ω, dt, dx, dy, swL, swR)

* @args: conservative variables and their left/right slopes
* @args: molecular and thermodynamic parameters
* @args: time step and cell size

"""
function flux_gks!(
    fw::X,
    wL::Y,
    wR::Y,
    γ::Real,
    inK::Real,
    μᵣ::Real,
    ω::Real,
    dt::Real,
    dx::Real,
    swL = zeros(eltype(fw), axes(wL))::X,
    swR = zeros(eltype(fw), axes(wR))::X,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,1},
}

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

    return nothing

end


function flux_gks!(
    fw::X,
    wL::Y,
    wR::Y,
    γ::Real,
    inK::Real,
    μᵣ::Real,
    ω::Real,
    dt::Real,
    dx::Real,
    dy::Real,
    swL = zeros(eltype(fw), axes(wL))::X,
    swR = zeros(eltype(fw), axes(wR))::X,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,1},
}

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

    return nothing

end

#--- 1D2V ---#
function flux_gks!(
    fw::T1,
    fh::T2,
    fb::T2,
    wL::T3,
    hL::T4,
    bL::T4,
    wR::T3,
    hR::T4,
    bR::T4,
    u::T5,
    ω::T5,
    inK,
    γ,
    visRef,
    visIdx,
    dt,
    dxL,
    dxR,
    swL = zeros(eltype(fw), axes(wL))::T1,
    swR = zeros(eltype(fw), axes(wR))::T1,
    shL = zeros(eltype(hL), axes(hL))::T4,
    sbL = zeros(eltype(bL), axes(bL))::T4,
    shR = zeros(eltype(hR), axes(hR))::T4,
    sbR = zeros(eltype(bR), axes(bR))::T4,
) where {
    T1<:AbstractArray{<:AbstractFloat,1},
    T2<:AbstractArray{<:AbstractFloat,1},
    T3<:AbstractArray{<:Real,1},
    T4<:AbstractArray{<:AbstractFloat,1},
    T5<:AbstractArray{<:AbstractFloat,1},
}

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
        Mt[2] *
        u^2 *
        (gaL[1] * H + gaL[2] * u * H + 0.5 * gaL[3] * (u^2 * H + B)) *
        δ +
        Mt[2] *
        u^2 *
        (gaR[1] * H + gaR[2] * u * H + 0.5 * gaR[3] * (u^2 * H + B)) *
        (1.0 - δ) +
        Mt[3] *
        u *
        (gaT[1] * H + gaT[2] * u * H + 0.5 * gaT[3] * (u^2 * H + B)) +
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
        Mt[3] *
        u *
        (gaT[1] * H + gaT[2] * u * H + 0.5 * gaT[3] * (u^2 * B + Mxi[2] * H)) +
        Mt[4] * u * BL * δ -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (
            faL[1] * BL +
            faL[2] * u * BL +
            0.5 * faL[3] * (u^2 * BL + Mxi[2] * HL)
        ) *
        δ + Mt[4] * u * BR * (1.0 - δ) -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (
            faR[1] * BR +
            faR[2] * u * BR +
            0.5 * faR[3] * (u^2 * BR + Mxi[2] * HR)
        ) *
        (1.0 - δ) -
        tau *
        Mt[4] *
        u *
        (
            faTL[1] * BL +
            faTL[2] * u * BL +
            0.5 * faTL[3] * (u^2 * BL + Mxi[2] * HL)
        ) *
        δ -
        tau *
        Mt[4] *
        u *
        (
            faTR[1] * BR +
            faTR[2] * u * BR +
            0.5 * faTR[3] * (u^2 * BR + Mxi[2] * HR)
        ) *
        (1.0 - δ)

    return nothing

end

#--- mixture ---#
function flux_gks!(
    fw::X,
    wL::Y,
    wR::Y,
    inK,
    γ,
    mi,
    ni,
    me,
    ne,
    Kn,
    dt,
    dxL,
    dxR,
    len,
    swL = zeros(eltype(fw), axes(wL))::X,
    swR = zeros(eltype(fw), axes(wR))::X,
) where {X<:AbstractArray{<:AbstractFloat,2},Y<:AbstractArray{<:Real,2}}

    primL = mixture_conserve_prim(wL, γ)
    primR = mixture_conserve_prim(wR, γ)

    Mu1, Mv1, Mw1, MuL1, MuR1 = mixture_gauss_moments(primL, inK)
    Mu2, Mv2, Mw2, MuL2, MuR2 = mixture_gauss_moments(primL, inK)

    w =
        primL[1] .* mixture_moments_conserve(MuL1, Mv1, Mw1, 0, 0, 0) .+
        primR[1] .* mixture_moments_conserve(MuR2, Mv2, Mw2, 0, 0, 0)
    prim = mixture_conserve_prim(w, γ)
    tau = aap_hs_collision_time(prim, mi, ni, me, ne, Kn)
    for i in eachindex(tau)
        tau[i] +=
            2.0 *
            dt *
            abs(primL[1, i] / primL[end, i] - primR[1, i] / primR[end, i]) /
            (primL[1, i] / primL[end, i] + primR[1, i] / primR[end, i])
    end
    prim = aap_hs_prim(prim, tau, mi, ni, me, ne, Kn) # pseudo primitive variables

    faL = mixture_pdf_slope(primL, swL, inK)
    sw = mixture_moments_conserve_slope(faL, Mu1, Mv1, Mw1, 1, 0, 0)
    for j in axes(sw, 2)
        sw[:, j] .*= -primL[1, j]
    end
    faTL = mixture_pdf_slope(primL, sw, inK)
    faR = mixture_pdf_slope(primR, swR, inK)
    sw = mixture_moments_conserve_slope(faR, Mu2, Mv2, Mw2, 1, 0, 0)
    for j in axes(sw, 2)
        sw[:, j] .*= -primR[1, j]
    end
    faTR = mixture_pdf_slope(primR, sw, inK)

    Mu, Mv, Mw, MuL, MuR = mixture_gauss_moments(prim, inK)
    sw0L = (w .- wL) ./ dxL
    sw0R = (wR .- w) ./ dxR
    gaL = mixture_pdf_slope(prim, sw0L, inK)
    gaR = mixture_pdf_slope(prim, sw0R, inK)
    sw =
        mixture_moments_conserve_slope(gaL, MuL, Mv, Mw, 1, 0, 0) .+
        mixture_moments_conserve_slope(gaR, MuR, Mv, Mw, 1, 0, 0)
    for j in axes(sw, 2)
        sw[:, j] .*= -prim[1, j]
    end
    gaT = mixture_pdf_slope(prim, sw, inK)

    # time-integration constants
    Mt = zeros(5, axes(fw, 2))
    for j in axes(Mt, 2)
        Mt[4, j] = tau[j] * (1.0 - exp(-dt / tau[j]))
        Mt[5, j] = -tau[j] * dt * exp(-dt / tau[j]) + tau[j] * Mt[4]
        Mt[1, j] = dt - Mt[4, j]
        Mt[2, j] = -tau[j] * Mt[1, j] + Mt[5, j]
        Mt[3, j] = 0.5 * dt^2 - tau[j] * Mt[1, j]
    end

    # flux related to central distribution
    Muv = mixture_moments_conserve(Mu, Mv, Mw, 1, 0, 0)
    MauL = mixture_moments_conserve_slope(gaL, MuL, Mv, Mw, 2, 0, 0)
    MauR = mixture_moments_conserve_slope(gaR, MuR, Mv, Mw, 2, 0, 0)
    MauT = mixture_moments_conserve_slope(gaT, Mu, Mv, Mw, 1, 0, 0)

    for j in axes(fw, 2)
        @. fw[:, j] =
            Mt[1, j] * prim[1, j] * Muv[:, j] +
            Mt[2, j] * prim[1, j] * (MauL[:, j] + MauR[:, j]) +
            Mt[3, j] * prim[1, j] * MauT[:, j]
    end

    # flux related to upwind distribution
    MuvL = mixture_moments_conserve(MuL1, Mv1, Mw1, 1, 0, 0)
    MauL = mixture_moments_conserve_slope(faL, MuL1, Mv1, Mw1, 2, 0, 0)
    MauLT = mixture_moments_conserve_slope(faTL, MuL1, Mv1, Mw1, 1, 0, 0)

    MuvR = mixture_moments_conserve(MuR2, Mv2, Mw2, 1, 0, 0)
    MauR = mixture_moments_conserve_slope(faR, MuR2, Mv2, Mw2, 2, 0, 0)
    MauRT = mixture_moments_conserve_slope(faTR, MuR2, Mv2, Mw2, 1, 0, 0)

    for j in axes(fw, 2)
        @. fw[:, j] +=
            Mt[4, j] * primL[1, j] * MuvL[:, j] -
            (Mt[5, j] + tau[j] * Mt[4, j]) * primL[1, j] * MauL[:, j] -
            tau[j] * Mt[4, j] * primL[1, j] * MauLT[:, j] +
            Mt[4, j] * primR[1, j] * MuvR[:, j] -
            (Mt[5, j] + tau[j] * Mt[4, j]) * primR[1, j] * MauR[:, j] -
            tau[j] * Mt[4, j] * primR[1, j] * MauRT[:, j]
    end

    @. fw .* len

    return nothing

end


"""
Unified gas kinetic scheme (UGKS)

* @args: particle distribution functions and their slopes at left/right sides of interface
* @args: particle velocity quadrature points and weights
* @args: time step

"""
function flux_ugks!(
    fw::T1,
    fh::T2,
    fb::T2,
    wL::T3,
    hL::T4,
    bL::T4,
    wR::T3,
    hR::T4,
    bR::T4,
    u::T5,
    ω::T5,
    inK,
    γ,
    visRef,
    visIdx,
    pr,
    dt,
    dxL,
    dxR,
    shL = zeros(eltype(hL), axes(hL))::T4,
    sbL = zeros(eltype(bL), axes(bL))::T4,
    shR = zeros(eltype(hR), axes(hR))::T4,
    sbR = zeros(eltype(bR), axes(bR))::T4,
) where {
    T1<:AbstractArray{<:AbstractFloat,1},
    T2<:AbstractArray{<:AbstractFloat,1},
    T3<:AbstractArray{<:Real,1},
    T4<:AbstractArray{<:AbstractFloat,1},
    T5<:AbstractArray{<:AbstractFloat,1},
} # 1D2F flux

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

    @. fw =
        Mt[1] * prim[1] * Muv +
        Mt[2] * prim[1] * (MauL + MauR) +
        Mt[3] * prim[1] * MauT

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
        Mt[2] *
        u^2 *
        (aL[1] * H + aL[2] * u * H + 0.5 * aL[3] * (u^2 * H + B)) *
        δ +
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
        Mt[3] *
        u *
        (aT[1] * B + aT[2] * u * B + 0.5 * aT[3] * (u^2 * B + Mxi[2] * H)) +
        Mt[4] * u * b - Mt[5] * u^2 * sb

    return nothing

end

function flux_ugks!(
    fw::T1,
    fh::T2,
    fb::T2,
    wL::T3,
    hL::T4,
    bL::T4,
    wR::T3,
    hR::T4,
    bR::T4,
    u::T5,
    v::T5,
    ω::T5,
    inK,
    γ,
    visRef,
    visIdx,
    pr,
    dt,
    dxL,
    dxR,
    len,
    shL = zeros(eltype(hL), axes(hL))::T4,
    sbL = zeros(eltype(bL), axes(bL))::T4,
    shR = zeros(eltype(hR), axes(hR))::T4,
    sbR = zeros(eltype(bR), axes(bR))::T4,
) where {
    T1<:AbstractArray{<:AbstractFloat,1},
    T2<:AbstractArray{<:AbstractFloat,2},
    T3<:AbstractArray{<:Real,1},
    T4<:AbstractArray{<:AbstractFloat,2},
    T5<:AbstractArray{<:AbstractFloat,2},
} # 2D2F flux

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
    fw[2] += Mt[4] * sum(ω .* u .^ 2 .* h) - Mt[5] * sum(ω .* u .^ 3 .* sh)
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

    return nothing

end

# ------------------------------------------------------------
# 3F2V with AAP model
# ------------------------------------------------------------
function flux_ugks!(
    fw::T1,
    fh0::T2,
    fh1::T2,
    fh2::T2,
    wL::T3,
    h0L::T4,
    h1L::T4,
    h2L::T4,
    wR::T3,
    h0R::T4,
    h1R::T4,
    h2R::T4,
    u::T5,
    v::T5,
    ω::T5,
    inK,
    γ,
    mi,
    ni,
    me,
    ne,
    Kn,
    dt,
    dxL,
    dxR,
    len,
    sh0L = zeros(eltype(h0L), axes(h0L))::T4,
    sh1L = zeros(eltype(h1L), axes(h1L))::T4,
    sh2L = zeros(eltype(h2L), axes(h2L))::T4,
    sh0R = zeros(eltype(h0R), axes(h0R))::T4,
    sh1R = zeros(eltype(h1R), axes(h1R))::T4,
    sh2R = zeros(eltype(h2R), axes(h2R))::T4,
) where {
    T1<:AbstractArray{<:AbstractFloat,2},
    T2<:AbstractArray{<:AbstractFloat,3},
    T3<:AbstractArray{<:Real,2},
    T4<:AbstractArray{<:AbstractFloat,3},
    T5<:AbstractArray{<:AbstractFloat,3},
}

    #--- reconstruct initial distribution ---#
    δ = heaviside.(u)

    h0 = @. h0L * δ + h0R * (1.0 - δ)
    h1 = @. h1L * δ + h1R * (1.0 - δ)
    h2 = @. h2L * δ + h2R * (1.0 - δ)

    sh0 = @. sh0L * δ + sh0R * (1.0 - δ)
    sh1 = @. sh1L * δ + sh1R * (1.0 - δ)
    sh2 = @. sh2L * δ + sh2R * (1.0 - δ)

    primL = mixture_conserve_prim(wL, γ)
    primR = mixture_conserve_prim(wR, γ)

    #--- construct interface distribution ---#
    Mu1, Mv1, Mxi1, MuL1, MuR1 = mixture_gauss_moments(primL, inK)
    Muv1 = mixture_moments_conserve(MuL1, Mv1, Mxi1, 0, 0, 0)
    Mu2, Mv2, Mxi2, MuL2, MuR2 = mixture_gauss_moments(primR, inK)
    Muv2 = mixture_moments_conserve(MuR2, Mv2, Mxi2, 0, 0, 0)

    w = similar(wL)
    for j in axes(w, 2)
        @. w[:, j] = primL[1, j] * Muv1[:, j] + primR[1, j] * Muv2[:, j]
    end
    prim = mixture_conserve_prim(w, γ)
    tau = aap_hs_collision_time(prim, mi, ni, me, ne, Kn)
    prim = aap_hs_prim(prim, tau, mi, ni, me, ne, Kn)

    a = mixture_pdf_slope(prim, (wR .- wL) ./ (dxR + dxL), inK)
    Mu, Mv, Mxi, MuL, MuR = mixture_gauss_moments(prim, inK)
    Mau = mixture_moments_conserve_slope(a, Mu, Mv, Mxi, 1, 0, 0)
    aT = mixture_pdf_slope(prim, -prim[1] .* Mau, inK)

    Mt = zeros(2, 2)
    @. Mt[2, :] = tau * (1.0 - exp(-dt / tau)) # f0
    @. Mt[1, :] = dt - Mt[2, :] # M0

    Mt = zeros(5, axes(w, 2))
    for j in axes(Mt, 2)
        Mt[4, j] = tau[j] * (1.0 - exp(-dt / tau[j])) # f0
        Mt[5, j] = -tau[j] * dt * exp(-dt / tau[j]) + tau[j] * Mt[4, j]
        Mt[1, j] = dt - Mt[4, j] # M0
        Mt[2, j] = -tau[j] * Mt[1, j] + Mt[5, j]
        Mt[3, j] = dt^2 / 2.0 - tau[j] * Mt[1, j]
    end

    #--- calculate interface flux ---#
    # flux from M0
    Muv = mixture_moments_conserve(Mu, Mv, Mxi, 1, 0, 0)
    Mau = mixture_moments_conserve_slope(a, Mu, Mv, Mxi, 2, 0, 0)
    MauT = mixture_moments_conserve_slope(aT, Mu, Mv, Mxi, 1, 0, 0)
    for j in axes(fw, 2)
        @. fw[:, j] =
            Mt[1, j] * prim[1, j] * Muv[:, j] +
            Mt[2, j] * prim[1, j] * Mau[:, j] +
            Mt[3, j] * prim[1, j] * MauT[:, j]
    end

    # flux from f0
    H0 = mixture_maxwellian(u, v, prim)
    H1 = similar(H0)
    H2 = similar(H0)
    for j in axes(H0, 3)
        H1[:, :, j] = H0[:, :, j] .* prim[4, j]
        H2[:, :, j] .= H0[:, :, j] .* (prim[4, j]^2 + 1.0 / (2.0 * prim[5, j]))
    end

    for j in axes(fw, 2)
        fw[1, j] +=
            Mt[4, j] * sum(ω[:, :, j] .* u[:, :, j] .* h0[:, :, j]) -
            Mt[5, j] * sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* sh0[:, :, j])
        fw[2, j] +=
            Mt[4, j] * sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* h0[:, :, j]) -
            Mt[5, j] * sum(ω[:, :, j] .* u[:, :, j] .^ 3 .* sh0[:, :, j])
        fw[3, j] +=
            Mt[4, j] *
            sum(ω[:, :, j] .* v[:, :, j] .* u[:, :, j] .* h0[:, :, j]) -
            Mt[5, j] *
            sum(ω[:, :, j] .* v[:, :, j] .* u[:, :, j] .^ 2 .* sh0[:, :, j])
        fw[4, j] +=
            Mt[4, j] * sum(ω[:, :, j] .* u[:, :, j] .* h1[:, :, j]) -
            Mt[5, j] * sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* sh1[:, :, j])
        fw[5, j] +=
            Mt[4, j] *
            0.5 *
            (
                sum(
                    ω[:, :, j] .* u[:, :, j] .*
                    (u[:, :, j] .^ 2 .+ v[:, :, j] .^ 2) .* h0[:, :, j],
                ) + sum(ω[:, :, j] .* u[:, :, j] .* h2[:, :, j])
            ) -
            Mt[5, j] *
            0.5 *
            (
                sum(
                    ω[:, :, j] .* u[:, :, j] .^ 2 .*
                    (u[:, :, j] .^ 2 .+ v[:, :, j] .^ 2) .* sh0[:, :, j],
                ) + sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* sh2[:, :, j])
            )

        @. fh0[:, :, j] =
            Mt[1, j] * u[:, :, j] * H0[:, :, j] +
            Mt[2, j] *
            u[:, :, j]^2 *
            (
                a[1, j] * H0[:, :, j] +
                a[2, j] * u[:, :, j] * H0[:, :, j] +
                a[3, j] * v[:, :, j] * H0[:, :, j] +
                a[4, j] * u[:, :, j] * H1[:, :, j] +
                0.5 *
                a[5, j] *
                ((u[:, :, j]^2 + v[:, :, j]^2) * H0[:, :, j] + H2[:, :, j])
            ) +
            Mt[3, j] *
            u[:, :, j] *
            (
                aT[1, j] * H0[:, :, j] +
                aT[2, j] * u[:, :, j] * H0[:, :, j] +
                aT[3, j] * v[:, :, j] * H0[:, :, j] +
                aT[4, j] * u[:, :, j] * H1[:, :, j] +
                0.5 *
                aT[5, j] *
                ((u[:, :, j]^2 + v[:, :, j]^2) * H0[:, :, j] + H2[:, :, j])
            ) +
            Mt[4, j] * u[:, :, j] * h0[:, :, j] -
            Mt[5, j] * u[:, :, j]^2 * sh0[:, :, j]
        @. fh1[:, :, j] =
            Mt[1, j] * u[:, :, j] * H1[:, :, j] +
            Mt[2, j] *
            u[:, :, j]^2 *
            (
                a[1, j] * H1[:, :, j] +
                a[2, j] * u[:, :, j] * H1[:, :, j] +
                a[3, j] * v[:, :, j] * H1[:, :, j] +
                a[4, j] * u[:, :, j] * H2[:, :, j] +
                0.5 *
                a[5, j] *
                (
                    (u[:, :, j]^2 + v[:, :, j]^2) * H1[:, :, j] +
                    Mxi[3, j] * H0[:, :, j]
                )
            ) +
            Mt[3, j] *
            u[:, :, j] *
            (
                aT[1, j] * H1[:, :, j] +
                aT[2, j] * u[:, :, j] * H1[:, :, j] +
                aT[3, j] * v[:, :, j] * H1[:, :, j] +
                aT[4, j] * u[:, :, j] * H2[:, :, j] +
                0.5 *
                aT[5, j] *
                (
                    (u[:, :, j]^2 + v[:, :, j]^2) * H1[:, :, j] +
                    Mxi[3, j] * H0[:, :, j]
                )
            ) +
            Mt[4, j] * u[:, :, j] * h1[:, :, j] -
            Mt[5, j] * u[:, :, j]^2 * sh1[:, :, j]
        @. fh2[:, :, j] =
            Mt[1, j] * u[:, :, j] * H2[:, :, j] +
            Mt[2, j] *
            u[:, :, j]^2 *
            (
                a[1, j] * H2[:, :, j] +
                a[2, j] * u[:, :, j] * H2[:, :, j] +
                a[3, j] * v[:, :, j] * H2[:, :, j] +
                a[4, j] * u[:, :, j] * Mxi[3, j] * H0[:, :, j] +
                0.5 *
                a[5, j] *
                (
                    (u[:, :, j]^2 + v[:, :, j]^2) * H2[:, :, j] +
                    Mxi[4, j] * H0[:, :, j]
                )
            ) +
            Mt[3, j] *
            u[:, :, j] *
            (
                aT[1, j] * H2[:, :, j] +
                aT[2, j] * u[:, :, j] * H2[:, :, j] +
                aT[3, j] * v[:, :, j] * H2[:, :, j] +
                aT[4, j] * u[:, :, j] * Mxi[3, j] * H0[:, :, j] +
                0.5 *
                aT[5, j] *
                (
                    (u[:, :, j]^2 + v[:, :, j]^2) * H2[:, :, j] +
                    Mxi[4, j] * H0[:, :, j]
                )
            ) +
            Mt[4, j] * u[:, :, j] * h2[:, :, j] -
            Mt[5, j] * u[:, :, j]^2 * sh2[:, :, j]
    end

    @. fw *= len
    @. fh0 *= len
    @. fh1 *= len
    @. fh2 *= len

    return nothing

end
