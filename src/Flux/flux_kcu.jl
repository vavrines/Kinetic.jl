"""
Kinetic central-upwind (KCU) method

    1D1F1V: flux_kcu!(fw, ff, wL, fL, wR, fR, u, ω, inK, γ, visRef, visIdx, Pr, dt)
    1D2F1V: flux_kcu!(fw, fh, fb, wL, hL, bL, wR, hR, bR, u, ω, inK, γ, visRef, visIdx, Pr, dt)
    1D4F1V: flux_kcu!(fw, fh0, fh1, fh2, fh3, wL, h0L, h1L, h2L, h3L, wR, h0R, h1R, h2R, h3R, u, ω, inK, γ, visRef, visIdx, Pr, dt)
    2D1F2V: flux_kcu!(fw, ff, wL, fL, wR, fR, u, v, ω, inK, γ, visRef, visIdx, Pr, dt, len)
    2D2F2V: flux_kcu!(fw, fh, fb, wL, hL, bL, wR, hR, bR, u, v, ω, inK, γ, visRef, visIdx, Pr, dt, len)
    2D3F2V: flux_kcu!(fw, fh0, fh1, fh2, wL, h0L, h1L, h2L, wR, h0R, h1R, h2R, u, v, ω, inK, γ, visRef, visIdx, Pr, dt, len)

- @args: particle distribution functions and their slopes at left/right sides of interface
- @args: particle velocity quadrature points and weights
- @args: time step and cell size

"""
function flux_kcu!(
    fw::X,
    ff::Y,
    wL::Z,
    fL::A,
    wR::Z,
    fR::A,
    u::B,
    ω::B,
    inK,
    γ,
    visRef,
    visIdx,
    Pr,
    dt,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,1},
    Z<:AbstractArray{<:Real,1},
    A<:AbstractArray{<:AbstractFloat,1},
    B<:AbstractArray{<:AbstractFloat,1},
} # 1D1F1V

    # upwind reconstruction
    δ = heaviside.(u)
    f = @. fL * δ + fR * (1.0 - δ)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    # construct interface distribution
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

    # calculate fluxes
    Mu, Mxi, MuL, MuR = gauss_moments(prim, inK)

    ## flux from M0
    Muv = moments_conserve(Mu, Mxi, 1, 0)
    @. fw = Mt[1] * prim[1] * Muv

    ## flux from f0
    g = maxwellian(u, prim)

    fw[1] += Mt[2] * sum(ω .* u .* f)
    fw[2] += Mt[2] * sum(ω .* u .^ 2 .* f)
    fw[3] += Mt[2] * 0.5 * (sum(ω .* u .^ 3 .* f))

    @. ff = Mt[1] * u * g + Mt[2] * u * f

    return nothing

end

#--- mixture ---#
function flux_kcu!(
    fw::X,
    ff::Y,
    wL::Z,
    fL::A,
    wR::Z,
    fR::A,
    u::B,
    ω::B,
    inK,
    γ,
    mi,
    ni,
    me,
    ne,
    Kn,
    dt,
) where {
    X<:AbstractArray{<:AbstractFloat,2},
    Y<:AbstractArray{<:AbstractFloat,2},
    Z<:AbstractArray{<:Real,2},
    A<:AbstractArray{<:AbstractFloat,2},
    B<:AbstractArray{<:AbstractFloat,2},
}

    # upwind reconstruction
    δ = heaviside.(u)
    f = @. fL * δ + fR * (1.0 - δ)

    primL = mixture_conserve_prim(wL, γ)
    primR = mixture_conserve_prim(wR, γ)

    # construct interface distribution
    Mu1, Mxi1, MuL1, MuR1 = mixture_gauss_moments(primL, inK)
    Mu2, Mxi2, MuL2, MuR2 = mixture_gauss_moments(primR, inK)
    Muv1 = mixture_moments_conserve(MuL1, Mxi1, 0, 0)
    Muv2 = mixture_moments_conserve(MuR2, Mxi2, 0, 0)

    w = similar(wL)
    for j in axes(w, 2)
        @. w[:, j] = primL[1, j] * Muv1[:, j] + primR[1, j] * Muv2[:, j]
    end
    prim = mixture_conserve_prim(w, γ)

    tau = aap_hs_collision_time(prim, mi, ni, me, ne, Kn)
    #@. tau +=
    #    abs(cellL.prim[1, :] / cellL.prim[end, :] - cellR.prim[1, :] / cellR.prim[end, :]) /
    #    (cellL.prim[1, :] / cellL.prim[end, :] + cellR.prim[1, :] / cellR.prim[end, :]) *
    #    dt *
    #    2.0
    prim = aap_hs_prim(prim, tau, mi, ni, me, ne, Kn)

    Mt = zeros(2, 2)
    @. Mt[2, :] = tau * (1.0 - exp(-dt / tau)) # f0
    @. Mt[1, :] = dt - Mt[2, :] # M0

    # calculate fluxes
    Mu, Mxi, MuL, MuR = mixture_gauss_moments(prim, inK)
    Muv = mixture_moments_conserve(Mu, Mxi, 1, 0)

    ## flux from M0
    for j in axes(fw, 2)
        @. fw[:, j] = Mt[1, j] * prim[1, j] * Muv[:, j]
    end

    ## flux from f0
    M = mixture_maxwellian(u, prim)

    for j in axes(fw, 2)
        fw[1, j] += Mt[2, j] * sum(ω[:, j] .* u[:, j] .* f[:, j])
        fw[2, j] += Mt[2, j] * sum(ω[:, j] .* u[:, j] .^ 2 .* f[:, j])
        fw[3, j] += Mt[2, j] * 0.5 * sum(ω[:, j] .* u[:, j] .^ 3 .* f[:, j])

        @. ff[:, j] =
            Mt[1, j] * u[:, j] * M[:, j] + Mt[2, j] * u[:, j] * f[:, j]
    end

    return nothing

end

# ------------------------------------------------------------
# 1D2F1V
# ------------------------------------------------------------
function flux_kcu!(
    fw::X,
    fh::Y,
    fb::Y,
    wL::Z,
    hL::A,
    bL::A,
    wR::Z,
    hR::A,
    bR::A,
    u::B,
    ω::B,
    inK,
    γ,
    visRef,
    visIdx,
    Pr,
    dt,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,1},
    Z<:AbstractArray{<:Real,1},
    A<:AbstractArray{<:AbstractFloat,1},
    B<:AbstractArray{<:AbstractFloat,1},
}

    # upwind reconstruction
    δ = heaviside.(u)
    h = @. hL * δ + hR * (1.0 - δ)
    b = @. bL * δ + bR * (1.0 - δ)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    # construct interface distribution
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

    # calculate fluxes
    Mu, Mxi, MuL, MuR = gauss_moments(prim, inK)

    ## flux from M0
    Muv = moments_conserve(Mu, Mxi, 1, 0)
    @. fw = Mt[1] * prim[1] * Muv

    ## flux from f0
    Mh = maxwellian(u, prim)
    Mb = Mh .* inK ./ (2.0 * prim[end])

    fw[1] += Mt[2] * sum(ω .* u .* h)
    fw[2] += Mt[2] * sum(ω .* u .^ 2 .* h)
    fw[3] += Mt[2] * 0.5 * (sum(ω .* u .^ 3 .* h) + sum(ω .* u .* b))

    @. fh = Mt[1] * u * Mh + Mt[2] * u * h
    @. fb = Mt[1] * u * Mb + Mt[2] * u * b

    return nothing

end

#--- mixture ---#
function flux_kcu!(
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
    mi,
    ni,
    me,
    ne,
    Kn,
    dt,
) where {
    T1<:AbstractArray{<:AbstractFloat,2},
    T2<:AbstractArray{<:AbstractFloat,2},
    T3<:AbstractArray{<:Real,2},
    T4<:AbstractArray{<:AbstractFloat,2},
    T5<:AbstractArray{<:AbstractFloat,2},
}

    # upwind reconstruction
    δ = heaviside.(u)

    h = @. hL * δ + hR * (1.0 - δ)
    b = @. bL * δ + bR * (1.0 - δ)

    primL = mixture_conserve_prim(wL, γ)
    primR = mixture_conserve_prim(wR, γ)

    # construct interface distribution
    Mu1, Mxi1, MuL1, MuR1 = mixture_gauss_moments(primL, inK)
    Mu2, Mxi2, MuL2, MuR2 = mixture_gauss_moments(primR, inK)
    Muv1 = mixture_moments_conserve(MuL1, Mxi1, 0, 0)
    Muv2 = mixture_moments_conserve(MuR2, Mxi2, 0, 0)

    w = similar(wL)
    for j in axes(w, 2)
        @. w[:, j] = primL[1, j] * Muv1[:, j] + primR[1, j] * Muv2[:, j]
    end
    prim = mixture_conserve_prim(w, γ)

    tau = aap_hs_collision_time(prim, mi, ni, me, ne, Kn)
    #@. tau +=
    #    abs(cellL.prim[1, :] / cellL.prim[end, :] - cellR.prim[1, :] / cellR.prim[end, :]) /
    #    (cellL.prim[1, :] / cellL.prim[end, :] + cellR.prim[1, :] / cellR.prim[end, :]) *
    #    dt *
    #    2.0
    prim = aap_hs_prim(prim, tau, mi, ni, me, ne, Kn)

    Mt = zeros(2, 2)
    @. Mt[2, :] = tau * (1.0 - exp(-dt / tau)) # f0
    @. Mt[1, :] = dt - Mt[2, :] # M0

    # calculate fluxes
    Mu, Mxi, MuL, MuR = mixture_gauss_moments(prim, inK)
    Muv = mixture_moments_conserve(Mu, Mxi, 1, 0)

    ## flux from M0
    for j in axes(fw, 2)
        @. fw[:, j] = Mt[1, j] * prim[1, j] * Muv[:, j]
    end

    ## flux from f0
    MH = mixture_maxwellian(u, prim)
    MB = similar(MH)
    for j in axes(MB, 2)
        MB[:, j] .= MH[:, j] .* inK ./ (2.0 * prim[end, j])
    end

    for j in axes(fw, 2)
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
            Mt[1, j] * u[:, j] * MH[:, j] + Mt[2, j] * u[:, j] * h[:, j]
        @. fb[:, j] =
            Mt[1, j] * u[:, j] * MB[:, j] + Mt[2, j] * u[:, j] * b[:, j]
    end

    return nothing

end

# ------------------------------------------------------------
# 1D4F1V
# ------------------------------------------------------------
function flux_kcu!(
    fw::X,
    fh0::Y,
    fh1::Y,
    fh2::Y,
    fh3::Y,
    wL::Z,
    h0L::A,
    h1L::A,
    h2L::A,
    h3L::A,
    wR::Z,
    h0R::A,
    h1R::A,
    h2R::A,
    h3R::A,
    u::B,
    ω::B,
    inK,
    γ,
    visRef,
    visIdx,
    Pr,
    dt,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,1},
    Z<:AbstractArray{<:Real,1},
    A<:AbstractArray{<:AbstractFloat,1},
    B<:AbstractArray{<:AbstractFloat,1},
}

    # upwind reconstruction
    δ = heaviside.(u)

    h0 = @. h0L * δ + h0R * (1.0 - δ)
    h1 = @. h1L * δ + h1R * (1.0 - δ)
    h2 = @. h2L * δ + h2R * (1.0 - δ)
    h3 = @. h3L * δ + h3R * (1.0 - δ)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    # construct interface distribution
    Mu1, Mv1, Mw1, MuL1, MuR1 = gauss_moments(primL, inK)
    Muv1 = moments_conserve(MuL1, Mv1, Mw1, 0, 0, 0)
    Mu2, Mv2, Mw2, MuL2, MuR2 = gauss_moments(primR, inK)
    Muv2 = moments_conserve(MuR2, Mv2, Mw2, 0, 0, 0)

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

    # calculate fluxes
    Mu, Mv, Mw, MuL, MuR = gauss_moments(prim, inK)

    ## flux from M0
    Muv = moments_conserve(Mu, Mv, Mw, 1, 0, 0)
    @. fw = Mt[1] * prim[1] * Muv

    ## flux from f0
    g0 = maxwellian(u, prim)
    g1 = Mv[1] .* g0
    g2 = Mw[1] .* g0
    g3 = (Mv[2] + Mw[2]) .* g0

    fw[1] += Mt[2] * sum(ω .* u .* h0)
    fw[2] += Mt[2] * sum(ω .* u .^ 2 .* h0)
    fw[3] += Mt[2] * sum(ω .* u .* h1)
    fw[4] += Mt[2] * sum(ω .* u .* h2)
    fw[5] += Mt[2] * 0.5 * (sum(ω .* u .^ 3 .* h0) + sum(ω .* u .* h3))

    @. fh0 = Mt[1] * u * g0 + Mt[2] * u * h0
    @. fh1 = Mt[1] * u * g1 + Mt[2] * u * h1
    @. fh2 = Mt[1] * u * g2 + Mt[2] * u * h2
    @. fh3 = Mt[1] * u * g3 + Mt[2] * u * h3

    return nothing

end

#--- mixture ---#
function flux_kcu!(
    fw::X,
    fh0::Y,
    fh1::Y,
    fh2::Y,
    fh3::Y,
    wL::Z,
    h0L::A,
    h1L::A,
    h2L::A,
    h3L::A,
    wR::Z,
    h0R::A,
    h1R::A,
    h2R::A,
    h3R::A,
    u::B,
    ω::B,
    inK,
    γ,
    mi,
    ni,
    me,
    ne,
    Kn,
    dt,
    isMHD = false::Bool,
) where {
    X<:AbstractArray{<:AbstractFloat,2},
    Y<:AbstractArray{<:AbstractFloat,2},
    Z<:AbstractArray{<:Real,2},
    A<:AbstractArray{<:AbstractFloat,2},
    B<:AbstractArray{<:AbstractFloat,2},
}

    #--- upwind reconstruction ---#
    δ = heaviside.(u)

    h0 = @. h0L * δ + h0R * (1.0 - δ)
    h1 = @. h1L * δ + h1R * (1.0 - δ)
    h2 = @. h2L * δ + h2R * (1.0 - δ)
    h3 = @. h3L * δ + h3R * (1.0 - δ)

    primL = mixture_conserve_prim(wL, γ)
    primR = mixture_conserve_prim(wR, γ)

    # construct interface distribution
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
    #@. tau +=
    #    abs(primL[1, :] / primL[end, :] - primR[1, :] / primR[end, :]) /
    #    (primL[1, :] / primL[end, :] + primR[1, :] / primR[end, :]) *
    #    dt *
    #    5.0
    if !isMHD
        prim = aap_hs_prim(prim, tau, mi, ni, me, ne, Kn)
    end

    Mt = zeros(2, 2)
    @. Mt[2, :] = tau * (1.0 - exp(-dt / tau)) # f0
    @. Mt[1, :] = dt - Mt[2, :] # M0

    # calculate fluxes
    Mu, Mv, Mw, MuL, MuR = mixture_gauss_moments(prim, inK)
    Muv = mixture_moments_conserve(Mu, Mv, Mw, 1, 0, 0)

    ## flux from M0
    for j in axes(fw, 2)
        @. fw[:, j] = Mt[1, j] * prim[1, j] * Muv[:, j]
    end

    ## flux from f0
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

        @. fh0[:, j] =
            Mt[1, j] * u[:, j] * g0[:, j] + Mt[2, j] * u[:, j] * h0[:, j]
        @. fh1[:, j] =
            Mt[1, j] * u[:, j] * g1[:, j] + Mt[2, j] * u[:, j] * h1[:, j]
        @. fh2[:, j] =
            Mt[1, j] * u[:, j] * g2[:, j] + Mt[2, j] * u[:, j] * h2[:, j]
        @. fh3[:, j] =
            Mt[1, j] * u[:, j] * g3[:, j] + Mt[2, j] * u[:, j] * h3[:, j]
    end

    return nothing

end

# ------------------------------------------------------------
# 2D1F2V
# ------------------------------------------------------------
function flux_kcu!(
    fw::X,
    ff::Y,
    wL::Z,
    fL::A,
    wR::Z,
    fR::A,
    u::B,
    v::B,
    ω::B,
    inK,
    γ,
    visRef,
    visIdx,
    Pr,
    dt,
    len,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,2},
    Z<:AbstractArray{<:Real,1},
    A<:AbstractArray{<:AbstractFloat,2},
    B<:AbstractArray{<:AbstractFloat,2},
}

    # reconstruct initial distribution
    δ = heaviside.(u)
    f = @. fL * δ + fR * (1.0 - δ)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    # construct interface distribution
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

    # calculate interface flux
    Mu, Mv, Mxi, MuL, MuR = gauss_moments(prim, inK)

    ## flux from M0
    Muv = moments_conserve(Mu, Mv, Mxi, 1, 0, 0)
    @. fw = Mt[1] * prim[1] * Muv * len

    ## flux from f0
    g = maxwellian(u, v, prim)

    fw[1] += Mt[2] * sum(ω .* u .* f) * len
    fw[2] += Mt[2] * sum(ω .* u .^ 2 .* f) * len
    fw[3] += Mt[2] * sum(ω .* v .* u .* f) * len
    fw[4] += Mt[2] * 0.5 * (sum(ω .* u .* (u .^ 2 .+ v .^ 2) .* f)) * len

    @. ff = (Mt[1] * u * g + Mt[2] * u * f) * len

    return nothing

end

# ------------------------------------------------------------
# 2D2F2V
# ------------------------------------------------------------
function flux_kcu!(
    fw::T1,
    fh::T2,
    fb::T2,
    wL::T3,
    hL::T4,
    bL::T4,
    wR::T3,
    hR::T1,
    bR::T1,
    u::T5,
    v::T5,
    ω::T5,
    inK,
    γ,
    visRef,
    visIdx,
    Pr,
    dt,
    len,
) where {
    T1<:AbstractArray{<:AbstractFloat,1},
    T2<:AbstractArray{<:AbstractFloat,2},
    T3<:AbstractArray{<:Real,1},
    T4<:AbstractArray{<:AbstractFloat,2},
    T5<:AbstractArray{<:AbstractFloat,2},
}

    # reconstruct initial distribution
    δ = heaviside.(u)
    h = @. hL * δ + hR * (1.0 - δ)
    b = @. bL * δ + bR * (1.0 - δ)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    # construct interface distribution
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

    # calculate interface flux
    Mu, Mv, Mxi, MuL, MuR = gauss_moments(prim, inK)

    ## flux from M0
    Muv = moments_conserve(Mu, Mv, Mxi, 1, 0, 0)
    @. fw = Mt[1] * prim[1] * Muv * len

    ## flux from f0
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

    return nothing

end

# ------------------------------------------------------------
# 2D3F2V
# ------------------------------------------------------------
function flux_kcu!(
    fw::X,
    fh0::Y,
    fh1::Y,
    fh2::Y,
    wL::Z,
    h0L::A,
    h1L::A,
    h2L::A,
    wR::Z,
    h0R::A,
    h1R::A,
    h2R::A,
    u::B,
    v::B,
    ω::B,
    inK,
    γ,
    visRef,
    visIdx,
    Pr,
    dt,
    len,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,2},
    Z<:AbstractArray{<:Real,1},
    A<:AbstractArray{<:AbstractFloat,2},
    B<:AbstractArray{<:AbstractFloat,2},
}

    # reconstruct initial distribution
    δ = heaviside.(u)
    h0 = @. h0L * δ + h0R * (1.0 - δ)
    h1 = @. h1L * δ + h1R * (1.0 - δ)
    h2 = @. h2L * δ + h2R * (1.0 - δ)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    # construct interface distribution
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

    # calculate interface flux
    Mu, Mv, Mxi, MuL, MuR = gauss_moments(prim, inK)

    ## flux from M0
    Muv = moments_conserve(Mu, Mv, Mxi, 1, 0, 0)
    @. fw = Mt[1] * prim[1] * Muv

    ## flux from f0
    H0 = maxwellian(u, v, prim)
    H1 = H0 .* prim[4]
    H2 = H0 .* (prim[4]^2 + 1.0 / (2.0 * prim[5]))

    fw[1] += Mt[2] * sum(ω .* u .* h0)
    fw[2] += Mt[2] * sum(ω .* u .^ 2 .* h0)
    fw[3] += Mt[2] * sum(ω .* v .* u .* h0)
    fw[4] += Mt[2] * sum(ω .* u .* h1)
    fw[5] +=
        Mt[2] *
        0.5 *
        (sum(ω .* u .* (u .^ 2 .+ v .^ 2) .* h0) + sum(ω .* u .* h2))

    @. fw *= len
    @. fh0 = (Mt[1] * u * H0 + Mt[2] * u * h0) * len
    @. fh1 = (Mt[1] * u * H1 + Mt[2] * u * h1) * len
    @. fh2 = (Mt[1] * u * H2 + Mt[2] * u * h2) * len

    return nothing

end

#--- mixture ---#
function flux_kcu!(
    fw::X,
    fh0::Y,
    fh1::Y,
    fh2::Y,
    wL::Z,
    h0L::A,
    h1L::A,
    h2L::A,
    wR::Z,
    h0R::A,
    h1R::A,
    h2R::A,
    u::B,
    v::B,
    ω::B,
    inK,
    γ,
    mi,
    ni,
    me,
    ne,
    Kn,
    dt,
    len,
    isMHD = false::Bool,
) where {
    X<:AbstractArray{<:AbstractFloat,2},
    Y<:AbstractArray{<:AbstractFloat,3},
    Z<:AbstractArray{<:Real,2},
    A<:AbstractArray{<:AbstractFloat,3},
    B<:AbstractArray{<:AbstractFloat,3},
}

    # reconstruct initial distribution
    δ = heaviside.(u)
    h0 = @. h0L * δ + h0R * (1.0 - δ)
    h1 = @. h1L * δ + h1R * (1.0 - δ)
    h2 = @. h2L * δ + h2R * (1.0 - δ)

    primL = mixture_conserve_prim(wL, γ)
    primR = mixture_conserve_prim(wR, γ)

    # construct interface distribution
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
    #@. tau +=
    #    abs(primL[1, :] / primL[end, :] - primR[1, :] / primR[end, :]) /
    #    (primL[1, :] / primL[end, :] + primR[1, :] / primR[end, :]) *
    #    dt *
    #    2.0
    if !isMHD
        prim = aap_hs_prim(prim, tau, mi, ni, me, ne, Kn)
    end

    Mt = zeros(2, 2)
    @. Mt[2, :] = tau * (1.0 - exp(-dt / tau)) # f0
    @. Mt[1, :] = dt - Mt[2, :] # M0

    # calculate interface flux
    Mu, Mv, Mxi, MuL, MuR = mixture_gauss_moments(prim, inK)

    ## flux from M0
    Muv = mixture_moments_conserve(Mu, Mv, Mxi, 1, 0, 0)
    for j in axes(fw, 2)
        @. fw[:, j] = Mt[1, j] * prim[1, j] * Muv[:, j]
    end

    ## flux from f0
    H0 = mixture_maxwellian(u, v, prim)
    H1 = similar(H0)
    H2 = similar(H0)
    for j in axes(H0, 3)
        H1[:, :, j] = H0[:, :, j] .* prim[4, j]
        H2[:, :, j] .= H0[:, :, j] .* (prim[4, j]^2 + 1.0 / (2.0 * prim[5, j]))
    end

    for j in axes(fw, 2)
        fw[1, j] += Mt[2, j] * sum(ω[:, :, j] .* u[:, :, j] .* h0[:, :, j])
        fw[2, j] += Mt[2, j] * sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* h0[:, :, j])
        fw[3, j] +=
            Mt[2, j] *
            sum(ω[:, :, j] .* v[:, :, j] .* u[:, :, j] .* h0[:, :, j])
        fw[4, j] += Mt[2, j] * sum(ω[:, :, j] .* u[:, :, j] .* h1[:, :, j])
        fw[5, j] +=
            Mt[2, j] *
            0.5 *
            (
                sum(
                    ω[:, :, j] .* u[:, :, j] .*
                    (u[:, :, j] .^ 2 .+ v[:, :, j] .^ 2) .* h0[:, :, j],
                ) + sum(ω[:, :, j] .* u[:, :, j] .* h2[:, :, j])
            )

        @. fh0[:, :, j] =
            Mt[1, j] * u[:, :, j] * H0[:, :, j] +
            Mt[2, j] * u[:, :, j] * h0[:, :, j]
        @. fh1[:, :, j] =
            Mt[1, j] * u[:, :, j] * H1[:, :, j] +
            Mt[2, j] * u[:, :, j] * h1[:, :, j]
        @. fh2[:, :, j] =
            Mt[1, j] * u[:, :, j] * H2[:, :, j] +
            Mt[2, j] * u[:, :, j] * h2[:, :, j]
    end

    @. fw *= len
    @. fh0 *= len
    @. fh1 *= len
    @. fh2 *= len

    return nothing

end
