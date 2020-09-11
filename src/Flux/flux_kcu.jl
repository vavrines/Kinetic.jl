"""
Kinetic central-upwind (KCU) method

* @arg: particle distribution functions and their slopes at left/right sides of interface
* @arg: particle velocity quadrature points and weights
* @arg: time step and cell size

"""
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
) # 1D1F flux

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
    fw[4] += Mt[2] * 0.5 * (sum(ω .* u .* (u .^ 2 .+ v .^ 2) .* h) + sum(ω .* u .* b)) * len

    @. fh = (Mt[1] * u * H + Mt[2] * u * h) * len
    @. fb = (Mt[1] * u * B + Mt[2] * u * b) * len

end


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
        Mu1[:, j], Mxi1[:, j], MuL1[:, j], MuR1[:, j] = gauss_moments(primL[:, j], inK)
        Muv1[:, j] = moments_conserve(MuL1[:, j], Mxi1[:, j], 0, 0)
        Mu2[:, j], Mxi2[:, j], MuL2[:, j], MuR2[:, j] = gauss_moments(primR[:, j], inK)
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
        abs(cellL.prim[1, :] / cellL.prim[end, :] - cellR.prim[1, :] / cellR.prim[end, :]) /
        (cellL.prim[1, :] / cellL.prim[end, :] + cellR.prim[1, :] / cellR.prim[end, :]) *
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
        Mu[:, j], Mxi[:, j], MuL[:, j], MuR[:, j] = gauss_moments(prim[:, j], inK)
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

        @. ff[:, j] = Mt[1, j] * u[:, j] * g[:, j] + Mt[2, j] * u[:, j] * f[:, j]
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
        Mu1[:, j], Mxi1[:, j], MuL1[:, j], MuR1[:, j] = gauss_moments(primL[:, j], inK)
        Muv1[:, j] = moments_conserve(MuL1[:, j], Mxi1[:, j], 0, 0)
        Mu2[:, j], Mxi2[:, j], MuL2[:, j], MuR2[:, j] = gauss_moments(primR[:, j], inK)
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
        abs(cellL.prim[1, :] / cellL.prim[end, :] - cellR.prim[1, :] / cellR.prim[end, :]) /
        (cellL.prim[1, :] / cellL.prim[end, :] + cellR.prim[1, :] / cellR.prim[end, :]) *
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
        Mu[:, j], Mxi[:, j], MuL[:, j], MuR[:, j] = gauss_moments(prim[:, j], inK)
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
            (sum(ω[:, j] .* u[:, j] .^ 3 .* h[:, j]) + sum(ω[:, j] .* u[:, j] .* b[:, j]))

        @. fh[:, j] = Mt[1, j] * u[:, j] * g0[:, j] + Mt[2, j] * u[:, j] * h[:, j]
        @. fb[:, j] = Mt[1, j] * u[:, j] * g1[:, j] + Mt[2, j] * u[:, j] * b[:, j]
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

        @. fh0[:, j] = Mt[1, j] * u[:, j] * g0[:, j] + Mt[2, j] * u[:, j] * h0[:, j]
        @. fh1[:, j] = Mt[1, j] * u[:, j] * g1[:, j] + Mt[2, j] * u[:, j] * h1[:, j]
        @. fh2[:, j] = Mt[1, j] * u[:, j] * g2[:, j] + Mt[2, j] * u[:, j] * h2[:, j]
        @. fh3[:, j] = Mt[1, j] * u[:, j] * g3[:, j] + Mt[2, j] * u[:, j] * h3[:, j]
    end

end


# ------------------------------------------------------------
# 2D3F flux with AAP model
# ------------------------------------------------------------
function flux_kcu!(
    fw::AbstractArray{<:AbstractFloat,2},
    fh0::AbstractArray{<:AbstractFloat,3},
    fh1::AbstractArray{<:AbstractFloat,3},
    fh2::AbstractArray{<:AbstractFloat,3},
    wL::AbstractArray{<:Real,2},
    h0L::AbstractArray{<:AbstractFloat,3},
    h1L::AbstractArray{<:AbstractFloat,3},
    h2L::AbstractArray{<:AbstractFloat,3},
    wR::AbstractArray{<:Real,2},
    h0R::AbstractArray{<:AbstractFloat,3},
    h1R::AbstractArray{<:AbstractFloat,3},
    h2R::AbstractArray{<:AbstractFloat,3},
    u::AbstractArray{<:AbstractFloat,3},
    v::AbstractArray{<:AbstractFloat,3},
    ω::AbstractArray{<:AbstractFloat,3},
    γ::Real,
    mi::Real,
    ni::Real,
    me::Real,
    ne::Real,
    kn::Real,
    dt::Real,
    len::Real,
)

    #--- reconstruct initial distribution ---#
    δ = heaviside.(u)
    h0 = @. h0L * δ + h0R * (1.0 - δ)
    h1 = @. h1L * δ + h1R * (1.0 - δ)
    h2 = @. h2L * δ + h2R * (1.0 - δ)

    primL = mixture_conserve_prim(wL, γ)
    primR = mixture_conserve_prim(wR, γ)

    #--- construct interface distribution ---#
    Mu1, Mv1, Mxi1, MuL1, MuR1 = mixture_gauss_moments(primL, inK)
    Muv1 = mixture_moments_conserve(MuL1, Mv1, Mxi1, 0, 0, 0)
    Mu2, Mv2, Mxi2, MuL2, MuR2 = mixture_gauss_moments(primR, inK)
    Muv2 = mixture_moments_conserve(MuR2, Mv2, Mxi2, 0, 0, 0)

    w = @. primL[1] * Muv1 + primR[1] * Muv2
    prim = mixture_conserve_prim(w, γ)
    tau = aap_hs_collision_time(prim, mi, ni, me, ne, kn)
    tau +=
        abs(primL[1, 2] / primL[end, 2] - primR[1, 2] / primR[end, 2]) /
        (primL[1, 2] / primL[end, 2] + primR[1, 2] / primR[end, 2]) *
        dt *
        2.0

    Mt = zeros(2, 2)
    @. Mt[2, :] = tau * (1.0 - exp(-dt / tau)) # f0
    @. Mt[1, :] = dt - Mt[2, :] # M0

    #--- calculate interface flux ---#
    Mu, Mv, Mxi, MuL, MuR = mixture_gauss_moments(prim, inK)

    # flux from M0
    Muv = mixture_moments_conserve(Mu, Mv, Mxi, 1, 0, 0)
    for j = 1:2
        @. fw[:, j] = Mt[1, j] * prim[1, j] * Muv[:, j]
    end

    # flux from f0
    H0 = mixture_maxwellian(u, v, prim)
    H1 = similar(H0)
    H2 = similar(H2)
    for j = 1:2
        H1[:, :, j] = H0[:, :, j] .* prim[4, j]
        H2[:, :, j] .= H0[:, :, j] .* (prim[4, j]^2 + 1.0 / (2.0 * prim[5, j]))
    end

    for j = 1:2
        fw[1, j] += Mt[2, j] * sum(ω[:, :, j] .* u[:, :, j] .* h0[:, :, j])
        fw[2, j] += Mt[2, j] * sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* h0[:, :, j])
        fw[3, j] += Mt[2, j] * sum(ω[:, :, j] .* v[:, :, j] .* u[:, :, j] .* h0[:, :, j])
        fw[4, j] += Mt[2, j] * sum(ω[:, :, j] .* u[:, :, j] .* h1[:, :, j])
        fw[5, j] +=
            Mt[2, j] *
            0.5 *
            (
                sum(
                    ω[:, :, j] .* u[:, :, j] .* (u[:, :, j] .^ 2 .+ v[:, :, j] .^ 2) .*
                    h0[:, :, j],
                ) + sum(ω[:, :, j] .* u[:, :, j] .* h2[:, :, j])
            )

        @. fh0[:, :, j] =
            Mt[1, j] * u[:, :, j] * H0[:, :, j] + Mt[2, j] * u[:, :, j] * h0[:, :, j]
        @. fh1[:, :, j] =
            Mt[1, j] * u[:, :, j] * H1[:, :, j] + Mt[2, j] * u[:, :, j] * h1[:, :, j]
        @. fh2[:, :, j] =
            Mt[1, j] * u[:, :, j] * H2[:, :, j] + Mt[2, j] * u[:, :, j] * h2[:, :, j]
    end

    @. fw *= len
    @. fh0 *= len
    @. fh1 *= len
    @. fh2 *= len

end
