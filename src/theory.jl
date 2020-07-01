# ============================================================
# Theory
# ============================================================


export gauss_moments,
    mixture_gauss_moments,
    moments_conserve,
    mixture_moments_conserve,
    pdf_slope,
    moments_conserve_slope,
    mixture_moments_conserve_slope,
    discrete_moments,
    heat_flux,
    maxwellian,
    mixture_maxwellian,
    reduce_distribution,
    full_distribution,
    conserve_prim,
    mixture_conserve_prim,
    prim_conserve,
    mixture_prim_conserve,
    heat_capacity_ratio,
    ref_vhs_vis,
    sound_speed,
    vhs_collision_time,
    aap_hs_collision_time,
    aap_hs_prim,
    aap_hs_diffeq,
    shift_pdf!,
    em_coefficients,
    hs_boltz_kn,
    kernel_mode,
    boltzmann_fft


"""
Velocity moments of particle distribution function
1. theoretical form

"""

# ------------------------------------------------------------
# Calculate directional velocity moments of Gaussian
# G = (λ / π)^(D / 2) * exp[-λ(c^2 + ξ^2)]
# ------------------------------------------------------------
function gauss_moments(prim::Array{<:Real,1}, inK::Real)

    MuL = OffsetArray{Float64}(undef, 0:6)
    MuR = similar(MuL)
    Mu = similar(MuL)

    MuL[0] = 0.5 * SpecialFunctions.erfc(-sqrt(prim[end]) * prim[2])
    MuL[1] = prim[2] * MuL[0] + 0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])
    MuR[0] = 0.5 * SpecialFunctions.erfc(sqrt(prim[end]) * prim[2])
    MuR[1] = prim[2] * MuR[0] - 0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])

    Threads.@threads for i = 2:6
        MuL[i] = prim[2] * MuL[i-1] + 0.5 * (i - 1) * MuL[i-2] / prim[end]
        MuR[i] = prim[2] * MuR[i-1] + 0.5 * (i - 1) * MuR[i-2] / prim[end]
    end

    @. Mu = MuL + MuR

    if length(prim) == 3

        Mxi = OffsetArray{Float64}(undef, 0:2)
        Mxi[0] = 1.0
        Mxi[1] = 0.5 * inK / prim[end]
        Mxi[2] = (inK^2 + 2.0 * inK) / (4.0 * prim[end]^2)

        return Mu, Mxi, MuL, MuR

    elseif length(prim) == 4

        Mv = OffsetArray{Float64}(undef, 0:6)
        Mv[0] = 1.0
        Mv[1] = prim[3]
        Threads.@threads for i = 2:6
            Mv[i] = prim[3] * Mv[i-1] + 0.5 * (i - 1) * Mv[i-2] / prim[end]
        end

        Mxi = OffsetArray{Float64}(undef, 0:2)
        Mxi[0] = 1.0
        Mxi[1] = 0.5 * inK / prim[end]
        Mxi[2] = (inK^2 + 2.0 * inK) / (4.0 * prim[end]^2)

        return Mu, Mv, Mxi, MuL, MuR

    elseif length(prim) == 5

        Mv = OffsetArray{Float64}(undef, 0:6)
        Mv[0] = 1.0
        Mv[1] = prim[3]
        Threads.@threads for i = 2:6
            Mv[i] = prim[3] * Mv[i-1] + 0.5 * (i - 1) * Mv[i-2] / prim[end]
        end

        Mw = OffsetArray{Float64}(undef, 0:6)
        Mw[0] = 1.0
        Mw[1] = prim[4]
        Threads.@threads for i = 2:6
            Mw[i] = prim[4] * Mw[i-1] + 0.5 * (i - 1) * Mw[i-2] / prim[end]
        end

        return Mu, Mv, Mw, MuL, MuR

    end

end


function mixture_gauss_moments(prim::Array{<:Real,2}, inK::Real)

    Mu = OffsetArray{Float64}(undef, 0:6, axes(prim, 2))
    MuL = similar(Mu)
    MuR = similar(Mu)

    if size(prim, 1) == 3

        Mxi = OffsetArray{Float64}(undef, 0:2, axes(prim, 2))
        for j in axes(prim, 2)
            _tu, _txi, _tuL, _tuR = Kinetic.gauss_moments(prim[:, j], inK)

            Mu[:, j] .= _tu
            Mxi[:, j] .= _txi
            MuL[:, j] .= _tuL
            MuR[:, j] .= _tuR
        end

        return Mu, Mxi, MuL, MuR

    elseif size(prim, 1) == 4

        Mv = OffsetArray{Float64}(undef, 0:6, axes(prim, 2))
        Mxi = OffsetArray{Float64}(undef, 0:2, axes(prim, 2))
        for j in axes(prim, 2)
            _tu, _tv, _txi, _tuL, _tuR = Kinetic.gauss_moments(prim[:, j], inK)

            Mu[:, j] .= _tu
            Mv[:, j] .= _tv
            Mxi[:, j] .= _txi
            MuL[:, j] .= _tuL
            MuR[:, j] .= _tuR
        end

        return Mu, Mv, Mxi, MuL, MuR

    elseif size(prim, 1) == 5

        Mv = OffsetArray{Float64}(undef, 0:6, axes(prim, 2))
        Mw = OffsetArray{Float64}(undef, 0:6, axes(prim, 2))

        for j in axes(prim, 2)
            _tu, _tv, _tw, _tuL, _tuR = Kinetic.gauss_moments(prim[:, j], inK)

            Mu[:, j] .= _tu
            Mv[:, j] .= _tv
            Mw[:, j] .= _tw
            MuL[:, j] .= _tuL
            MuR[:, j] .= _tuR
        end

        return Mu, Mv, Mw, MuL, MuR

    end

end


# ------------------------------------------------------------
# Calculate conservative moments
# ------------------------------------------------------------
function moments_conserve(
    Mu::OffsetArray{<:Real,1},
    Mxi::OffsetArray{<:Real,1},
    alpha::Int,
    delta::Int,
)

    uv = zeros(3)
    uv[1] = Mu[alpha] * Mxi[delta÷2]
    uv[2] = Mu[alpha+1] * Mxi[delta÷2]
    uv[3] = 0.5 * (Mu[alpha+2] * Mxi[delta÷2] + Mu[alpha] * Mxi[(delta+2)÷2])

    return uv

end


function moments_conserve(
    Mu::OffsetArray{<:Real,1},
    Mv::OffsetArray{<:Real,1},
    Mw::OffsetArray{<:Real,1},
    alpha::Int,
    beta::Int,
    delta::Int,
)

    if length(Mw) == 3 # internal motion

        uv = zeros(4)
        uv[1] = Mu[alpha] * Mv[beta] * Mw[delta÷2]
        uv[2] = Mu[alpha+1] * Mv[beta] * Mw[delta÷2]
        uv[3] = Mu[alpha] * Mv[beta+1] * Mw[delta÷2]
        uv[4] =
            0.5 * (
                Mu[alpha+2] * Mv[beta] * Mw[delta÷2] +
                Mu[alpha] * Mv[beta+2] * Mw[delta÷2] +
                Mu[alpha] * Mv[beta] * Mw[(delta+2)÷2]
            )

    else

        uv = zeros(5)
        uv[1] = Mu[alpha] * Mv[beta] * Mw[delta]
        uv[2] = Mu[alpha+1] * Mv[beta] * Mw[delta]
        uv[3] = Mu[alpha] * Mv[beta+1] * Mw[delta]
        uv[4] = Mu[alpha] * Mv[beta] * Mw[delta+1]
        uv[5] =
            0.5 * (
                Mu[alpha+2] * Mv[beta] * Mw[delta] +
                Mu[alpha] * Mv[beta+2] * Mw[delta] +
                Mu[alpha] * Mv[beta] * Mw[delta+2]
            )

    end

    return uv

end


function mixture_moments_conserve(
    Mu::OffsetArray{<:Real,2},
    Mxi::OffsetArray{<:Real,2},
    alpha::Int,
    delta::Int,
)

    Muv = zeros(3, axes(Mu, 2))
    for j in axes(Muv, 2)
        Muv[:, j] .= moments_conserve(Mu[:, j], Mxi[:, j], alpha, delta)
    end

    return Muv

end


function mixture_moments_conserve(
    Mu::OffsetArray{<:Real,2},
    Mv::OffsetArray{<:Real,2},
    Mw::OffsetArray{<:Real,2},
    alpha::Int,
    beta::Int,
    delta::Int,
)

    Muv = ifelse(length(Mw) == 3, zeros(4, axes(Mu, 2)), zeros(5, axes(Mu, 2)))
    for j in axes(Muv, 2)
        Muv[:, j] .= moments_conserve(Mu[:, j], Mv[:, j], Mw[:, j], alpha, beta, delta)
    end

    return Muv

end


# ------------------------------------------------------------
# Calculate slope of particle distribution function
# a = a1 + u * a2 + 0.5 * u^2 * a3
# ------------------------------------------------------------
function pdf_slope(prim::Array{<:Real,1}, sw::Array{<:Real,1}, inK::Real)

    sl = zeros(axes(prim))

    if length(prim) == 3

        sl[3] =
            4.0 * prim[3]^2 / (inK + 1.0) / prim[1] * (
                2.0 * sw[3] - 2.0 * prim[2] * sw[2] +
                sw[1] * (prim[2]^2 - 0.5 * (inK + 1.0) / prim[3])
            )
        sl[2] = 2.0 * prim[3] / prim[1] * (sw[2] - prim[2] * sw[1]) - prim[2] * sl[3]
        sl[1] =
            sw[1] / prim[1] - prim[2] * sl[2] -
            0.5 * (prim[2]^2 + 0.5 * (inK + 1.0) / prim[3]) * sl[3]

    elseif length(prim) == 4

        sl[4] =
            4.0 * prim[4]^2 / (inK + 2.0) / prim[1] * (
                2.0 * sw[4] - 2.0 * prim[2] * sw[2] - 2.0 * prim[3] * sw[3] +
                sw[1] * (prim[2]^2 + prim[3]^2 - 0.5 * (inK + 2.0) / prim[4])
            )
        sl[3] = 2.0 * prim[4] / prim[1] * (sw[3] - prim[3] * sw[1]) - prim[3] * sl[4]
        sl[2] = 2.0 * prim[4] / prim[1] * (sw[2] - prim[2] * sw[1]) - prim[2] * sl[4]
        sl[1] =
            sw[1] / prim[1] - prim[2] * sl[2] - prim[3] * sl[3] -
            0.5 * (prim[2]^2 + prim[3]^2 + 0.5 * (inK + 2.0) / prim[4]) * sl[4]

    elseif length(prim) == 5

        sl[5] =
            4.0 * prim[5]^2 / (inK + 3.0) / prim[1] * (
                2.0 * sw[5] - 2.0 * prim[2] * sw[2] - 2.0 * prim[3] * sw[3] -
                2.0 * prim[4] * sw[4] +
                sw[1] * (prim[2]^2 + prim[3]^2 + prim[4]^2 - 0.5 * (inK + 3.0) / prim[5])
            )
        sl[4] = 2.0 * prim[5] / prim[1] * (sw[4] - prim[4] * sw[1]) - prim[4] * sl[5]
        sl[3] = 2.0 * prim[5] / prim[1] * (sw[3] - prim[3] * sw[1]) - prim[3] * sl[5]
        sl[2] = 2.0 * prim[5] / prim[1] * (sw[2] - prim[2] * sw[1]) - prim[2] * sl[5]
        sl[1] =
            sw[1] / prim[1] - prim[2] * sl[2] - prim[3] * sl[3] - prim[4] * sl[4] -
            0.5 * (prim[2]^2 + prim[3]^2 + prim[4]^2 + 0.5 * (inK + 3.0) / prim[5]) * sl[5]

    end

    return sl

end


# ------------------------------------------------------------
# Calculate slope-related conservative moments
# a = a1 + u * a2 + 0.5 * u^2 * a3
# ------------------------------------------------------------
moments_conserve_slope(
    a::Array{<:Real,1},
    Mu::OffsetArray{<:Real,1},
    Mxi::OffsetArray{<:Real,1},
    alpha::Int,
) =
    a[1] .* moments_conserve(Mu, Mxi, alpha + 0, 0) .+
    a[2] .* moments_conserve(Mu, Mxi, alpha + 1, 0) .+
    0.5 * a[3] .* moments_conserve(Mu, Mxi, alpha + 2, 0) .+
    0.5 * a[3] .* moments_conserve(Mu, Mxi, alpha + 0, 2)


function moments_conserve_slope(
    a::Array{<:Real,1},
    Mu::OffsetArray{<:Real,1},
    Mv::OffsetArray{<:Real,1},
    Mxi::OffsetArray{<:Real,1},
    alpha::Int,
    beta::Int,
)

    au =
        a[1] .* moments_conserve(Mu, Mv, Mxi, alpha + 0, beta + 0, 0) .+
        a[2] .* moments_conserve(Mu, Mv, Mxi, alpha + 1, beta + 0, 0) .+
        a[3] .* moments_conserve(Mu, Mv, Mxi, alpha + 0, beta + 1, 0) .+
        0.5 * a[4] .* moments_conserve(Mu, Mv, Mxi, alpha + 2, beta + 0, 0) .+
        0.5 * a[4] .* moments_conserve(Mu, Mv, Mxi, alpha + 0, beta + 2, 0) .+
        0.5 * a[4] .* moments_conserve(Mu, Mv, Mxi, alpha + 0, beta + 0, 2)

    return au

end


function moments_conserve_slope(
    a::Array{<:Real,1},
    Mu::OffsetArray{<:Real,1},
    Mv::OffsetArray{<:Real,1},
    Mw::OffsetArray{<:Real,1},
    alpha::Int,
    beta::Int,
    delta::Int,
)

    au =
        a[1] .* moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 0, delta + 0) .+
        a[2] .* moments_conserve(Mu, Mv, Mw, alpha + 1, beta + 0, delta + 0) .+
        a[3] .* moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 1, delta + 0) .+
        a[4] .* moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 0, delta + 1) .+
        0.5 * a[5] .* moments_conserve(Mu, Mv, Mw, alpha + 2, beta + 0, delta + 0) .+
        0.5 * a[5] .* moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 2, delta + 0) .+
        0.5 * a[5] .* moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 0, delta + 2)

    return au

end


function mixture_moments_conserve_slope(
    a::Array{<:Real,2},
    Mu::OffsetArray{<:Real,2},
    Mxi::OffsetArray{<:Real,2},
    alpha::Int,
)

    au = zeros(3, axes(a, 2))
    for j in axes(au, 2)
        au[:, j] .= moments_conserve_slope(a[:, j], Mu[:, j], Mxi[:, j], alpha)
    end

    return au

end


function mixture_moments_conserve_slope(
    a::Array{<:Real,2},
    Mu::OffsetArray{<:Real,2},
    Mv::OffsetArray{<:Real,2},
    Mxi::OffsetArray{<:Real,2},
    alpha::Int,
    beta::Int,
)

    au = zeros(4, axes(a, 2))
    for j in axes(au, 2)
        au[:, j] .=
            moments_conserve_slope(a[:, j], Mu[:, j], Mv[:, j], Mxi[:, j], alpha, beta)
    end

    return au

end


function mixture_moments_conserve_slope(
    a::Array{<:Real,2},
    Mu::OffsetArray{<:Real,2},
    Mv::OffsetArray{<:Real,2},
    Mw::OffsetArray{<:Real,2},
    alpha::Int,
    beta::Int,
    delta::Int,
)

    au = zeros(typeof(a[1]), 5, axes(a, 2))
    for j in axes(au, 2)
        au[:, j] .= moments_conserve_slope(
            a[:, j],
            Mu[:, j],
            Mv[:, j],
            Mw[:, j],
            alpha,
            beta,
            delta,
        )
    end

    return au

end


"""
Velocity moments of particle distribution function
2. discrete form
"""

# ------------------------------------------------------------
# Velocity moments of order n
# ------------------------------------------------------------
# --- 1D ---#
discrete_moments(
    f::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    ω::AbstractArray{<:AbstractFloat,1},
    n::Int,
) = sum(@. ω * u^n * f)


# --- 2D ---#
discrete_moments(
    f::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
    n::Int,
) = sum(@. ω * u^n * f)


# --- 3D ---#
discrete_moments(
    f::AbstractArray{<:AbstractFloat,3},
    u::AbstractArray{<:AbstractFloat,3},
    ω::AbstractArray{<:AbstractFloat,3},
    n::Int,
) = sum(@. ω * u^n * f)


# ------------------------------------------------------------
# Conservative moments
# ------------------------------------------------------------
# --- 1D ---#
moments_conserve(
    f::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    ω::AbstractArray{<:AbstractFloat,1},
) = [
    discrete_moments(f, u, ω, 0),
    discrete_moments(f, u, ω, 1),
    0.5 * discrete_moments(f, u, ω, 2),
]


moments_conserve(
    h::AbstractArray{<:AbstractFloat,1},
    b::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    ω::AbstractArray{<:AbstractFloat,1},
) = [
    discrete_moments(h, u, ω, 0),
    discrete_moments(h, u, ω, 1),
    0.5 * (discrete_moments(h, u, ω, 2) + discrete_moments(b, u, ω, 0)),
]


function mixture_moments_conserve(
    f::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
)
    w = zeros(eltype(f), 3, size(f, 2))
    for j in axes(w, 2)
        w[:, j] .= moments_conserve(f[:, j], u, ω)
    end

    return w
end


function mixture_moments_conserve(
    h::AbstractArray{<:AbstractFloat,2},
    b::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
)
    w = zeros(eltype(h), 3, size(h, 2))
    for j in axes(w, 2)
        w[:, j] .= moments_conserve(h[:, j], b[:, j], u, ω)
    end

    return w
end


# --- 2D ---#
moments_conserve(
    f::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
) = [
    discrete_moments(f, u, ω, 0),
    discrete_moments(f, u, ω, 1),
    discrete_moments(f, v, ω, 1),
    0.5 * (discrete_moments(f, u, ω, 2) + discrete_moments(f, v, ω, 2)),
]


moments_conserve(
    h::AbstractArray{<:AbstractFloat,2},
    b::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
) = [
    discrete_moments(h, u, ω, 0),
    discrete_moments(h, u, ω, 1),
    discrete_moments(h, v, ω, 1),
    0.5 * (
        discrete_moments(h, u, ω, 2) +
        discrete_moments(h, v, ω, 2) +
        discrete_moments(b, u, ω, 0)
    ),
]


function mixture_moments_conserve(
    f::AbstractArray{<:AbstractFloat,3},
    u::AbstractArray{<:AbstractFloat,3},
    v::AbstractArray{<:AbstractFloat,3},
    ω::AbstractArray{<:AbstractFloat,3},
)
    w = zeros(eltype(f), 4, size(f, 3))
    for j in axes(w, 2)
        w[:, j] .= moments_conserve(f[:, :, j], u[:, :, j], v[:, :, j], ω[:, :, j])
    end

    return w
end


function mixture_moments_conserve(
    h::AbstractArray{<:AbstractFloat,2},
    b::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
)
    w = zeros(eltype(h), 4, size(h, 3))
    for j in axes(w, 2)
        w[:, j] .=
            moments_conserve(h[:, :, j], b[:, :, j], u[:, :, j], v[:, :, j], ω[:, :, j])
    end

    return w
end


# --- 3D ---#
moments_conserve(
    f::AbstractArray{<:AbstractFloat,3},
    u::AbstractArray{<:AbstractFloat,3},
    v::AbstractArray{<:AbstractFloat,3},
    w::AbstractArray{<:AbstractFloat,3},
    ω::AbstractArray{<:AbstractFloat,3},
) = [
    discrete_moments(f, u, ω, 0),
    discrete_moments(f, u, ω, 1),
    discrete_moments(f, v, ω, 1),
    discrete_moments(f, w, ω, 1),
    0.5 * (
        discrete_moments(f, u, ω, 2) +
        discrete_moments(f, v, ω, 2) +
        discrete_moments(f, w, ω, 2)
    ),
]


moments_conserve(
    h0::AbstractArray{<:AbstractFloat,1},
    h1::AbstractArray{<:AbstractFloat,1},
    h2::AbstractArray{<:AbstractFloat,1},
    h3::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    ω::AbstractArray{<:AbstractFloat,1},
) = [
    discrete_moments(h0, u, ω, 0),
    discrete_moments(h0, u, ω, 1),
    discrete_moments(h1, u, ω, 0),
    discrete_moments(h2, u, ω, 0),
    0.5 * discrete_moments(h0, u, ω, 2) + 0.5 * discrete_moments(h3, u, ω, 0),
]


function mixture_moments_conserve(
    f::AbstractArray{<:Real,4},
    u::AbstractArray{<:Real,4},
    v::AbstractArray{<:Real,4},
    w::AbstractArray{<:Real,4},
    ω::AbstractArray{<:Real,4},
)
    w = zeros(eltype(f), 5, size(f, 4))
    for j in axes(w, 2)
        w[:, j] .= moments_conserve(
            f[:, :, :, j],
            u[:, :, :, j],
            v[:, :, :, j],
            w[:, :, :, j],
            ω[:, :, :, j],
        )
    end

    return w
end


function mixture_moments_conserve(
    h0::AbstractArray{<:AbstractFloat,2},
    h1::AbstractArray{<:AbstractFloat,2},
    h2::AbstractArray{<:AbstractFloat,2},
    h3::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:Real,2},
    ω::AbstractArray{<:Real,2},
)
    w = zeros(eltype(h0), 5, size(h0, 2))
    for j in axes(w, 2)
        w[:, j] .=
            moments_conserve(h0[:, j], h1[:, j], h2[:, j], h3[:, j], u[:, j], ω[:, j])
    end

    return w
end


"""
Stress tensor from particle distribution function

"""

function stress_tensor(
    f::AbstractArray{<:AbstractFloat,2},
    prim::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
)
    P = zeros(eltype(prim), 2, 2)

    P[1, 1] = sum(@. ω * (u - prim[2]) * (u - prim[2]) * f)
    P[1, 2] = sum(@. ω * (u - prim[2]) * (v - prim[3]) * f)
    P[2, 1] = P[1, 2]
    P[1, 2] = sum(@. ω * (v - prim[3]) * (v - prim[3]) * f)

    return P
end


"""
Heat flux from particle distribution function

"""

# --- 1D ---#
heat_flux(
    h::AbstractArray{<:AbstractFloat,1},
    prim::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    ω::AbstractArray{<:AbstractFloat,1},
) = 0.5 * sum(@. ω * (u - prim[2]) * (u - prim[2])^2 * h)


heat_flux(
    h::AbstractArray{<:AbstractFloat,1},
    b::AbstractArray{<:AbstractFloat,1},
    prim::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    ω::AbstractArray{<:AbstractFloat,1},
) = 0.5 * (sum(@. ω * (u - prim[2]) * (u - prim[2])^2 * h) + sum(@. ω * (u - prim[2]) * b))


# --- 2D ---#
function heat_flux(
    h::AbstractArray{<:AbstractFloat,2},
    prim::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
)
    q = zeros(eltype(prim), 2)

    q[1] = 0.5 * sum(@. ω * (u - prim[2]) * ((u - prim[2])^2 + (v - prim[3])^2) * h)
    q[2] = 0.5 * sum(@. ω * (v - prim[3]) * ((u - prim[2])^2 + (v - prim[3])^2) * h)

    return q
end


function heat_flux(
    h::AbstractArray{<:AbstractFloat,2},
    b::AbstractArray{<:AbstractFloat,2},
    prim::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
)
    q = zeros(eltype(prim), 2)

    q[1] =
        0.5 * (
            sum(@. ω * (u - prim[2]) * ((u - prim[2])^2 + (v - prim[3])^2) * h) +
            sum(@. ω * (u - prim[2]) * b)
        )
    q[2] =
        0.5 * (
            sum(@. ω * (v - prim[3]) * ((u - prim[2])^2 + (v - prim[3])^2) * h) +
            sum(@. ω * (v - prim[3]) * b)
        )

    return q
end


"""
Equilibrium in discrete form
1. Gas: Maxwellian

>@param[in] : particle velocity quadrature points
>@param[in] : density, velocity and inverse of temperature
>@return : Maxwellian distribution function

"""

# --- 1D ---#
maxwellian(u::AbstractArray{<:Real,1}, ρ::Real, U::Real, λ::Real) =
    @. ρ * (λ / π)^0.5 * exp(-λ * (u - U)^2)


maxwellian(u::AbstractArray{<:Real,1}, prim::Array{<:Real,1}) =
    maxwellian(u, prim[1], prim[2], prim[end]) # in case of input with length 4/5


function mixture_maxwellian(u::AbstractArray{<:Real,2}, prim::Array{<:Real,2})
    mixM = similar(u)
    for j in axes(mixM, 2)
        mixM[:, j] .= maxwellian(u[:, j], prim[:, j])
    end

    return mixM
end


# --- 2D ---#
maxwellian(
    u::AbstractArray{<:Real,2},
    v::AbstractArray{<:Real,2},
    ρ::Real,
    U::Real,
    V::Real,
    λ::Real,
) = @. ρ * (λ / π) * exp(-λ * ((u - U)^2 + (v - V)^2))


maxwellian(u::AbstractArray{<:Real,2}, v::AbstractArray{<:Real,2}, prim::Array{<:Real,1}) =
    maxwellian(u, v, prim[1], prim[2], prim[3], prim[end]) # in case of input with length 5


function mixture_maxwellian(
    u::AbstractArray{<:Real,3},
    v::AbstractArray{<:Real,3},
    prim::Array{<:Real,2},
)
    mixM = zeros(axes(u))
    for k in axes(mixM, 3)
        mixM[:, :, k] .= maxwellian(u[:, :, k], v[:, :, k], prim[:, k])
    end

    return mixM
end


# --- 3D ---#
maxwellian(
    u::AbstractArray{<:Real,3},
    v::AbstractArray{<:Real,3},
    w::AbstractArray{<:Real,3},
    ρ::Real,
    U::Real,
    V::Real,
    W::Real,
    λ::Real,
) = @. ρ * (λ / π)^1.5 * exp(-λ * ((u - U)^2 + (v - V)^2 + (w - W)^2))


maxwellian(
    u::AbstractArray{<:Real,3},
    v::AbstractArray{<:Real,3},
    w::AbstractArray{<:Real,3},
    prim::Array{<:Real,1},
) = maxwellian(u, v, w, prim[1], prim[2], prim[3], prim[4], prim[5])


function mixture_maxwellian(
    u::AbstractArray{<:Real,4},
    v::AbstractArray{<:Real,4},
    w::AbstractArray{<:Real,4},
    prim::Array{<:Real,2},
)
    mixM = zeros(axes(u))
    for l in axes(mixM, 4)
        mixM[:, :, :, l] .=
            maxwellian(u[:, :, :, l], v[:, :, :, l], w[:, :, :, l], prim[:, l])
    end

    return mixM
end


"""
Reduced distribution function

>@param[in] : particle distribution function with 3D velocity space
>@param[in] : quadrature weights with 2D setting (v & w by default)

"""

function reduce_distribution(
    f::AbstractArray{<:AbstractFloat,3},
    weights::AbstractArray{<:AbstractFloat,2},
    dim = 1::Int,
)

    if dim == 1
        h = similar(f, axes(f, 1))
        for i in eachindex(h)
            h[i] = sum(@. weights * f[i, :, :])
        end
    elseif dim == 2
        h = similar(f, axes(f, 2))
        for j in eachindex(h)
            h[j] = sum(@. weights * f[:, j, :])
        end
    elseif dim == 3
        h = similar(f, axes(f, 3))
        for k in eachindex(h)
            h[k] = sum(@. weights * f[:, :, k])
        end
    else
    end

    return h
end


function reduce_distribution(
    f::AbstractArray{<:AbstractFloat,3},
    v::AbstractArray{<:AbstractFloat,3},
    w::AbstractArray{<:AbstractFloat,3},
    weights::AbstractArray{<:AbstractFloat,2},
    dim = 1::Int,
)

    if dim == 1
        h = similar(f, axes(f, 1))
        b = similar(h)
        for i in eachindex(h)
            h[i] = sum(@. weights * f[i, :, :])
            b[i] = sum(@. weights * (v[i, :, :]^2 + w[i, :, :]^2) * f[i, :, :])
        end
    elseif dim == 2
        h = similar(f, axes(f, 2))
        b = similar(h)
        for j in eachindex(h)
            h[j] = sum(@. weights * f[:, j, :])
            b[j] = sum(@. weights * (v[:, j, :]^2 + w[:, j, :]^2) * f[:, j, :])
        end
    elseif dim == 3
        h = similar(f, axes(f, 3))
        b = similar(h)
        for k in eachindex(h)
            h[k] = sum(@. weights * f[:, :, k])
            b[k] = sum(@. weights * (v[:, :, k]^2 + w[:, :, k]^2) * f[:, :, k])
        end
    else
    end

    return h, b
end


"""
Recover full distribution function from reduced ones

>@param[in] h & b : reduced particle distribution function with 1D velocity space
>@param[in] u : quadrature nodes in 1D velocity space
>@param[in] weights : quadrature weights in 1D velocity space
>@param[in] v & w : quadrature nodes in the rest velocity space (with 3D setting)

>@param[out] f : particle distribution function with 3D velocity space

"""

function full_distribution(
    h::AbstractArray{<:AbstractFloat,1},
    b::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    weights::AbstractArray{<:AbstractFloat,1},
    v::AbstractArray{<:AbstractFloat,3},
    w::AbstractArray{<:AbstractFloat,3},
    ρ::Real,
    γ = 5 / 3::Real,
)

    @assert length(h) == size(v, 1) throw(DimensionMismatch("reduced and full distribution function mismatch"))

    Ei = 0.5 * discrete_moments(b, u, weights, 0)
    λi = 0.5 * ρ / (γ - 1.0) / Ei / 3.0 * 2.0

    f = similar(v)
    for k in axes(f, 3), j in axes(f, 2), i in axes(f, 1)
        f[i, j, k] = h[i] * (λi / π) * exp(-λi * v[i, j, k]^2) * exp(-λi * w[i, j, k]^2)
    end

    return f
end


full_distribution(
    h::AbstractArray{<:AbstractFloat,1},
    b::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    weights::AbstractArray{<:AbstractFloat,1},
    v::AbstractArray{<:AbstractFloat,3},
    w::AbstractArray{<:AbstractFloat,3},
    prim::Array{<:Real,1},
    γ = 5 / 3::Real,
) = full_distribution(h, b, u, weights, v, w, prim[1], γ)


"""
Transforms between conservative and primitive variables

"""

# ------------------------------------------------------------
# primitive -> conservative
# ------------------------------------------------------------
function prim_conserve(prim::Array{<:Real,1}, γ::Real)

    W = zeros(axes(prim))

    if length(prim) == 3 # 1D
        W[1] = prim[1]
        W[2] = prim[1] * prim[2]
        W[3] = 0.5 * prim[1] / prim[3] / (γ - 1.0) + 0.5 * prim[1] * prim[2]^2
    elseif length(prim) == 4 # 2D
        W[1] = prim[1]
        W[2] = prim[1] * prim[2]
        W[3] = prim[1] * prim[3]
        W[4] = 0.5 * prim[1] / prim[4] / (γ - 1.0) + 0.5 * prim[1] * (prim[2]^2 + prim[3]^2)
    elseif length(prim) == 5 # 3D
        W[1] = prim[1]
        W[2] = prim[1] * prim[2]
        W[3] = prim[1] * prim[3]
        W[4] = prim[1] * prim[4]
        W[5] =
            0.5 * prim[1] / prim[5] / (γ - 1.0) +
            0.5 * prim[1] * (prim[2]^2 + prim[3]^2 + prim[4]^2)
    else
        println("prim -> w : dimension error")
    end

    return W

end


prim_conserve(ρ::Real, U::Real, λ::Real, γ::Real) = prim_conserve([ρ, U, λ], γ)


prim_conserve(ρ::Real, U::Real, V::Real, λ::Real, γ::Real) = prim_conserve([ρ, U, V, λ], γ)


function mixture_prim_conserve(prim::Array{<:Real,2}, γ::Real)
    w = zeros(axes(prim))
    for j in axes(w, 2)
        w[:, j] .= prim_conserve(prim[:, j], γ)
    end

    return w
end


# ------------------------------------------------------------
# conservative -> primitive
# ------------------------------------------------------------
function conserve_prim(W::Array{<:Real,1}, γ::Real)

    prim = zeros(axes(W))

    if length(W) == 3 # 1D
        prim[1] = W[1]
        prim[2] = W[2] / W[1]
        prim[3] = 0.5 * W[1] / (γ - 1.0) / (W[3] - 0.5 * W[2]^2 / W[1])
    elseif length(W) == 4 # 2D
        prim[1] = W[1]
        prim[2] = W[2] / W[1]
        prim[3] = W[3] / W[1]
        prim[4] = 0.5 * W[1] / (γ - 1.0) / (W[4] - 0.5 * (W[2]^2 + W[3]^2) / W[1])
    elseif length(W) == 5 # 3D
        prim[1] = W[1]
        prim[2] = W[2] / W[1]
        prim[3] = W[3] / W[1]
        prim[4] = W[4] / W[1]
        prim[5] = 0.5 * W[1] / (γ - 1.0) / (W[5] - 0.5 * (W[2]^2 + W[3]^2 + W[4]^2) / W[1])
    else
        println("w -> prim : dimension error")
    end

    return prim

end


conserve_prim(ρ::Real, M::Real, E::Real, gamma::Real) = conserve_prim([ρ, M, E], gamma)


conserve_prim(ρ::Real, MX::Real, MY::Real, E::Real, gamma::Real) =
    conserve_prim([ρ, MX, MY, E], gamma)


function mixture_conserve_prim(w::Array{<:Real,2}, γ::Real)
    prim = zeros(axes(w))
    for j in axes(prim, 2)
        prim[:, j] .= conserve_prim(w[:, j], γ)
    end

    return prim
end


"""
Thermodynamical properties
"""

# ------------------------------------------------------------
# Calculate heat capacity ratio
# ------------------------------------------------------------
function heat_capacity_ratio(K::Real, D::Int)

    if D == 1
        γ = (K + 3.0) / (K + 1.0)
    elseif D == 2
        γ = (K + 4.0) / (K + 2.0)
    elseif D == 3
        γ = (K + 5.0) / (K + 3.0)
    end

    return γ

end


# ------------------------------------------------------------
# Calculate speed of sound
# ------------------------------------------------------------
sound_speed(λ::Real, γ::Real) = (0.5 * γ / λ)^0.5


sound_speed(prim::Array{<:Real,1}, γ::Real) = sound_speed(prim[end], γ)


function sound_speed(prim::Array{<:Real,2}, γ::Real)
    c = zeros(axes(prim, 2))
    for j in eachindex(c)
        c[j] = sound_speed(prim[end, j], γ)
    end

    return maximum(c)
end


"""
Single component gas models

"""

# ------------------------------------------------------------
# Calculate reference viscosity
# 1. variable hard sphere (VHS) model
# ------------------------------------------------------------
ref_vhs_vis(Kn::Real, alpha::Real, omega::Real) =
    5.0 * (alpha + 1.0) * (alpha + 2.0) * √π /
    (4.0 * alpha * (5.0 - 2.0 * omega) * (7.0 - 2.0 * omega)) * Kn


# ------------------------------------------------------------
# Calculate collision time
# 1. variable hard sphere (VHS) model
# ------------------------------------------------------------
vhs_collision_time(prim::Array{<:Real,1}, muRef::Real, omega::Real) =
    muRef * 2.0 * prim[end]^(1.0 - omega) / prim[1]


# ------------------------------------------------------------
# Calculate effective Knudsen number for fast spectral method
# 1. hard sphere (HS) model
# ------------------------------------------------------------
hs_boltz_kn(mu_ref, alpha) =
    64 * sqrt(2.0)^alpha / 5.0 * gamma((alpha + 3) / 2) * gamma(2.0) * sqrt(pi) * mu_ref


```
Fast spectral method for Boltzmann collision operator

```

# ------------------------------------------------------------
# Calculate collision kernel
# ------------------------------------------------------------
function kernel_mode(
    M::Int,
    umax::Real,
    vmax::Real,
    wmax::Real,
    du::Real,
    dv::Real,
    dw::Real,
    unum::Int,
    vnum::Int,
    wnum::Int,
    alpha::Real;
    quad_num = 64,
)

    supp = sqrt(2.0) * 2.0 * max(umax, vmax, wmax) / (3.0 + sqrt(2.0))

    fre_vx = range(-π / du, (unum ÷ 2 - 1) * 2.0 * π / unum / du, length = unum)
    fre_vy = range(-π / dv, (vnum ÷ 2 - 1) * 2.0 * π / vnum / dv, length = vnum)
    fre_vz = range(-π / dw, (wnum ÷ 2 - 1) * 2.0 * π / wnum / dw, length = wnum)

    # abscissa, gweight = gausslegendre(quad_num)
    # @. abscissa = (0. * (1. - abscissa) + supp * (1. + abscissa)) / 2
    # @. gweight *= (supp - 0.) / 2

    abscissa, gweight = lgwt(quad_num, 0, supp)

    phi = zeros(unum, vnum, wnum, M * (M - 1))
    psi = zeros(unum, vnum, wnum, M * (M - 1))
    phipsi = zeros(unum, vnum, wnum)
    for loop = 1:M-1
        theta = π / M * loop
        for loop2 = 1:M
            theta2 = π / M * loop2
            idx = (loop - 1) * M + loop2
            for k = 1:wnum, j = 1:vnum, i = 1:unum
                s =
                    fre_vx[i] * sin(theta) * cos(theta2) +
                    fre_vy[j] * sin(theta) * sin(theta2) +
                    fre_vz[k] * cos(theta)
                # phi
                int_temp = 0.0
                for id = 1:quad_num
                    int_temp +=
                        2.0 * gweight[id] * cos(s * abscissa[id]) * (abscissa[id]^alpha)
                end
                phi[i, j, k, idx] = int_temp * sin(theta)
                # psi
                s = fre_vx[i]^2 + fre_vy[j]^2 + fre_vz[k]^2 - s^2
                if s <= 0.0
                    psi[i, j, k, idx] = π * supp^2
                else
                    s = sqrt(s)
                    bel = supp * s
                    bessel = besselj(1, bel)
                    psi[i, j, k, idx] = 2.0 * π * supp * bessel / s
                end
                # phipsi
                phipsi[i, j, k] += phi[i, j, k, idx] * psi[i, j, k, idx]
            end
        end
    end

    return phi, psi, phipsi

end


# ------------------------------------------------------------
# calculate collision operator with FFT
# ------------------------------------------------------------
function boltzmann_fft(
    f::AbstractArray{<:Real,3},
    Kn::Real,
    M::Int,
    ϕ::AbstractArray{<:Real,4},
    ψ::AbstractArray{<:Real,4},
    phipsi::AbstractArray{<:Real,3},
)

    f_spec = f .+ 0im
    bfft!(f_spec)
    f_spec ./= size(f, 1) * size(f, 2) * size(f, 3)
    f_spec .= fftshift(f_spec)

    # --- gain term ---#
    f_temp = zeros(axes(f_spec)) .+ 0im
    for i = 1:M*(M-1)
        fg1 = f_spec .* ϕ[:, :, :, i]
        fg2 = f_spec .* ψ[:, :, :, i]
        fg11 = fft(fg1)
        fg22 = fft(fg2)
        f_temp .+= fg11 .* fg22
    end

    # --- loss term ---#
    fl1 = f_spec .* phipsi
    fl2 = f_spec
    fl11 = fft(fl1)
    fl22 = fft(fl2)
    f_temp .-= fl11 .* fl22

    Q = @. 4.0 * π^2 / Kn / M^2 * real(f_temp)

    return Q

end


"""
Multiple component gas models
"""

# ------------------------------------------------------------
# Calculate mixture collision time from AAP model
# ------------------------------------------------------------
function aap_hs_collision_time(
    prim::Array{<:Real,2},
    mi::Real,
    ni::Real,
    me::Real,
    ne::Real,
    kn::Real,
)

    ν = zeros(axes(prim, 2))

    ν[1] =
        prim[1, 1] / (mi * (ni + ne)) * 4.0 * sqrt(π) / 3.0 *
        sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 1]) / (sqrt(2.0) * π * kn) +
        prim[1, 2] / (me * (ni + ne)) * 4.0 * sqrt(π) / 3.0 *
        sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2]) / (sqrt(2.0) * π * kn)
    ν[2] =
        prim[1, 1] / (mi * (ni + ne)) * 4.0 * sqrt(π) / 3.0 *
        sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2]) / (sqrt(2.0) * π * kn) +
        prim[1, 2] / (me * (ni + ne)) * 4.0 * sqrt(π) / 3.0 *
        sqrt(1.0 / prim[end, 2] + 1.0 / prim[end, 2]) / (sqrt(2.0) * π * kn)

    return 1.0 ./ ν

end


# ------------------------------------------------------------
# Calculate mixture primitive variables from AAP model
# ------------------------------------------------------------
function aap_hs_prim(
    prim::Array{<:Real,2},
    tau::Array{<:Real,1},
    mi::Real,
    ni::Real,
    me::Real,
    ne::Real,
    kn::Real,
)

    mixprim = similar(prim)

    if size(prim, 1) == 3
        mixprim[1, :] = deepcopy(prim[1, :])
        mixprim[2, 1] =
            prim[2, 1] +
            tau[1] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[2, 2] - prim[2, 1])
        mixprim[2, 2] =
            prim[2, 2] +
            tau[2] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[2, 1] - prim[2, 2])
        mixprim[3, 1] =
            1.0 / (
                1.0 / prim[end, 1] - 2.0 / 3.0 * (mixprim[2, 1] - prim[2, 1])^2 +
                tau[1] / kn * 2.0 * mi / (mi + me) *
                (
                    4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                    sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
                ) *
                (
                    1.0 / prim[end, 2] * me / mi - 1.0 / prim[end, 1] +
                    2.0 / 3.0 * me / mi * (prim[2, 2] - prim[2, 1])^2
                )
            )
        mixprim[3, 2] =
            1.0 / (
                1.0 / prim[end, 2] - 2.0 / 3.0 * (mixprim[2, 2] - prim[2, 2])^2 +
                tau[2] / kn * 2.0 * me / (mi + me) *
                (
                    4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                    sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
                ) *
                (
                    1.0 / prim[end, 1] * mi / me - 1.0 / prim[end, 2] +
                    2.0 / 3.0 * mi / me * (prim[2, 1] - prim[2, 2])^2
                )
            )
    elseif size(prim, 1) == 4
        mixprim[1, :] = deepcopy(prim[1, :])
        mixprim[2, 1] =
            prim[2, 1] +
            tau[1] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[2, 2] - prim[2, 1])
        mixprim[2, 2] =
            prim[2, 2] +
            tau[2] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[2, 1] - prim[2, 2])
        mixprim[3, 1] =
            prim[3, 1] +
            tau[1] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[3, 2] - prim[3, 1])
        mixprim[3, 2] =
            prim[3, 2] +
            tau[2] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[3, 1] - prim[3, 2])
        mixprim[4, 1] =
            1.0 / (
                1.0 / prim[end, 1] - 2.0 / 3.0 * (mixprim[2, 1] - prim[2, 1])^2 -
                2.0 / 3.0 * (mixprim[3, 1] - prim[3, 1])^2 +
                tau[1] / kn * 2.0 * mi / (mi + me) *
                (
                    4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                    sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
                ) *
                (
                    1.0 / prim[end, 2] * me / mi - 1.0 / prim[end, 1] +
                    2.0 / 3.0 * me / mi * (prim[2, 2] - prim[2, 1])^2 +
                    2.0 / 3.0 * me / mi * (prim[3, 2] - prim[3, 1])^2
                )
            )
        mixprim[4, 2] =
            1.0 / (
                1.0 / prim[end, 2] - 2.0 / 3.0 * (mixprim[2, 2] - prim[2, 2])^2 -
                2.0 / 3.0 * (mixprim[3, 2] - prim[3, 2])^2 +
                tau[2] / kn * 2.0 * me / (mi + me) *
                (
                    4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                    sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
                ) *
                (
                    1.0 / prim[end, 1] * mi / me - 1.0 / prim[end, 2] +
                    2.0 / 3.0 * mi / me * (prim[2, 1] - prim[2, 2])^2 +
                    2.0 / 3.0 * mi / me * (prim[3, 1] - prim[3, 2])^2
                )
            )
    elseif size(prim, 1) == 5
        mixprim[1, :] = deepcopy(prim[1, :])
        mixprim[2, 1] =
            prim[2, 1] +
            tau[1] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[2, 2] - prim[2, 1])
        mixprim[2, 2] =
            prim[2, 2] +
            tau[2] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[2, 1] - prim[2, 2])
        mixprim[3, 1] =
            prim[3, 1] +
            tau[1] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[3, 2] - prim[3, 1])
        mixprim[3, 2] =
            prim[3, 2] +
            tau[2] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[3, 1] - prim[3, 2])
        mixprim[4, 1] =
            prim[4, 1] +
            tau[1] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[4, 2] - prim[4, 1])
        mixprim[4, 2] =
            prim[4, 2] +
            tau[2] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[4, 1] - prim[4, 2])
        mixprim[5, 1] =
            1.0 / (
                1.0 / prim[end, 1] - 2.0 / 3.0 * (mixprim[2, 1] - prim[2, 1])^2 -
                2.0 / 3.0 * (mixprim[3, 1] - prim[3, 1])^2 -
                2.0 / 3.0 * (mixprim[4, 1] - prim[4, 1])^2 +
                tau[1] / kn * 2.0 * mi / (mi + me) *
                (
                    4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                    sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
                ) *
                (
                    1.0 / prim[end, 2] * me / mi - 1.0 / prim[end, 1] +
                    2.0 / 3.0 * me / mi * (prim[2, 2] - prim[2, 1])^2 +
                    2.0 / 3.0 * me / mi * (prim[3, 2] - prim[3, 1])^2 +
                    2.0 / 3.0 * me / mi * (prim[4, 2] - prim[4, 1])^2
                )
            )
        mixprim[5, 2] =
            1.0 / (
                1.0 / prim[end, 2] - 2.0 / 3.0 * (mixprim[2, 2] - prim[2, 2])^2 -
                2.0 / 3.0 * (mixprim[3, 2] - prim[3, 2])^2 -
                2.0 / 3.0 * (mixprim[4, 2] - prim[4, 2])^2 +
                tau[2] / kn * 2.0 * me / (mi + me) *
                (
                    4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                    sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
                ) *
                (
                    1.0 / prim[end, 1] * mi / me - 1.0 / prim[end, 2] +
                    2.0 / 3.0 * mi / me * (prim[2, 1] - prim[2, 2])^2 +
                    2.0 / 3.0 * mi / me * (prim[3, 1] - prim[3, 2])^2 +
                    2.0 / 3.0 * mi / me * (prim[4, 1] - prim[4, 2])^2
                )
            )
    else
        println("AAP mixture : dimension error")
    end

    return mixprim

end


# ------------------------------------------------------------
# Mixture source term function for DifferentialEquations.jl
# ------------------------------------------------------------
function aap_hs_diffeq(du, u, p, t)

    I₁, I₂, I₃, I₄, I₅, E₁, E₂, E₃, E₄, E₅ = u
    τᵢ, τₑ, mi, ni, me, ne, kn, γ = p

    τ = [τᵢ, τₑ]
    w = [
        I₁ E₁
        I₂ E₂
        I₃ E₃
        I₄ E₄
        I₅ E₅
    ]

    # modified variables
    prim = mixture_conserve_prim(w, γ)
    mixprim = aap_hs_prim(prim, τ, mi, ni, me, ne, kn)
    mixw = mixture_conserve_prim(prim, γ)

    du[1] = (mixw[1, 1] - I₁) / τᵢ
    du[2] = (mixw[2, 1] - I₂) / τᵢ
    du[3] = (mixw[3, 1] - I₃) / τᵢ
    du[4] = (mixw[4, 1] - I₄) / τᵢ
    du[5] = (mixw[5, 1] - I₅) / τᵢ
    du[6] = (mixw[1, 2] - E₁) / τₑ
    du[7] = (mixw[2, 2] - E₂) / τₑ
    du[8] = (mixw[3, 2] - E₃) / τₑ
    du[9] = (mixw[4, 2] - E₄) / τₑ
    du[10] = (mixw[5, 2] - E₅) / τₑ

    nothing

end


"""
Shift distribution function by external force

"""

# ------------------------------------------------------------
# Single-component gas
# ------------------------------------------------------------
function shift_pdf!(f::AbstractArray{<:Real,1}, a::Real, du::Real, dt::Real)

    q0 = eachindex(f) |> first # for OffsetArray
    q1 = eachindex(f) |> last

    if a > 0
        shift = Int(floor(a * dt / du)) # only for uniform velocity grid
        for k = q1:-1:q0+shift
            f[k] = f[k-shift]
        end
        for k = q0:shift+q0-1
            f[k] = 0.0
        end

        for k = q0+1:q1
            f[k] += (dt * a - du * shift) * (f[k-1] - f[k]) / du
        end
    else
        shift = Int(floor(-a * dt / du))
        for k = q0:q1-shift
            f[k] = f[k+shift]
        end
        for k = q1-shift+1:q1
            f[k] = 0.0
        end

        for k = q0:q1-1
            f[k] += (dt * a + du * shift) * (f[k] - f[k+1]) / du
        end
    end

    f[q0] = f[q0+1]
    f[q1] = f[q1-1]

end


# ------------------------------------------------------------
# Multi-component gas
# ------------------------------------------------------------
function shift_pdf!(
    f::AbstractArray{<:Real,2},
    a::Array{<:Real,1},
    du::AbstractArray{<:Real,1},
    dt::Real,
)
    for j in axes(f, 2)
        shift_pdf!(f[:, j], a[j], du[j], dt)
    end
end


function em_coefficients(
    prim::Array{<:Real,2},
    E::Array{<:Real,1},
    B::Array{<:Real,1},
    mr::Real,
    lD::Real,
    rL::Real,
    dt::Real,
)

    A = zeros(9, 9)
    A[1, 1] = -1.0 / (2.0 * rL)
    A[2, 2] = -1.0 / (2.0 * rL)
    A[3, 3] = -1.0 / (2.0 * rL)
    A[4, 1] = mr / (2.0 * rL)
    A[5, 2] = mr / (2.0 * rL)
    A[6, 3] = mr / (2.0 * rL)
    A[7, 1] = 1.0 / (dt)
    A[8, 2] = 1.0 / (dt)
    A[9, 3] = 1.0 / (dt)

    A[1, 4] = 1.0 / (dt)
    A[1, 5] = -B[3] / (2.0 * rL)
    A[1, 6] = B[2] / (2.0 * rL)
    A[2, 4] = B[3] / (2.0 * rL)
    A[2, 5] = 1.0 / (dt)
    A[2, 6] = -B[1] / (2.0 * rL)
    A[3, 4] = -B[2] / (2.0 * rL)
    A[3, 5] = B[1] / (2.0 * rL)
    A[3, 6] = 1.0 / (dt)

    A[4, 7] = 1.0 / (dt)
    A[4, 8] = mr * B[3] / (2.0 * rL)
    A[4, 9] = -mr * B[2] / (2.0 * rL)
    A[5, 7] = -mr * B[3] / (2.0 * rL)
    A[5, 8] = 1.0 / (dt)
    A[5, 9] = mr * B[1] / (2.0 * rL)
    A[6, 7] = mr * B[2] / (2.0 * rL)
    A[6, 8] = -mr * B[1] / (2.0 * rL)
    A[6, 9] = 1.0 / (dt)

    A[7, 4] = prim[1, 1] / (2.0 * rL * lD^2)
    A[8, 5] = prim[1, 1] / (2.0 * rL * lD^2)
    A[9, 6] = prim[1, 1] / (2.0 * rL * lD^2)
    A[7, 7] = -(prim[1, 2] * mr) / (2.0 * rL * lD^2)
    A[8, 8] = -(prim[1, 2] * mr) / (2.0 * rL * lD^2)
    A[9, 9] = -(prim[1, 2] * mr) / (2.0 * rL * lD^2)

    b = zeros(9)
    b[1] =
        prim[2, 1] / (dt) + E[1] / (2.0 * rL) - B[2] * prim[4, 1] / (2.0 * rL) +
        B[3] * prim[3, 1] / (2.0 * rL)
    b[2] =
        prim[3, 1] / (dt) + E[2] / (2.0 * rL) - B[3] * prim[2, 1] / (2.0 * rL) +
        B[1] * prim[4, 1] / (2.0 * rL)
    b[3] =
        prim[4, 1] / (dt) + E[3] / (2.0 * rL) - B[1] * prim[3, 1] / (2.0 * rL) +
        B[2] * prim[2, 1] / (2.0 * rL)
    b[4] =
        prim[2, 2] / (dt) - mr * E[1] / (2.0 * rL) + mr * B[2] * prim[4, 2] / (2.0 * rL) -
        mr * B[3] * prim[3, 2] / (2.0 * rL)
    b[5] =
        prim[3, 2] / (dt) - mr * E[2] / (2.0 * rL) + mr * B[3] * prim[2, 2] / (2.0 * rL) -
        mr * B[1] * prim[4, 2] / (2.0 * rL)
    b[6] =
        prim[4, 2] / (dt) - mr * E[3] / (2.0 * rL) + mr * B[1] * prim[3, 2] / (2.0 * rL) -
        mr * B[2] * prim[2, 2] / (2.0 * rL)
    b[7] =
        E[1] / (dt) - prim[1, 1] * prim[2, 1] / (2.0 * rL * lD^2) +
        prim[1, 2] * prim[2, 2] * mr / (2.0 * rL * lD^2)
    b[8] =
        E[2] / (dt) - prim[1, 1] * prim[3, 1] / (2.0 * rL * lD^2) +
        prim[1, 2] * prim[3, 2] * mr / (2.0 * rL * lD^2)
    b[9] =
        E[3] / (dt) - prim[1, 1] * prim[4, 1] / (2.0 * rL * lD^2) +
        prim[1, 2] * prim[4, 2] * mr / (2.0 * rL * lD^2)

    return A, b

end
