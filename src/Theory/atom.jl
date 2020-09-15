# ============================================================
# Atomistic Theory
# ============================================================

"""
Calculate moments of Gaussian distribution G = (λ / π)^(D / 2) * exp[-λ(c^2 + ξ^2)]

* with internality: `gauss_moments(prim::AbstractArray{<:Real,1}, inK::Real)`
* without internality: `gauss_moments(prim::AbstractArray{<:Real,1})`

"""
function gauss_moments(prim::AbstractArray{<:Real,1})

    if eltype(prim) <: Int
        MuL = OffsetArray{Float64}(undef, 0:6)
    else
        MuL = OffsetArray{eltype(prim)}(undef, 0:6)
    end
    MuR = similar(MuL)
    Mu = similar(MuL)

    MuL[0] = 0.5 * SpecialFunctions.erfc(-sqrt(prim[end]) * prim[2])
    MuL[1] = prim[2] * MuL[0] + 0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])
    MuR[0] = 0.5 * SpecialFunctions.erfc(sqrt(prim[end]) * prim[2])
    MuR[1] = prim[2] * MuR[0] - 0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])

    for i = 2:6
        MuL[i] = prim[2] * MuL[i-1] + 0.5 * (i - 1) * MuL[i-2] / prim[end]
        MuR[i] = prim[2] * MuR[i-1] + 0.5 * (i - 1) * MuR[i-2] / prim[end]
    end

    @. Mu = MuL + MuR

    if length(prim) == 3

        return Mu, MuL, MuR

    elseif length(prim) == 4

        Mv = similar(MuL)
        Mv[0] = 1.0
        Mv[1] = prim[3]
        for i = 2:6
            Mv[i] = prim[3] * Mv[i-1] + 0.5 * (i - 1) * Mv[i-2] / prim[end]
        end

        return Mu, Mv, MuL, MuR

    elseif length(prim) == 5

        Mv = similar(MuL)
        Mv[0] = 1.0
        Mv[1] = prim[3]
        for i = 2:6
            Mv[i] = prim[3] * Mv[i-1] + 0.5 * (i - 1) * Mv[i-2] / prim[end]
        end

        Mw = similar(MuL)
        Mw[0] = 1.0
        Mw[1] = prim[4]
        for i = 2:6
            Mw[i] = prim[4] * Mw[i-1] + 0.5 * (i - 1) * Mw[i-2] / prim[end]
        end

        return Mu, Mv, Mw, MuL, MuR

    end

end


function gauss_moments(prim::AbstractArray{<:Real,1}, inK::Real)

    if eltype(prim) <: Int
        MuL = OffsetArray{Float64}(undef, 0:6)
    else
        MuL = OffsetArray{eltype(prim)}(undef, 0:6)
    end
    MuR = similar(MuL)
    Mu = similar(MuL)

    MuL[0] = 0.5 * SpecialFunctions.erfc(-sqrt(prim[end]) * prim[2])
    MuL[1] = prim[2] * MuL[0] + 0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])
    MuR[0] = 0.5 * SpecialFunctions.erfc(sqrt(prim[end]) * prim[2])
    MuR[1] = prim[2] * MuR[0] - 0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])

    for i = 2:6
        MuL[i] = prim[2] * MuL[i-1] + 0.5 * (i - 1) * MuL[i-2] / prim[end]
        MuR[i] = prim[2] * MuR[i-1] + 0.5 * (i - 1) * MuR[i-2] / prim[end]
    end

    @. Mu = MuL + MuR

    if length(prim) == 3

        Mxi = similar(MuL, 0:2)
        Mxi[0] = 1.0
        Mxi[1] = 0.5 * inK / prim[end]
        Mxi[2] = (inK^2 + 2.0 * inK) / (4.0 * prim[end]^2)

        return Mu, Mxi, MuL, MuR

    elseif length(prim) == 4

        Mv = similar(MuL)
        Mv[0] = 1.0
        Mv[1] = prim[3]
        for i = 2:6
            Mv[i] = prim[3] * Mv[i-1] + 0.5 * (i - 1) * Mv[i-2] / prim[end]
        end

        Mxi = similar(MuL, 0:2)
        Mxi[0] = 1.0
        Mxi[1] = 0.5 * inK / prim[end]
        Mxi[2] = (inK^2 + 2.0 * inK) / (4.0 * prim[end]^2)

        return Mu, Mv, Mxi, MuL, MuR

    elseif length(prim) == 5

        Mv = similar(MuL)
        Mv[0] = 1.0
        Mv[1] = prim[3]
        for i = 2:6
            Mv[i] = prim[3] * Mv[i-1] + 0.5 * (i - 1) * Mv[i-2] / prim[end]
        end

        Mw = similar(MuL)
        Mw[0] = 1.0
        Mw[1] = prim[4]
        for i = 2:6
            Mw[i] = prim[4] * Mw[i-1] + 0.5 * (i - 1) * Mw[i-2] / prim[end]
        end

        return Mu, Mv, Mw, MuL, MuR

    end

end


"""
Calculate moments of Gaussian distribution in multi-component gas

`mixture_gauss_moments(prim::AbstractArray{<:Real,2}, inK::Real)`

"""
function mixture_gauss_moments(prim::AbstractArray{<:Real,2}, inK::Real)

    if eltype(prim) <: Int
        Mu = OffsetArray{Float64}(undef, 0:6, axes(prim, 2))
    else
        Mu = OffsetArray{eltype(prim)}(undef, 0:6, axes(prim, 2))
    end

    MuL = similar(Mu)
    MuR = similar(Mu)

    if size(prim, 1) == 3

        Mxi = similar(Mu, 0:2, axes(prim, 2))
        for j in axes(prim, 2)
            _tu, _txi, _tuL, _tuR = gauss_moments(prim[:, j], inK)

            Mu[:, j] .= _tu
            Mxi[:, j] .= _txi
            MuL[:, j] .= _tuL
            MuR[:, j] .= _tuR
        end

        return Mu, Mxi, MuL, MuR

    elseif size(prim, 1) == 4

        Mv = similar(Mu)
        Mxi = similar(Mu, 0:2, axes(prim, 2))
        for j in axes(prim, 2)
            _tu, _tv, _txi, _tuL, _tuR = gauss_moments(prim[:, j], inK)

            Mu[:, j] .= _tu
            Mv[:, j] .= _tv
            Mxi[:, j] .= _txi
            MuL[:, j] .= _tuL
            MuR[:, j] .= _tuR
        end

        return Mu, Mv, Mxi, MuL, MuR

    elseif size(prim, 1) == 5

        Mv = similar(Mu)
        Mw = similar(Mu)

        for j in axes(prim, 2)
            _tu, _tv, _tw, _tuL, _tuR = gauss_moments(prim[:, j], inK)

            Mu[:, j] .= _tu
            Mv[:, j] .= _tv
            Mw[:, j] .= _tw
            MuL[:, j] .= _tuL
            MuR[:, j] .= _tuR
        end

        return Mu, Mv, Mw, MuL, MuR

    end

end


"""
Calculate conservative moments of particle distribution

`moments_conserve(Mu::OffsetArray{<:Real,1}, Mxi::OffsetArray{<:Real,1}, 
    alpha::Int, delta::Int)`
`moments_conserve(Mu::OffsetArray{<:Real,1}, Mv::OffsetArray{<:Real,1},
    Mw::OffsetArray{<:Real,1}, alpha::Int, beta::Int, delta::Int)`

"""
moments_conserve(Mu::OffsetArray{<:AbstractFloat,1}, alpha::Int) = Mu[alpha]


function moments_conserve(
    Mu::OffsetArray{<:AbstractFloat,1},
    Mxi::OffsetArray{<:AbstractFloat,1},
    alpha::Int,
    delta::Int,
)

    uv = zeros(eltype(Mu), 3)
    uv[1] = Mu[alpha] * Mxi[delta÷2]
    uv[2] = Mu[alpha+1] * Mxi[delta÷2]
    uv[3] = 0.5 * (Mu[alpha+2] * Mxi[delta÷2] + Mu[alpha] * Mxi[(delta+2)÷2])

    return uv

end


function moments_conserve(
    Mu::OffsetArray{<:AbstractFloat,1},
    Mv::OffsetArray{<:AbstractFloat,1},
    Mw::OffsetArray{<:AbstractFloat,1},
    alpha::Int,
    beta::Int,
    delta::Int,
)

    if length(Mw) == 3 # internal motion

        uv = zeros(eltype(Mu), 4)
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

        uv = zeros(eltype(Mu), 5)
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


"""
Calculate conservative moments of particle distribution in multi-component gas

"""
function mixture_moments_conserve(
    Mu::OffsetArray{<:AbstractFloat,2},
    Mxi::OffsetArray{<:AbstractFloat,2},
    alpha::Int,
    delta::Int,
)

    Muv = zeros(eltype(Mu), 3, size(Mu, 2))
    for j in axes(Muv, 2)
        Muv[:, j] .= moments_conserve(Mu[:, j], Mxi[:, j], alpha, delta)
    end

    return Muv

end


function mixture_moments_conserve(
    Mu::OffsetArray{<:AbstractFloat,2},
    Mv::OffsetArray{<:AbstractFloat,2},
    Mw::OffsetArray{<:AbstractFloat,2},
    alpha::Int,
    beta::Int,
    delta::Int,
)

    Muv = ifelse(
        length(Mw) == 3,
        zeros(eltype(Mu), 4, size(Mu, 2)),
        zeros(eltype(Mu), 5, size(Mu, 2)),
    )
    for j in axes(Muv, 2)
        Muv[:, j] .= moments_conserve(Mu[:, j], Mv[:, j], Mw[:, j], alpha, beta, delta)
    end

    return Muv

end


"""
Calculate slope of particle distribution function, 
assuming a = a1 + u * a2 + 0.5 * u^2 * a3

`pdf_slope(u::Real, Δ::Real)`

`pdf_slope(prim::AbstractArray{<:Real,1}, sw::AbstractArray{<:Real,1}, inK::Real)`

"""
pdf_slope(u::Real, Δ::Real) = Δ / u


function pdf_slope(prim::AbstractArray{<:Real,1}, sw::AbstractArray{<:Real,1}, inK::Real)

    sl = zeros(eltype(sw), axes(prim))

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


"""
Calculate slope of multi-component particle distribution function, 
assuming `a = a1 + u * a2 + 0.5 * u^2 * a3`

`mixture_pdf_slope(prim::AbstractArray{<:Real,2}, sw::AbstractArray{<:Real,2}, inK::Real)`

"""
function mixture_pdf_slope(prim::AbstractArray{<:Real,2}, sw::AbstractArray{<:Real,2}, inK::Real)

    sl = zeros(eltype(sw), axes(prim))
    for j in axes(sl, 2)
        sl[:, j] .= pdf_slope(prim[:, j], sw[:, j], inK)
    end

    return sl

end


"""
Calculate slope-related conservative moments,
assuming a = a1 + u * a2 + 0.5 * u^2 * a3

"""
moments_conserve_slope(a::Real, Mu::OffsetArray{<:Real,1}, alpha::Int) = 
    a * moments_conserve(Mu, alpha)


moments_conserve_slope(
    a::AbstractArray{<:Real,1},
    Mu::OffsetArray{<:Real,1},
    Mxi::OffsetArray{<:Real,1},
    alpha::Int,
) =
    a[1] .* moments_conserve(Mu, Mxi, alpha + 0, 0) .+
    a[2] .* moments_conserve(Mu, Mxi, alpha + 1, 0) .+
    0.5 * a[3] .* moments_conserve(Mu, Mxi, alpha + 2, 0) .+
    0.5 * a[3] .* moments_conserve(Mu, Mxi, alpha + 0, 2)


function moments_conserve_slope(
    a::AbstractArray{<:Real,1},
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
    a::AbstractArray{<:Real,1},
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
    a::AbstractArray{<:Real,2},
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
    a::AbstractArray{<:Real,2},
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
    a::AbstractArray{<:Real,2},
    Mu::OffsetArray{<:Real,2},
    Mv::OffsetArray{<:Real,2},
    Mw::OffsetArray{<:Real,2},
    alpha::Int,
    beta::Int,
    delta::Int,
)

    au = zeros(eltype(a), 5, axes(a, 2))
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
2) discrete form

"""

"""
Discrete moments of particle distribution

* `discrete_moments(f, ω)` : direct quadrature
* `discrete_moments(f, u, ω, n)` : velocity moments

"""
discrete_moments(f::AbstractArray{<:AbstractFloat,1}, ω::AbstractArray{<:AbstractFloat,1}) =
    sum(@. ω * f)


#--- 1D DVM ---#
discrete_moments(
    f::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    ω::AbstractArray{<:AbstractFloat,1},
    n::Int,
) = sum(@. ω * u^n * f)


#--- 2D DVM ---#
discrete_moments(
    f::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
    n::Int,
) = sum(@. ω * u^n * f)


#--- 3D DVM ---#
discrete_moments(
    f::AbstractArray{<:AbstractFloat,3},
    u::AbstractArray{<:AbstractFloat,3},
    ω::AbstractArray{<:AbstractFloat,3},
    n::Int,
) = sum(@. ω * u^n * f)


# ------------------------------------------------------------
# Conservative moments
# ------------------------------------------------------------
#--- 1D ---#
function moments_conserve(
    f::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    ω::AbstractArray{<:AbstractFloat,1},
)
    w = zeros(eltype(f), 3)
    w[1] = discrete_moments(f, u, ω, 0)
    w[2] = discrete_moments(f, u, ω, 1)
    w[3] = 0.5 * discrete_moments(f, u, ω, 2)

    return w
end


function moments_conserve(
    h::AbstractArray{<:AbstractFloat,1},
    b::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    ω::AbstractArray{<:AbstractFloat,1},
)
    w = zeros(eltype(f), 3)
    w[1] = discrete_moments(h, u, ω, 0)
    w[2] = discrete_moments(h, u, ω, 1)
    w[3] = 0.5 * (discrete_moments(h, u, ω, 2) + discrete_moments(b, u, ω, 0))

    return w
end


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


#--- 2D ---#
function moments_conserve(
    f::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
)
    w = zeros(eltype(f), 4)
    w[1] = discrete_moments(f, u, ω, 0)
    w[2] = discrete_moments(f, u, ω, 1)
    w[3] = discrete_moments(f, v, ω, 1)
    w[4] = 0.5 * (discrete_moments(f, u, ω, 2) + discrete_moments(f, v, ω, 2))

    return w
end

function moments_conserve(
    h::AbstractArray{<:AbstractFloat,2},
    b::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
)
    w = zeros(eltype(h), 4)
    w[1] = discrete_moments(h, u, ω, 0)
    w[2] = discrete_moments(h, u, ω, 1)
    w[3] = discrete_moments(h, v, ω, 1)
    w[4] = 0.5 * (
        discrete_moments(h, u, ω, 2) +
        discrete_moments(h, v, ω, 2) +
        discrete_moments(b, u, ω, 0)
    )

    return w
end


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


#--- 3D ---#
function moments_conserve(
    f::AbstractArray{<:AbstractFloat,3},
    u::AbstractArray{<:AbstractFloat,3},
    v::AbstractArray{<:AbstractFloat,3},
    w::AbstractArray{<:AbstractFloat,3},
    ω::AbstractArray{<:AbstractFloat,3},
)
    moments = zeros(eltype(f), 5)

    moments[1] = discrete_moments(f, u, ω, 0)
    moments[2] = discrete_moments(f, u, ω, 1)
    moments[3] = discrete_moments(f, v, ω, 1)
    moments[4] = discrete_moments(f, w, ω, 1)
    moments[5] = 0.5 * (
        discrete_moments(f, u, ω, 2) +
        discrete_moments(f, v, ω, 2) +
        discrete_moments(f, w, ω, 2)
    )

    return moments
end


function moments_conserve(
    h0::AbstractArray{<:AbstractFloat,1},
    h1::AbstractArray{<:AbstractFloat,1},
    h2::AbstractArray{<:AbstractFloat,1},
    h3::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    ω::AbstractArray{<:AbstractFloat,1},
)
    moments = zeros(eltype(h0), 5)

    moments[1] = discrete_moments(h0, u, ω, 0)
    moments[2] = discrete_moments(h0, u, ω, 1)
    moments[3] = discrete_moments(h1, u, ω, 0)
    moments[4] = discrete_moments(h2, u, ω, 0)
    moments[5] = 0.5 * discrete_moments(h0, u, ω, 2) + 0.5 * discrete_moments(h3, u, ω, 0)

    return w
end


function mixture_moments_conserve(
    f::AbstractArray{<:AbstractFloat,4},
    u::AbstractArray{<:AbstractFloat,4},
    v::AbstractArray{<:AbstractFloat,4},
    w::AbstractArray{<:AbstractFloat,4},
    ω::AbstractArray{<:AbstractFloat,4},
)
    moments = zeros(eltype(f), 5, size(f, 4))

    for j in axes(w, 2)
        moments[:, j] .= moments_conserve(
            f[:, :, :, j],
            u[:, :, :, j],
            v[:, :, :, j],
            w[:, :, :, j],
            ω[:, :, :, j],
        )
    end

    return moments
end


function mixture_moments_conserve(
    h0::AbstractArray{<:AbstractFloat,2},
    h1::AbstractArray{<:AbstractFloat,2},
    h2::AbstractArray{<:AbstractFloat,2},
    h3::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
)
    moments = zeros(eltype(h0), 5, size(h0, 2))
    for j in axes(w, 2)
        moments[:, j] .=
            moments_conserve(h0[:, j], h1[:, j], h2[:, j], h3[:, j], u[:, j], ω[:, j])
    end

    return moments
end


"""
Calculate stress tensor from particle distribution function

"""
function stress(
    f::AbstractArray{<:AbstractFloat,1},
    prim::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    ω::AbstractArray{<:AbstractFloat,1},
)
    return sum(@. ω * (u - prim[2]) * (u - prim[2]) * f)
end


function stress(
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
Evaluate heat flux from particle distribution function

"""
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


#--- 2D ---#
function heat_flux(
    h::AbstractArray{<:AbstractFloat,2},
    prim::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
)
    q = zeros(eltype(f), 2)

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
    q = zeros(eltype(f), 2)

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


#--- 3D ---#
function heat_flux(
    f::AbstractArray{<:AbstractFloat,3},
    prim::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,3},
    v::AbstractArray{<:AbstractFloat,3},
    w::AbstractArray{<:AbstractFloat,3},
    ω::AbstractArray{<:AbstractFloat,3},
)
    q = zeros(eltype(f), 3)

    q[1] =
        0.5 * sum(@. ω *
               (u - prim[2]) *
               ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) *
               f)
    q[2] =
        0.5 * sum(@. ω *
               (v - prim[3]) *
               ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) *
               f)
    q[3] =
        0.5 * sum(@. ω *
               (w - prim[4]) *
               ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) *
               f)

    return q
end


"""
Maxwellian in discrete form

>@param[in] : particle velocity quadrature points
>@param[in] : density, velocity and inverse of temperature
>@return : Maxwellian distribution function

"""

#--- 1D ---#
maxwellian(u::AbstractArray{<:AbstractFloat,1}, ρ::Real, U::Real, λ::Real) =
    @. ρ * sqrt(λ / π) * exp(-λ * (u - U)^2)


maxwellian(u::AbstractArray{<:AbstractFloat,1}, prim::AbstractArray{<:Real,1}) =
    maxwellian(u, prim[1], prim[2], prim[end]) # in case of input with length 4/5


function mixture_maxwellian(u::AbstractArray{<:AbstractFloat,2}, prim::AbstractArray{<:Real,2})
    mixM = similar(u)
    for j in axes(mixM, 2)
        mixM[:, j] .= maxwellian(u[:, j], prim[:, j])
    end

    return mixM
end


#--- 2D ---#
maxwellian(
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    ρ::Real,
    U::Real,
    V::Real,
    λ::Real,
) = @. ρ * (λ / π) * exp(-λ * ((u - U)^2 + (v - V)^2))


maxwellian(
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    prim::AbstractArray{<:Real,1},
) = maxwellian(u, v, prim[1], prim[2], prim[3], prim[end]) # in case of input with length 5


function mixture_maxwellian(
    u::AbstractArray{<:AbstractFloat,3},
    v::AbstractArray{<:AbstractFloat,3},
    prim::AbstractArray{<:Real,2},
)
    mixM = zeros(axes(u))
    for k in axes(mixM, 3)
        mixM[:, :, k] .= maxwellian(u[:, :, k], v[:, :, k], prim[:, k])
    end

    return mixM
end


#--- 3D ---#
maxwellian(
    u::AbstractArray{<:AbstractFloat,3},
    v::AbstractArray{<:AbstractFloat,3},
    w::AbstractArray{<:AbstractFloat,3},
    ρ::Real,
    U::Real,
    V::Real,
    W::Real,
    λ::Real,
) = @. ρ * sqrt((λ / π)^3) * exp(-λ * ((u - U)^2 + (v - V)^2 + (w - W)^2))


maxwellian(
    u::AbstractArray{<:AbstractFloat,3},
    v::AbstractArray{<:AbstractFloat,3},
    w::AbstractArray{<:AbstractFloat,3},
    prim::AbstractArray{<:Real,1},
) = maxwellian(u, v, w, prim[1], prim[2], prim[3], prim[4], prim[5])


function mixture_maxwellian(
    u::AbstractArray{<:AbstractFloat,4},
    v::AbstractArray{<:AbstractFloat,4},
    w::AbstractArray{<:AbstractFloat,4},
    prim::AbstractArray{<:Real,2},
)
    mixM = zeros(eltype(u), axes(u))
    for l in axes(mixM, 4)
        mixM[:, :, :, l] .=
            maxwellian(u[:, :, :, l], v[:, :, :, l], w[:, :, :, l], prim[:, l])
    end

    return mixM
end


"""
Shakhov non-equilibrium part

> @param[in] : particle velocity quadrature points
> @param[in] : discrete Maxwellian
> @param[in] : primitive variables, Prandtl number, heat flux, inner degree of freedom

"""
#--- 1F1V ---#
function shakhov(
    u::AbstractArray{<:AbstractFloat,1},
    M::AbstractArray{<:AbstractFloat,1},
    q::Real,
    prim::AbstractArray{<:Real,1},
    Pr::Real,
)

    M_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       (u - prim[2]) *
       q *
       (2.0 * prim[end] * (u - prim[2])^2 - 5.0) *
       M

    return M_plus

end


#--- 2F1V ---#
function shakhov(
    u::AbstractArray{<:AbstractFloat,1},
    H::AbstractArray{<:AbstractFloat,1},
    B::AbstractArray{<:AbstractFloat,1},
    q::Real,
    prim::AbstractArray{<:Real,1},
    Pr::Real,
    K::Real,
)

    H_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       (u - prim[2]) *
       q *
       (2.0 * prim[end] * (u - prim[2])^2 + K - 5.0) *
       H
    B_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       (u - prim[2]) *
       q *
       (2.0 * prim[end] * (u - prim[2])^2 + K - 3.0) *
       B

    return H_plus, B_plus

end


#--- 1F2V ---#
function shakhov(
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    M::AbstractArray{<:AbstractFloat,2},
    q::AbstractArray{<:AbstractFloat,1},
    prim::AbstractArray{<:Real,1},
    Pr::Real,
)

    M_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       ((u - prim[2]) * q[1] + (v - prim[3]) * q[2]) *
       (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2) - 5.0) *
       M

    return M_plus

end


#--- 2F2V ---#
function shakhov(
    u::AbstractArray{<:AbstractFloat,2},
    v::AbstractArray{<:AbstractFloat,2},
    H::AbstractArray{<:AbstractFloat,2},
    B::AbstractArray{<:AbstractFloat,2},
    q::AbstractArray{<:AbstractFloat,1},
    prim::AbstractArray{<:Real,1},
    Pr::Real,
    K::Real,
)

    H_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       ((u - prim[2]) * q[1] + (v - prim[3]) * q[2]) *
       (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 5.0) *
       H
    B_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       ((u - prim[2]) * q[1] + (v - prim[3]) * q[2]) *
       (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 3.0) *
       B

    return H_plus, B_plus

end


#--- 1F3V ---#
function shakhov(
    u::AbstractArray{<:AbstractFloat,3},
    v::AbstractArray{<:AbstractFloat,3},
    w::AbstractArray{<:AbstractFloat,3},
    M::AbstractArray{<:AbstractFloat,3},
    q::AbstractArray{<:AbstractFloat,1},
    prim::AbstractArray{<:Real,1},
    Pr::Real,
)

    M_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       ((u - prim[2]) * q[1] + (v - prim[3]) * q[2] + (w - prim[4]) * q[3]) *
       (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) - 5.0) *
       M

    return M_plus

end


"""
Reduced distribution function

@arg : particle distribution function with full velocity space
@arg : quadrature weights with reduced velocity setting (v & w by default)

"""
function reduce_distribution(
    f::AbstractArray{<:AbstractFloat,2},
    weights::AbstractArray{<:AbstractFloat,1},
    dim = 1::Int,
)

    if dim == 1
        h = similar(f, axes(f, 1))
        for i in eachindex(h)
            h[i] = sum(@. weights * f[i, :])
        end
    elseif dim == 2
        h = similar(f, axes(f, 2))
        for j in eachindex(h)
            h[j] = sum(@. weights * f[:, j])
        end
    else
        throw("dimension dismatch")
    end

    return h

end

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
    prim::AbstractArray{<:Real,1},
    γ = 5 / 3::Real,
) = full_distribution(h, b, u, weights, v, w, prim[1], γ)


"""
Calculate reference viscosity
* variable hard sphere (VHS) model

"""
ref_vhs_vis(Kn::Real, alpha::Real, omega::Real) =
    5.0 * (alpha + 1.0) * (alpha + 2.0) * √π /
    (4.0 * alpha * (5.0 - 2.0 * omega) * (7.0 - 2.0 * omega)) * Kn


"""
Calculate collision time
* variable hard sphere (VHS) model

"""
vhs_collision_time(prim::Array{<:Real,1}, muRef::Real, omega::Real) =
    muRef * 2.0 * prim[end]^(1.0 - omega) / prim[1]


"""
# Calculate effective Knudsen number for fast spectral method
* hard sphere (HS) model
"""
hs_boltz_kn(mu_ref::Real, alpha::Real) =
    64 * sqrt(2.0)^alpha / 5.0 * gamma((alpha + 3) / 2) * gamma(2.0) * sqrt(pi) * mu_ref


"""
Calculate collision kernel for fast spectral method

"""
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


"""
Calculate collision operator with FFT-based fast spectral method

"""
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


function boltzmann_fft!(
    Q::AbstractArray{<:Real,3},
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

    @. Q = 4.0 * π^2 / Kn / M^2 * real(f_temp)

end


"""
Calculate mixture collision time from AAP model

"""
function aap_hs_collision_time(
    prim::AbstractArray{<:Real,2},
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


"""
Calculate mixture primitive variables from AAP model

"""
function aap_hs_prim(
    prim::AbstractArray{<:Real,2},
    tau::AbstractArray{<:Real,1},
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

        throw("AAP mixture : dimension dismatch")

    end

    return mixprim

end


"""
Source term of AAP model in DifferentialEquations.jl

"""
function aap_hs_diffeq!(du, u, p, t)

    if length(u) == 6
        I₁, I₂, I₃, E₁, E₂, E₃ = u
        w = [
            I₁ E₁
            I₂ E₂
            I₃ E₃
        ]
    elseif length(u) == 8
        I₁, I₂, I₃, I₄, E₁, E₂, E₃, E₄ = u
        w = [
            I₁ E₁
            I₂ E₂
            I₃ E₃
            I₄ E₄
        ]
    elseif length(u) == 10
        I₁, I₂, I₃, I₄, I₅, E₁, E₂, E₃, E₄, E₅ = u
        w = [
            I₁ E₁
            I₂ E₂
            I₃ E₃
            I₄ E₄
            I₅ E₅
        ]
    else
    end

    τᵢ, τₑ, mi, ni, me, ne, kn, γ = p
    τ = [τᵢ, τₑ]

    # modified variables
    prim = mixture_conserve_prim(w, γ)
    mixprim = aap_hs_prim(prim, τ, mi, ni, me, ne, kn)
    mixw = mixture_conserve_prim(mixprim, γ)

    if length(u) == 6
        du[1] = (mixw[1, 1] - I₁) / τᵢ
        du[2] = (mixw[2, 1] - I₂) / τᵢ
        du[3] = (mixw[3, 1] - I₃) / τᵢ
        du[4] = (mixw[1, 2] - E₁) / τₑ
        du[5] = (mixw[2, 2] - E₂) / τₑ
        du[6] = (mixw[3, 2] - E₃) / τₑ
    elseif length(u) == 8
        du[1] = (mixw[1, 1] - I₁) / τᵢ
        du[2] = (mixw[2, 1] - I₂) / τᵢ
        du[3] = (mixw[3, 1] - I₃) / τᵢ
        du[4] = (mixw[4, 1] - I₄) / τᵢ
        du[5] = (mixw[1, 2] - E₁) / τₑ
        du[6] = (mixw[2, 2] - E₂) / τₑ
        du[7] = (mixw[3, 2] - E₃) / τₑ
        du[8] = (mixw[4, 2] - E₄) / τₑ
    elseif length(u) == 10
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
    end

    nothing

end


"""
Shift distribution function by external force

"""
function shift_pdf!(
    f::AbstractArray{<:AbstractFloat,1},
    a::Real,
    du::AbstractFloat,
    dt::Real,
)

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


#--- Multi-component gas ---#
function shift_pdf!(
    f::AbstractArray{<:AbstractFloat,2},
    a::AbstractArray{<:Real,1},
    du::AbstractArray{<:AbstractFloat,1},
    dt::Real,
)
    for j in axes(f, 2)
        _f = @view f[:, j]
        shift_pdf!(_f, a[j], du[j], dt)
    end
end
