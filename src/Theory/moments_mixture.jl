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
        size(Mw, 1) == 3,
        zeros(eltype(Mu), 4, size(Mu, 2)),
        zeros(eltype(Mu), 5, size(Mu, 2)),
    )
    for j in axes(Muv, 2)
        Muv[:, j] .= moments_conserve(Mu[:, j], Mv[:, j], Mw[:, j], alpha, beta, delta)
    end

    return Muv

end


"""
Calculate slope-related conservative moments 
`a = a1 + u * a2 + 0.5 * u^2 * a3`

"""
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


# ------------------------------------------------------------
# Discrete moments of conservative variables
# ------------------------------------------------------------
#--- 1F1V ---#
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

#--- 2F1V ---#
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

#--- 4F1V ---#
function mixture_moments_conserve(
    h0::AbstractArray{<:AbstractFloat,2},
    h1::AbstractArray{<:AbstractFloat,2},
    h2::AbstractArray{<:AbstractFloat,2},
    h3::AbstractArray{<:AbstractFloat,2},
    u::AbstractArray{<:AbstractFloat,2},
    ω::AbstractArray{<:AbstractFloat,2},
)
    moments = zeros(eltype(h0), 5, size(h0, 2))
    for j in axes(moments, 2)
        moments[:, j] .=
            moments_conserve(h0[:, j], h1[:, j], h2[:, j], h3[:, j], u[:, j], ω[:, j])
    end

    return moments
end

#--- 1F2V ---#
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

#--- 2F2V ---#
function mixture_moments_conserve(
    h::AbstractArray{<:AbstractFloat,3},
    b::AbstractArray{<:AbstractFloat,3},
    u::AbstractArray{<:AbstractFloat,3},
    v::AbstractArray{<:AbstractFloat,3},
    ω::AbstractArray{<:AbstractFloat,3},
)
    w = zeros(eltype(h), 4, size(h, 3))
    for j in axes(w, 2)
        w[:, j] .=
            moments_conserve(h[:, :, j], b[:, :, j], u[:, :, j], v[:, :, j], ω[:, :, j])
    end

    return w
end

#--- 3F2V ---#
function mixture_moments_conserve(
    h0::AbstractArray{<:AbstractFloat,3},
    h1::AbstractArray{<:AbstractFloat,3},
    h2::AbstractArray{<:AbstractFloat,3},
    u::AbstractArray{<:AbstractFloat,3},
    v::AbstractArray{<:AbstractFloat,3},
    ω::AbstractArray{<:AbstractFloat,3},
)
    w = zeros(eltype(h), 5, size(h0, 3))
    for j in axes(w, 2)
        w[:, j] .=
            moments_conserve(
                h0[:, :, j], 
                h1[:, :, j], 
                h2[:, :, j], 
                u[:, :, j], 
                v[:, :, j], 
                ω[:, :, j]
            )
    end

    return w
end

#--- 1F3V ---#
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
