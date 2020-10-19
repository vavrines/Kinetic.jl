"""
Calculate moments of Gaussian distribution in multi-component gas

    mixture_gauss_moments(prim::T, inK) where {T<:AbstractArray{<:Real,2}}

"""
function mixture_gauss_moments(prim::T, inK) where {T<:AbstractArray{<:Real,2}}

    if eltype(prim) <: Int
        Mu = OffsetArray(
            similar(prim, Float64, 7, axes(prim, 2)),
            0:6,
            axes(prim, 2),
        )
    else
        Mu = OffsetArray(similar(prim, 7, axes(prim, 2)), 0:6, axes(prim, 2))
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

    mixture_moments_conserve(
        Mu::T,
        Mxi::T,
        alpha::I,
        delta::I,
    ) where {T<:OffsetArray{<:AbstractFloat,2},I<:Int}

    function mixture_moments_conserve(
        Mu::T,
        Mv::T,
        Mw::T,
        alpha::I,
        beta::I,
        delta::I,
    ) where {T<:OffsetArray{<:AbstractFloat,2},I<:Int}

"""
function mixture_moments_conserve(
    Mu::T,
    Mxi::T,
    alpha::I,
    delta::I,
) where {T<:OffsetArray{<:AbstractFloat,2},I<:Int}

    Muv = similar(Mu, 3, size(Mu, 2))
    for j in axes(Muv, 2)
        Muv[:, j] .= moments_conserve(Mu[:, j], Mxi[:, j], alpha, delta)
    end

    return Muv

end

function mixture_moments_conserve(
    Mu::T,
    Mv::T,
    Mw::T,
    alpha::I,
    beta::I,
    delta::I,
) where {T<:OffsetArray{<:AbstractFloat,2},I<:Int}

    Muv = ifelse(
        size(Mw, 1) == 3,
        similar(Mu, 4, size(Mu, 2)),
        similar(Mu, 5, size(Mu, 2)),
    )
    for j in axes(Muv, 2)
        Muv[:, j] .=
            moments_conserve(Mu[:, j], Mv[:, j], Mw[:, j], alpha, beta, delta)
    end

    return Muv

end

# ------------------------------------------------------------
# Discrete moments of conservative variables
# ------------------------------------------------------------
#--- 1F1V ---#
function mixture_moments_conserve(
    f::X,
    u::T,
    ω::T,
) where {X<:AbstractArray{<:AbstractFloat,2},T<:AbstractArray{<:AbstractFloat,2}}

    w = similar(f, 3, size(f, 2))
    for j in axes(w, 2)
        w[:, j] .= moments_conserve(f[:, j], u, ω)
    end

    return w

end

#--- 2F1V ---#
function mixture_moments_conserve(
    h::X,
    b::X,
    u::T,
    ω::T,
) where {X<:AbstractArray{<:AbstractFloat,2},T<:AbstractArray{<:AbstractFloat,2}}

    w = similar(h, 3, size(h, 2))
    for j in axes(w, 2)
        w[:, j] .= moments_conserve(h[:, j], b[:, j], u, ω)
    end

    return w

end

#--- 4F1V ---#
function mixture_moments_conserve(
    h0::X,
    h1::X,
    h2::X,
    h3::X,
    u::T,
    ω::T,
) where {X<:AbstractArray{<:AbstractFloat,2},T<:AbstractArray{<:AbstractFloat,2}}

    moments = similar(h0, 5, size(h0, 2))
    for j in axes(moments, 2)
        moments[:, j] .= moments_conserve(
            h0[:, j],
            h1[:, j],
            h2[:, j],
            h3[:, j],
            u[:, j],
            ω[:, j],
        )
    end

    return moments

end

#--- 1F2V ---#
function mixture_moments_conserve(
    f::X,
    u::T,
    v::T,
    ω::T,
) where {X<:AbstractArray{<:AbstractFloat,3},T<:AbstractArray{<:AbstractFloat,3}}

    w = similar(f, 4, size(f, 3))
    for j in axes(w, 2)
        w[:, j] .=
            moments_conserve(f[:, :, j], u[:, :, j], v[:, :, j], ω[:, :, j])
    end

    return w

end

#--- 2F2V ---#
function mixture_moments_conserve(
    h::X,
    b::X,
    u::T,
    v::T,
    ω::T,
) where {X<:AbstractArray{<:AbstractFloat,3},T<:AbstractArray{<:AbstractFloat,3}}

    w = similar(h, 4, size(f, 3))
    for j in axes(w, 2)
        w[:, j] .= moments_conserve(
            h[:, :, j],
            b[:, :, j],
            u[:, :, j],
            v[:, :, j],
            ω[:, :, j],
        )
    end

    return w

end

#--- 3F2V ---#
function mixture_moments_conserve(
    h0::X,
    h1::X,
    h2::X,
    u::T,
    v::T,
    ω::T,
) where {X<:AbstractArray{<:AbstractFloat,3},T<:AbstractArray{<:AbstractFloat,3}}
    
    w = similar(h0, 5, size(h0, 3))
    for j in axes(w, 2)
        w[:, j] .= moments_conserve(
            h0[:, :, j],
            h1[:, :, j],
            h2[:, :, j],
            u[:, :, j],
            v[:, :, j],
            ω[:, :, j],
        )
    end

    return w

end

#--- 1F3V ---#
function mixture_moments_conserve(
    f::X,
    u::T,
    v::T,
    w::T,
    ω::T,
) where {X<:AbstractArray{<:AbstractFloat,4},T<:AbstractArray{<:AbstractFloat,4}}

    moments = similar(f, 5, size(f, 4))
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


"""
Calculate slope-related conservative moments
`a = a1 + u * a2 + 0.5 * u^2 * a3`

"""
function mixture_moments_conserve_slope(
    a::X,
    Mu::Y,
    Mxi::Y,
    alpha::I,
) where {X<:AbstractArray{<:Real,2},Y<:OffsetArray{<:Real,2},I<:Int}

    au = similar(a, 3, axes(a, 2))
    for j in axes(au, 2)
        au[:, j] .= moments_conserve_slope(a[:, j], Mu[:, j], Mxi[:, j], alpha)
    end

    return au

end

function mixture_moments_conserve_slope(
    a::X,
    Mu::Y,
    Mv::Y,
    Mxi::Y,
    alpha::I,
    beta::I,
) where {X<:AbstractArray{<:Real,2},Y<:OffsetArray{<:Real,2},I<:Int}

    au = similar(a, 4, axes(a, 2))
    for j in axes(au, 2)
        au[:, j] .= moments_conserve_slope(
            a[:, j],
            Mu[:, j],
            Mv[:, j],
            Mxi[:, j],
            alpha,
            beta,
        )
    end

    return au

end

function mixture_moments_conserve_slope(
    a::X,
    Mu::Y,
    Mv::Y,
    Mw::Y,
    alpha::I,
    beta::I,
    delta::I,
) where {X<:AbstractArray{<:Real,2},Y<:OffsetArray{<:Real,2},I<:Int}

    au = similar(a, 5, axes(a, 2))
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
