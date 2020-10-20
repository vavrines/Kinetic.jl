"""
Calculate moments of Gaussian distribution `G = (λ / π)^(D / 2) * exp[-λ(c^2 + ξ^2)]`

* internality: `gauss_moments(prim::T) where {T<:AbstractArray{<:Real,1}}`
* no internality: `gauss_moments(prim::T, inK) where {T<:AbstractArray{<:Real,1}}`

"""
function gauss_moments(prim::T) where {T<:AbstractArray{<:Real,1}}

    if eltype(prim) <: Int
        MuL = OffsetArray(similar(prim, Float64, 7), 0:6)
    else
        MuL = OffsetArray(similar(prim, 7), 0:6)
    end
    MuR = similar(MuL)
    Mu = similar(MuL)

    MuL[0] = 0.5 * SpecialFunctions.erfc(-sqrt(prim[end]) * prim[2])
    MuL[1] =
        prim[2] * MuL[0] +
        0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])
    MuR[0] = 0.5 * SpecialFunctions.erfc(sqrt(prim[end]) * prim[2])
    MuR[1] =
        prim[2] * MuR[0] -
        0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])
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

# ------------------------------------------------------------
# A more general function
# deal with absent internality by setting inK = 0
# ------------------------------------------------------------
function gauss_moments(prim::T, inK) where {T<:AbstractArray{<:Real,1}}

    if eltype(prim) <: Int
        MuL = OffsetArray(similar(prim, Float64, 7), 0:6)
    else
        MuL = OffsetArray(similar(prim, 7), 0:6)
    end
    MuR = similar(MuL)
    Mu = similar(MuL)

    MuL[0] = 0.5 * SpecialFunctions.erfc(-sqrt(prim[end]) * prim[2])
    MuL[1] =
        prim[2] * MuL[0] +
        0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])
    MuR[0] = 0.5 * SpecialFunctions.erfc(sqrt(prim[end]) * prim[2])
    MuR[1] =
        prim[2] * MuR[0] -
        0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])
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
Calculate conservative moments of particle distribution

`moments_conserve(Mu::OffsetArray{<:AbstractFloat,1}, alpha::Int)`

`moments_conserve(Mu::OffsetArray{<:Real,1}, Mxi::OffsetArray{<:Real,1},
    alpha::Int, delta::Int)`

`moments_conserve(Mu::OffsetArray{<:Real,1}, Mv::OffsetArray{<:Real,1},
    Mw::OffsetArray{<:Real,1}, alpha::Int, beta::Int, delta::Int)`

"""
moments_conserve(
    Mu::T,
    alpha::I,
) where {T<:OffsetArray{<:AbstractFloat,1},I<:Int} = Mu[alpha]

function moments_conserve(
    Mu::T,
    Mxi::T,
    alpha::I,
    delta::I,
) where {T<:OffsetArray{<:AbstractFloat,1},I<:Int}

    uv = similar(Mu, 3)
    uv[1] = Mu[alpha] * Mxi[delta÷2]
    uv[2] = Mu[alpha+1] * Mxi[delta÷2]
    uv[3] = 0.5 * (Mu[alpha+2] * Mxi[delta÷2] + Mu[alpha] * Mxi[(delta+2)÷2])

    return uv

end

function moments_conserve(
    Mu::T,
    Mv::T,
    Mw::T,
    alpha::I,
    beta::I,
    delta::I,
) where {T<:OffsetArray{<:AbstractFloat,1},I<:Int}

    if length(Mw) == 3 # internal motion
        uv = similar(Mu, 4)
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
        uv = similar(Mu, 5)
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

# ------------------------------------------------------------
# Discrete moments of conservative variables
# ------------------------------------------------------------
#--- 1F1V ---#
function moments_conserve(
    f::X,
    u::T,
    ω::T,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    T<:AbstractArray{<:AbstractFloat,1},
}
    w = similar(f, 3)
    w[1] = discrete_moments(f, u, ω, 0)
    w[2] = discrete_moments(f, u, ω, 1)
    w[3] = 0.5 * discrete_moments(f, u, ω, 2)

    return w
end

#--- 2F1V ---#
function moments_conserve(
    h::X,
    b::X,
    u::T,
    ω::T,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    T<:AbstractArray{<:AbstractFloat,1},
}
    w = similar(h, 3)
    w[1] = discrete_moments(h, u, ω, 0)
    w[2] = discrete_moments(h, u, ω, 1)
    w[3] = 0.5 * (discrete_moments(h, u, ω, 2) + discrete_moments(b, u, ω, 0))

    return w
end

#--- 1F2V ---#
function moments_conserve(
    f::X,
    u::T,
    v::T,
    ω::T,
) where {
    X<:AbstractArray{<:AbstractFloat,2},
    T<:AbstractArray{<:AbstractFloat,2},
}
    w = similar(f, 4)
    w[1] = discrete_moments(f, u, ω, 0)
    w[2] = discrete_moments(f, u, ω, 1)
    w[3] = discrete_moments(f, v, ω, 1)
    w[4] = 0.5 * (discrete_moments(f, u, ω, 2) + discrete_moments(f, v, ω, 2))

    return w
end

#--- 2F2V ---#
function moments_conserve(
    h::X,
    b::X,
    u::T,
    v::T,
    ω::T,
) where {
    X<:AbstractArray{<:AbstractFloat,2},
    T<:AbstractArray{<:AbstractFloat,2},
}
    w = similar(h, 4)
    w[1] = discrete_moments(h, u, ω, 0)
    w[2] = discrete_moments(h, u, ω, 1)
    w[3] = discrete_moments(h, v, ω, 1)
    w[4] =
        0.5 * (
            discrete_moments(h, u, ω, 2) +
            discrete_moments(h, v, ω, 2) +
            discrete_moments(b, u, ω, 0)
        )

    return w
end

#--- 3F2V ---#
function moments_conserve(
    h0::X,
    h1::X,
    h2::X,
    u::T,
    v::T,
    ω::T,
) where {
    X<:AbstractArray{<:AbstractFloat,2},
    T<:AbstractArray{<:AbstractFloat,2},
}
    w = similar(h0, 5)
    w[1] = discrete_moments(h0, u, ω, 0)
    w[2] = discrete_moments(h0, u, ω, 1)
    w[3] = discrete_moments(h0, v, ω, 1)
    w[4] = discrete_moments(h1, u, ω, 0)
    w[5] =
        0.5 * (
            discrete_moments(h0, u, ω, 2) +
            discrete_moments(h0, v, ω, 2) +
            discrete_moments(h2, u, ω, 0)
        )

    return w
end

#--- 1F3V ---#
function moments_conserve(
    f::X,
    u::T,
    v::T,
    w::T,
    ω::T,
) where {
    X<:AbstractArray{<:AbstractFloat,3},
    T<:AbstractArray{<:AbstractFloat,3},
}
    moments = similar(f, 5)

    moments[1] = discrete_moments(f, u, ω, 0)
    moments[2] = discrete_moments(f, u, ω, 1)
    moments[3] = discrete_moments(f, v, ω, 1)
    moments[4] = discrete_moments(f, w, ω, 1)
    moments[5] =
        0.5 * (
            discrete_moments(f, u, ω, 2) +
            discrete_moments(f, v, ω, 2) +
            discrete_moments(f, w, ω, 2)
        )

    return moments
end

#--- 4F1V ---#
function moments_conserve(
    h0::X,
    h1::X,
    h2::X,
    h3::X,
    u::T,
    ω::T,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    T<:AbstractArray{<:AbstractFloat,1},
}
    moments = similar(h0, 5)

    moments[1] = discrete_moments(h0, u, ω, 0)
    moments[2] = discrete_moments(h0, u, ω, 1)
    moments[3] = discrete_moments(h1, u, ω, 0)
    moments[4] = discrete_moments(h2, u, ω, 0)
    moments[5] =
        0.5 * discrete_moments(h0, u, ω, 2) +
        0.5 * discrete_moments(h3, u, ω, 0)

    return moments
end


"""
Calculate slope-related conservative moments
`a = a1 + u * a2 + 0.5 * u^2 * a3`

"""
moments_conserve_slope(
    a,
    Mu::T,
    alpha::I,
) where {T<:OffsetArray{<:Real,1},I<:Int} = a * moments_conserve(Mu, alpha)

moments_conserve_slope(
    a::X,
    Mu::Y,
    Mxi::Y,
    alpha::I,
) where {X<:AbstractArray{<:Real,1},Y<:OffsetArray{<:Real,1},I<:Int} =
    a[1] .* moments_conserve(Mu, Mxi, alpha + 0, 0) .+
    a[2] .* moments_conserve(Mu, Mxi, alpha + 1, 0) .+
    0.5 * a[3] .* moments_conserve(Mu, Mxi, alpha + 2, 0) .+
    0.5 * a[3] .* moments_conserve(Mu, Mxi, alpha + 0, 2)

function moments_conserve_slope(
    a::X,
    Mu::Y,
    Mv::Y,
    Mxi::Y,
    alpha::I,
    beta::I,
) where {X<:AbstractArray{<:Real,1},Y<:OffsetArray{<:Real,1},I<:Int}

    return a[1] .* moments_conserve(Mu, Mv, Mxi, alpha + 0, beta + 0, 0) .+
           a[2] .* moments_conserve(Mu, Mv, Mxi, alpha + 1, beta + 0, 0) .+
           a[3] .* moments_conserve(Mu, Mv, Mxi, alpha + 0, beta + 1, 0) .+
           0.5 * a[4] .*
           moments_conserve(Mu, Mv, Mxi, alpha + 2, beta + 0, 0) .+
           0.5 * a[4] .*
           moments_conserve(Mu, Mv, Mxi, alpha + 0, beta + 2, 0) .+
           0.5 * a[4] .* moments_conserve(Mu, Mv, Mxi, alpha + 0, beta + 0, 2)

end

function moments_conserve_slope(
    a::X,
    Mu::Y,
    Mv::Y,
    Mw::Y,
    alpha::I,
    beta::I,
    delta::I,
) where {X<:AbstractArray{<:Real,1},Y<:OffsetArray{<:Real,1},I<:Int}

    return a[1] .*
           moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 0, delta + 0) .+
           a[2] .*
           moments_conserve(Mu, Mv, Mw, alpha + 1, beta + 0, delta + 0) .+
           a[3] .*
           moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 1, delta + 0) .+
           a[4] .*
           moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 0, delta + 1) .+
           0.5 * a[5] .*
           moments_conserve(Mu, Mv, Mw, alpha + 2, beta + 0, delta + 0) .+
           0.5 * a[5] .*
           moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 2, delta + 0) .+
           0.5 * a[5] .*
           moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 0, delta + 2)

end


"""
Discrete moments of particle distribution

* `discrete_moments(f, ω)`: direct quadrature
* `discrete_moments(f, u, ω, n)`: velocity moments

"""
discrete_moments(
    f::X,
    ω::T,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    T<:AbstractArray{<:AbstractFloat,1},
} = sum(@. ω * f)

#--- 1V ---#
discrete_moments(
    f::X,
    u::T,
    ω::T,
    n,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    T<:AbstractArray{<:AbstractFloat,1},
} = sum(@. ω * u^n * f)

#--- 2V ---#
discrete_moments(
    f::X,
    u::T,
    ω::T,
    n,
) where {
    X<:AbstractArray{<:AbstractFloat,2},
    T<:AbstractArray{<:AbstractFloat,2},
} = sum(@. ω * u^n * f)

#--- 3V ---#
discrete_moments(
    f::X,
    u::T,
    ω::T,
    n,
) where {
    X<:AbstractArray{<:AbstractFloat,3},
    T<:AbstractArray{<:AbstractFloat,3},
} = sum(@. ω * u^n * f)


"""
Calculate stress tensor from particle distribution function

"""
stress(
    f::X,
    prim::Y,
    u::Z,
    ω::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:Real,1},
    Z<:AbstractArray{<:AbstractFloat,1},
} = sum(@. ω * (u - prim[2]) * (u - prim[2]) * f)

function stress(
    f::X,
    prim::Y,
    u::Z,
    v::Z,
    ω::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,2},
    Y<:AbstractArray{<:Real,1},
    Z<:AbstractArray{<:AbstractFloat,2},
}

    P = similar(prim, 2, 2)

    P[1, 1] = sum(@. ω * (u - prim[2]) * (u - prim[2]) * f)
    P[1, 2] = sum(@. ω * (u - prim[2]) * (v - prim[3]) * f)
    P[2, 1] = P[1, 2]
    P[1, 2] = sum(@. ω * (v - prim[3]) * (v - prim[3]) * f)

    return P

end


"""
Calculate heat flux from particle distribution function

"""
heat_flux(
    h::X,
    prim::Y,
    u::Z,
    ω::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:Real,1},
    Z<:AbstractArray{<:AbstractFloat,1},
} = 0.5 * sum(@. ω * (u - prim[2]) * (u - prim[2])^2 * h) # 1F1V

#--- 2F1V ---#
heat_flux(
    h::X,
    b::X,
    prim::Y,
    u::Z,
    ω::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:Real,1},
    Z<:AbstractArray{<:AbstractFloat,1},
} =
    0.5 * (
        sum(@. ω * (u - prim[2]) * (u - prim[2])^2 * h) +
        sum(@. ω * (u - prim[2]) * b)
    )

#--- 1F2V ---#
function heat_flux(
    h::X,
    prim::Y,
    u::Z,
    v::Z,
    ω::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,2},
    Y<:AbstractArray{<:Real,1},
    Z<:AbstractArray{<:AbstractFloat,2},
}

    q = similar(h, 2)
    q[1] =
        0.5 *
        sum(@. ω * (u - prim[2]) * ((u - prim[2])^2 + (v - prim[3])^2) * h)
    q[2] =
        0.5 *
        sum(@. ω * (v - prim[3]) * ((u - prim[2])^2 + (v - prim[3])^2) * h)

    return q

end

#--- 2F2V ---#
function heat_flux(
    h::X,
    b::X,
    prim::Y,
    u::Z,
    v::Z,
    ω::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,2},
    Y<:AbstractArray{<:Real,1},
    Z<:AbstractArray{<:AbstractFloat,2},
}

    q = similar(h, 2)

    q[1] =
        0.5 * (
            sum(@. ω *
                   (u - prim[2]) *
                   ((u - prim[2])^2 + (v - prim[3])^2) *
                   h) + sum(@. ω * (u - prim[2]) * b)
        )
    q[2] =
        0.5 * (
            sum(@. ω *
                   (v - prim[3]) *
                   ((u - prim[2])^2 + (v - prim[3])^2) *
                   h) + sum(@. ω * (v - prim[3]) * b)
        )

    return q

end

#--- 1F3V ---#
function heat_flux(
    f::X,
    prim::Y,
    u::Z,
    v::Z,
    w::Z,
    ω::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,3},
    Y<:AbstractArray{<:Real,1},
    Z<:AbstractArray{<:AbstractFloat,3},
}

    q = similar(f, 3)

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
