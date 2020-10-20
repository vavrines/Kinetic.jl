# ============================================================
# Reconstruction and Slope Limiters
# ============================================================

export vanleer, minmod, superbee, vanalbaba
export reconstruct2, reconstruct2!, reconstruct3, reconstruct3!
export weno5

# ------------------------------------------------------------
# Slope limiter functions
# ------------------------------------------------------------

"""
van Leer limiter

    vanleer(sL::Real, sR::Real)

"""
vanleer(sL::T, sR::T) where {T} =
    (fortsign(1.0, sL) + fortsign(1.0, sR)) * abs(sL) * abs(sR) /
    (abs(sL) + abs(sR) + 1.e-7)


"""
minmod limiter

    minmod(sL::Real, sR::Real)

"""
minmod(sL::T, sR::T) where {T} =
    0.5 * (fortsign(1.0, sL) + fortsign(1.0, sR)) * min(abs(sR), abs(sL))


"""
superbee limiter

    superbee(sL::Real, sR::Real)

"""
function superbee(sL::T, sR::T) where {T}

    if sR >= 0.5 * sL && sR <= 2.0 * sL
        return 0.5 *
               (fortsign(1.0, sL) + fortsign(1.0, sR)) *
               max(abs(sL), abs(sR))
    elseif sR < 0.5 * sL && sR > 2.0 * sL
        return (fortsign(1.0, sL) + fortsign(1.0, sR)) * min(abs(sL), abs(sR))
    else
        return 0.0
    end

end


"""
van Albaba limiter

    vanalbaba(sL::Real, sR::Real)

"""
vanalbaba(sL::T, sR::T) where {T} =
    (sL^2 * sR + sL * sR^2) / (sL^2 + sR^2 + 1.e-7)

# ------------------------------------------------------------
# Reconstruction methodologies
# ------------------------------------------------------------

"""
Two-cell reconstruction

"""
reconstruct2(wL::X, wR::X, Δx::Y) where {X,Y} = (wR - wL) / Δx


reconstruct2(wL::T, wR::T, Δx) where {T<:AbstractArray{<:Real,1}} =
    (wR .- wL) ./ Δx

function reconstruct2(wL::T, wR::T, Δx) where {T<:AbstractArray{<:Real,2}}

    s = zeros(axes(wL))
    for j in axes(s, 2)
        s[:, j] .= reconstruct2(wL[:, j], wR[:, j], Δx)
    end

    return s

end

function reconstruct2(wL::T, wR::T, Δx) where {T<:AbstractArray{<:Real,3}}

    s = zeros(axes(wL))
    for k in axes(s, 3), j in axes(s, 2)
        s[:, j, k] .= reconstruct2(wL[:, j, k], wR[:, j, k], Δx)
    end

    return s

end


"""
Two-cell reconstruction

"""
function reconstruct2!(
    sw::X,
    wL::Y,
    wR::Y,
    Δx,
) where {X<:AbstractArray{<:AbstractFloat,1},Y<:AbstractArray{<:Real,1}}
    sw .= (wR .- wL) ./ Δx
end

function reconstruct2!(
    sw::X,
    wL::Y,
    wR::Y,
    Δx,
) where {X<:AbstractArray{<:AbstractFloat,2},Y<:AbstractArray{<:Real,2}}

    for j in axes(sw, 2)
        swj = @view sw[:, j]
        reconstruct2!(swj, wL[:, j], wR[:, j], Δx)
    end

end

function reconstruct2!(
    sw::X,
    wL::Y,
    wR::Y,
    Δx,
) where {X<:AbstractArray{<:AbstractFloat,3},Y<:AbstractArray{<:Real,3}}

    for k in axes(sw, 3), j in axes(sw, 2)
        swjk = @view sw[:, j, k]
        reconstruct2!(swjk, wL[:, j, k], wR[:, j, k], Δx)
    end

end


"""
Three-cell reconstruction

"""
function reconstruct3(
    wL::T,
    wN::T,
    wR::T,
    ΔxL::T,
    ΔxR::T,
    limiter = :vanleer::Symbol,
) where {T}

    sL = (wN - wL) / ΔxL
    sR = (wR - wN) / ΔxR

    if limiter == :linear
        return 0.5 * (sL + sR)
    elseif limiter == :vanleer
        return vanleer(sL, sR)
    elseif limiter == :minmod
        return minmod(sL, sR)
    else
    end

end

function reconstruct3(
    wL::T,
    wN::T,
    wR::T,
    ΔxL,
    ΔxR,
    limiter = :vanleer::Symbol,
) where {T<:AbstractArray{<:Real,1}}

    sL = (wN .- wL) ./ ΔxL
    sR = (wR .- wN) ./ ΔxR

    if limiter == :linear
        return 0.5 .* (sL .+ sR)
    elseif limiter == :vanleer
        return vanleer.(sL, sR)
    elseif limiter == :minmod
        return minmod.(sL, sR)
    else
    end

end

function reconstruct3(
    wL::T,
    wN::T,
    wR::T,
    ΔxL,
    ΔxR,
    limiter = :vanleer::Symbol,
) where {T<:AbstractArray{<:Real,2}}

    s = zeros(axes(wL))
    for j in axes(s, 2)
        s[:, j] .= reconstruct3(wL[:, j], wN[:, j], wR[:, j], ΔxL, ΔxR, limiter)
    end

    return s

end

function reconstruct3(
    wL::T,
    wN::T,
    wR::T,
    ΔxL,
    ΔxR,
    limiter = :vanleer::Symbol,
) where {T<:AbstractArray{<:Real,3}}

    s = zeros(axes(wL))
    for k in axes(s, 3), j in axes(s, 2)
        s[:, j, k] .= reconstruct3(
            wL[:, j, k],
            wN[:, j, k],
            wR[:, j, k],
            ΔxL,
            ΔxR,
            limiter,
        )
    end

    return s

end


"""
Three-cell reconstruction

"""
function reconstruct3!(
    sw::X,
    wL::Y,
    wN::Y,
    wR::Y,
    ΔxL,
    ΔxR,
    limiter = :vanleer::Symbol,
) where {X<:AbstractArray{<:AbstractFloat,1},Y<:AbstractArray{<:Real,1}}

    sL = (wN .- wL) ./ ΔxL
    sR = (wR .- wN) ./ ΔxR

    if limiter == :linear
        sw .= 0.5 .* (sL .+ sR)
    elseif limiter == :vanleer
        sw .= vanleer.(sL, sR)
    elseif limiter == :minmod
        sw .= minmod.(sL, sR)
    else
    end

end

function reconstruct3!(
    sw::X,
    wL::Y,
    wN::Y,
    wR::Y,
    ΔxL,
    ΔxR,
    limiter = :vanleer::Symbol,
) where {X<:AbstractArray{<:AbstractFloat,2},Y<:AbstractArray{<:Real,2}}

    for j in axes(sw, 2)
        swj = @view sw[:, j]
        reconstruct3!(swj, wL[:, j], wN[:, j], wR[:, j], ΔxL, ΔxR, limiter)
    end

end

function reconstruct3!(
    sw::X,
    wL::Y,
    wN::Y,
    wR::Y,
    ΔxL,
    ΔxR,
    limiter = :vanleer::Symbol,
) where {X<:AbstractArray{<:AbstractFloat,3},Y<:AbstractArray{<:Real,3}}

    for k in axes(sw, 3), j in axes(sw, 2)
        swjk = @view sw[:, j, k]
        reconstruct3!(
            swjk,
            wL[:, j, k],
            wN[:, j, k],
            wR[:, j, k],
            ΔxL,
            ΔxR,
            limiter,
        )
    end

end

function reconstruct3!(
    sw::X,
    wL::Y,
    wN::Y,
    wR::Y,
    ΔxL,
    ΔxR,
    limiter = :vanleer::Symbol,
) where {X<:AbstractArray{<:AbstractFloat,4},Y<:AbstractArray{<:Real,4}}

    for l in axes(sw, 4), k in axes(sw, 3), j in axes(sw, 2)
        sjkl = @view sw[:, j, k, l]
        reconstruct3!(
            sjkl,
            wL[:, j, k, l],
            wN[:, j, k, l],
            wR[:, j, k, l],
            ΔxL,
            ΔxR,
            limiter,
        )
    end

end


"""
5th-order WENO-JS interpolation

    weno5(wL2::T, wL1::T, wN::T, wR1::T, wR2::T) where {T}

"""
function weno5(wL2::T, wL1::T, wN::T, wR1::T, wR2::T) where {T}

    ϵ = 1e-6

    β0 =
        13.0 / 12.0 * (wN - 2.0 * wR1 + wR2)^2 +
        1.0 / 4.0 * (3.0 * wN - 4.0 * wR1 + wR2)^2
    β1 = 13.0 / 12.0 * (wL1 - 2.0 * wN + wR1)^2 + 1.0 / 4.0 * (wL1 - wR1)^2
    β2 =
        13.0 / 12.0 * (wL2 - 2.0 * wL1 + wN)^2 +
        1.0 / 4.0 * (wL2 - 4.0 * wL1 + 3.0 * wN)^2

    #--- right interface ---#
    dr0 = 0.3
    dr1 = 0.6
    dr2 = 0.1

    αr0 = dr0 / (ϵ + β0)^2
    αr1 = dr1 / (ϵ + β1)^2
    αr2 = dr2 / (ϵ + β2)^2

    ωr0 = αr0 / (αr0 + αr1 + αr2)
    ωr1 = αr1 / (αr0 + αr1 + αr2)
    ωr2 = αr2 / (αr0 + αr1 + αr2)

    qr0 = 1.0 / 3.0 * wN + 5.0 / 6.0 * wR1 - 1.0 / 6.0 * wR2
    qr1 = -1.0 / 6.0 * wL1 + 5.0 / 6.0 * wN + 1.0 / 3.0 * wR1
    qr2 = 1.0 / 3.0 * wL2 - 7.0 / 6.0 * wL1 + 11.0 / 6.0 * wN

    wR = ωr0 * qr0 + ωr1 * qr1 + ωr2 * qr2

    #--- left interface ---#
    dl0 = 0.1
    dl1 = 0.6
    dl2 = 0.3

    αl0 = dl0 / (ϵ + β0)^2
    αl1 = dl1 / (ϵ + β1)^2
    αl2 = dl2 / (ϵ + β2)^2

    ωl0 = αl0 / (αl0 + αl1 + αl2)
    ωl1 = αl1 / (αl0 + αl1 + αl2)
    ωl2 = αl2 / (αl0 + αl1 + αl2)

    ql0 = 11.0 / 6.0 * wN - 7.0 / 6.0 * wR1 + 1.0 / 3.0 * wR2
    ql1 = 1.0 / 3.0 * wL1 + 5.0 / 6.0 * wN - 1.0 / 6.0 * wR1
    ql2 = -1.0 / 6.0 * wL2 + 5.0 / 6.0 * wL1 + 1.0 / 3.0 * wN

    wL = αl0 * ql0 + αl1 * ql1 + αl2 * ql2

    return wL, wR

end
