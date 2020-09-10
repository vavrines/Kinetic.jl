# ============================================================
# Reconstruction and Slope Limiters
# ============================================================

export vanleer, minmod, superbee, vanalbaba
export reconstruct2, reconstruct2!, reconstruct3, reconstruct3!
export weno5

# ------------------------------------------------------------
# Slope limiter functions
# ------------------------------------------------------------

vanleer(sL::Real, sR::Real) =
    (fortsign(1.0, sL) + fortsign(1.0, sR)) * abs(sL) * abs(sR) /
    (abs(sL) + abs(sR) + 1.e-7)


minmod(sL::Real, sR::Real) =
    0.5 * (fortsign(1.0, sL) + fortsign(1.0, sR)) * min(abs(sR), abs(sL))


function superbee(sL::Real, sR::Real)

    if sR >= 0.5 * sL && sR <= 2.0 * sL
        return 0.5 * (fortsign(1.0, sL) + fortsign(1.0, sR)) * max(abs(sL), abs(sR))
    elseif sR < 0.5 * sL && sR > 2.0 * sL
        return (fortsign(1.0, sL) + fortsign(1.0, sR)) * min(abs(sL), abs(sR))
    else
        return 0.0
    end

end


vanalbaba(sL, sR) = (sL^2 * sR + sL * sR^2) / (sL^2 + sR^2 + 1.e-7)


# ------------------------------------------------------------
# Reconstruction methodologies
# ------------------------------------------------------------

"""
Two-cell reconstruction

"""
reconstruct2(wL::Real, wR::Real, Δx::Real) = (wR - wL) / Δx


reconstruct2(wL::AbstractArray{<:Real,1}, wR::AbstractArray{<:Real,1}, Δx::Real) =
    (wR .- wL) ./ Δx


function reconstruct2!(
    sw::AbstractArray{<:AbstractFloat,1},
    wL::AbstractArray{<:Real,1},
    wR::AbstractArray{<:Real,1},
    Δx::Real,
)
    sw .= (wR .- wL) ./ Δx
end


function reconstruct2(wL::AbstractArray{<:Real,2}, wR::AbstractArray{<:Real,2}, Δx::Real)

    s = zeros(axes(wL))
    for j in axes(s, 2)
        s[:, j] .= reconstruct2(wL[:, j], wR[:, j], Δx)
    end

    return s

end


function reconstruct2!(
    sw::AbstractArray{<:AbstractFloat,2},
    wL::AbstractArray{<:Real,2},
    wR::AbstractArray{<:Real,2},
    Δx::Real,
)

    for j in axes(sw, 2)
        swj = @view sw[:, j]
        reconstruct2!(swj, wL[:, j], wR[:, j], Δx)
    end

end


function reconstruct2(wL::AbstractArray{<:Real,3}, wR::AbstractArray{<:Real,3}, Δx::Real)

    s = zeros(axes(wL))
    for k in axes(s, 3), j in axes(s, 2)
        s[:, j, k] .= reconstruct2(wL[:, j, k], wR[:, j, k], Δx)
    end

    return s

end


function reconstruct2!(
    sw::AbstractArray{<:AbstractFloat,3},
    wL::AbstractArray{<:Real,3},
    wR::AbstractArray{<:Real,3},
    Δx::Real,
)

    for k in axes(sw, 3), j in axes(sw, 2)
        swjk = @view sw[:, j, k]
        reconstruct2!(swjk, wL[:, j, k], wR[:, j, k], Δx)
    end

end


"""
Three-cell reconstruction

"""
function reconstruct3(
    wL::Real,
    wN::Real,
    wR::Real,
    ΔxL::Real,
    ΔxR::Real,
    limiter = :vanleer::Symbol,
)

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
    wL::AbstractArray{<:Real,1},
    wN::AbstractArray{<:Real,1},
    wR::AbstractArray{<:Real,1},
    ΔxL::Real,
    ΔxR::Real,
    limiter = :vanleer::Symbol,
)

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


function reconstruct3!(
    sw::AbstractArray{<:AbstractFloat,1},
    wL::AbstractArray{<:Real,1},
    wN::AbstractArray{<:Real,1},
    wR::AbstractArray{<:Real,1},
    ΔxL::Real,
    ΔxR::Real,
    limiter = :vanleer::Symbol,
)

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


function reconstruct3(
    wL::AbstractArray{<:Real,2},
    wN::AbstractArray{<:Real,2},
    wR::AbstractArray{<:Real,2},
    ΔxL::Real,
    ΔxR::Real,
    limiter = :vanleer::Symbol,
)

    s = zeros(axes(wL))
    for j in axes(s, 2)
        s[:, j] .= reconstruct3(wL[:, j], wN[:, j], wR[:, j], ΔxL, ΔxR, limiter)
    end

    return s

end


function reconstruct3!(
    sw::AbstractArray{<:AbstractFloat,2},
    wL::AbstractArray{<:Real,2},
    wN::AbstractArray{<:Real,2},
    wR::AbstractArray{<:Real,2},
    ΔxL::Real,
    ΔxR::Real,
    limiter = :vanleer::Symbol,
)

    for j in axes(sw, 2)
        swj = @view sw[:, j]
        reconstruct3!(swj, wL[:, j], wN[:, j], wR[:, j], ΔxL, ΔxR, limiter)
    end

end


function reconstruct3(
    wL::AbstractArray{<:Real,3},
    wN::AbstractArray{<:Real,3},
    wR::AbstractArray{<:Real,3},
    ΔxL::Real,
    ΔxR::Real,
    limiter = :vanleer::Symbol,
)

    s = zeros(axes(wL))
    for k in axes(s, 3), j in axes(s, 2)
        s[:, j, k] .= reconstruct3(wL[:, j, k], wN[:, j, k], wR[:, j, k], ΔxL, ΔxR, limiter)
    end

    return s

end


function reconstruct3!(
    sw::AbstractArray{<:AbstractFloat,3},
    wL::AbstractArray{<:Real,3},
    wN::AbstractArray{<:Real,3},
    wR::AbstractArray{<:Real,3},
    ΔxL::Real,
    ΔxR::Real,
    limiter = :vanleer::Symbol,
)

    for k in axes(sw, 3), j in axes(sw, 2)
        swjk = @view sw[:, j, k]
        reconstruct3!(swjk, wL[:, j, k], wN[:, j, k], wR[:, j, k], ΔxL, ΔxR, limiter)
    end

end

function reconstruct3!(
    sw::AbstractArray{<:AbstractFloat,4},
    wL::AbstractArray{<:Real,4},
    wN::AbstractArray{<:Real,4},
    wR::AbstractArray{<:Real,4},
    ΔxL::Real,
    ΔxR::Real,
    limiter = :vanleer::Symbol,
)

    for l in axes(sw, 4), k in axes(sw, 3), j in axes(sw, 2)
        sjkl = @view sw[:, j, k, l]
        reconstruct3!(sjkl, wL[:, j, k, l], wN[:, j, k, l], wR[:, j, k, l], ΔxL, ΔxR, limiter)
    end

end


function weno5(wL2::Real, wL1::Real, wN::Real, wR1::Real, wR2::Real)

    ϵ = 1e-6

    β0 = 13.0 / 12.0 * (wN - 2.0 * wR1 + wR2)^2 + 1.0 / 4.0 * (3.0 * wN - 4.0 * wR1 + wR2)^2
    β1 = 13.0 / 12.0 * (wL1 - 2.0 * wN + wR1)^2 + 1.0 / 4.0 * (wL1 - wR1)^2
    β2 = 13.0 / 12.0 * (wL2 - 2.0 * wL1 + wN)^2 + 1.0 / 4.0 * (wL2 - 4.0 * wL1 + 3.0 * wN)^2

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
