# ------------------------------------------------------------
# Reconstruction and Slope Limiters
# ------------------------------------------------------------


export vanleer, minmod, superbee, vanalbaba
export reconstruct2, reconstruct3


vanleer(sL::Real, sR::Real) =
    (fortsign(1., sL) + fortsign(1., sR)) * abs(sL) * abs(sR) /
    (abs(sL) + abs(sR) + 1.e-7)


minmod(sL::Real, sR::Real) =
    0.5 * (fortsign(1., sL) + fortsign(1., sR)) * min(abs(sR), abs(sL))


function superbee(sL::Real, sR::Real)

    if sR >= 0.5 * sL && sR <= 2. * sL
        return 0.5 * (fortsign(1., sL) + fortsign(1., sR)) * max(abs(sL), abs(sR))
    elseif sR < 0.5 * sL && sR > 2. * sL
        return (fortsign(1., sL) + fortsign(1., sR)) * min(abs(sL), abs(sR))
    else
        return 0.
    end

end


vanalbaba(sL, sR) = (sL^2 * sR + sL * sR^2) / (sL^2 + sR^2 + 1.e-7)


reconstruct2(wL::Real, wR::Real, Δx::Real) = (wR - wL) / Δx

reconstruct2(wL::AbstractArray{<:Real,1}, wR::AbstractArray{<:Real,1}, Δx::Real) =
    (wR .- wL) ./ Δx

function reconstruct2(
    wL::AbstractArray{<:Real,2},
    wR::AbstractArray{<:Real,2},
    Δx::Real,
)
    s = zeros(axes(wL))
    for j in axes(s, 2)
        s[:, j] .= reconstruct2(wL[:, j], wR[:, j], Δx)
    end

    return s
end

function reconstruct2(
    wL::AbstractArray{<:Real,3},
    wR::AbstractArray{<:Real,3},
    Δx::Real,
)
    s = zeros(axes(wL))
    for k in axes(s, 3), j in axes(s, 2)
        s[:, j, k] .= reconstruct2(wL[:, j, k], wR[:, j, k], Δx)
    end

    return s
end


function reconstruct3(
    wL::Real,
    wN::Real,
    wR::Real,
    ΔxL::Real,
    ΔxR::Real,
    limiter = "vanleer"::AbstractString,
)

    sL = (wN - wL) / ΔxL
    sR = (wR - wN) / ΔxR

    if limiter == "linear"
        s = 0.5 * (sL + sR)
    elseif limiter == "vanleer"
        s = vanleer(sL, sR)
    elseif limiter == "minmod"
        s = minmod(sL, sR)
    else
    end

    return s

end

function reconstruct3(
    wL::AbstractArray{<:Real,1},
    wN::AbstractArray{<:Real,1},
    wR::AbstractArray{<:Real,1},
    ΔxL::Real,
    ΔxR::Real,
    limiter = "vanleer"::AbstractString,
)

    sL = (wN .- wL) ./ ΔxL
    sR = (wR .- wN) ./ ΔxR

    if limiter == "linear"
        s = 0.5 .* (sL .+ sR)
    elseif limiter == "vanleer"
        s = vanleer.(sL, sR)
    elseif limiter == "minmod"
        s = minmod.(sL, sR)
    else
    end

    return s

end

function reconstruct3(
    wL::AbstractArray{<:Real,2},
    wN::AbstractArray{<:Real,2},
    wR::AbstractArray{<:Real,2},
    ΔxL::Real,
    ΔxR::Real,
    limiter = "vanleer"::AbstractString,
)

    s = zeros(axes(wL))

    for j in axes(s, 2)
        s[:, j] .= reconstruct3(wL[:, j], wN[:, j], wR[:, j], ΔxL, ΔxR, limiter)
    end

    return s

end

function reconstruct3(
    wL::AbstractArray{<:Real,3},
    wN::AbstractArray{<:Real,3},
    wR::AbstractArray{<:Real,3},
    ΔxL::Real,
    ΔxR::Real,
    limiter = "vanleer"::AbstractString,
)

    s = zeros(axes(wL))

    for k in axes(s, 3), j in axes(s, 2)
        s[:, j, k] .= reconstruct3(wL[:, j, k], wN[:, j, k], wR[:, j, k], ΔxL, ΔxR, limiter)
    end

    return s

end
