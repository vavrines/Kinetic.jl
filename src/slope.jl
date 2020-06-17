# ------------------------------------------------------------
# Reconstruction and Slope Limiters
# ------------------------------------------------------------


export vanleer, minmod, reconstruct2, reconstruct3


vanleer(sL, sR) =
    (fortsign(1.0, sL) + fortsign(1.0, sR)) * abs(sL) * abs(sR) /
    (abs(sL) + abs(sR) + 1.e-7)


minmod(sL, sR) = 0.5 * (fortsign(1.0, sL) + fortsign(1.0, sR)) * min(abs(sR), abs(sL))


reconstruct2(wL::Float64, wR::Float64, Δx::Float64) = (wR - wL) / Δx

reconstruct2(wL::AbstractArray{Float64,1}, wR::AbstractArray{Float64,1}, Δx::Float64) =
    (wR .- wL) ./ Δx

function reconstruct2(
    wL::AbstractArray{Float64,2},
    wR::AbstractArray{Float64,2},
    Δx::Float64,
)
    s = zeros(axes(wL))
    for j in axes(s, 2)
        s[:, j] .= reconstruct2(wL[:, j], wR[:, j], Δx)
    end

    return s
end

function reconstruct2(
    wL::AbstractArray{Float64,3},
    wR::AbstractArray{Float64,3},
    Δx::Float64,
)
    s = zeros(axes(wL))
    for k in axes(s, 3), j in axes(s, 2)
        s[:, j, k] .= reconstruct2(wL[:, j, k], wR[:, j, k], Δx)
    end

    return s
end


function reconstruct3(
    wL::Float64,
    wN::Float64,
    wR::Float64,
    ΔxL::Float64,
    ΔxR::Float64,
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
    wL::AbstractArray{Float64,1},
    wN::AbstractArray{Float64,1},
    wR::AbstractArray{Float64,1},
    ΔxL::Float64,
    ΔxR::Float64,
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
    wL::AbstractArray{Float64,2},
    wN::AbstractArray{Float64,2},
    wR::AbstractArray{Float64,2},
    ΔxL::Float64,
    ΔxR::Float64,
    limiter = "vanleer"::AbstractString,
)

    s = zeros(axes(wL))

    for j in axes(s, 2)
        s[:, j] .= reconstruct3(wL[:, j], wN[:, j], wR[:, j], ΔxL, ΔxR, limiter)
    end

    return s

end

function reconstruct3(
    wL::AbstractArray{Float64,3},
    wN::AbstractArray{Float64,3},
    wR::AbstractArray{Float64,3},
    ΔxL::Float64,
    ΔxR::Float64,
    limiter = "vanleer"::AbstractString,
)

    s = zeros(axes(wL))

    for k in axes(s, 3), j in axes(s, 2)
        s[:, j, k] .= reconstruct3(wL[:, j, k], wN[:, j, k], wR[:, j, k], ΔxL, ΔxR, limiter)
    end

    return s

end
