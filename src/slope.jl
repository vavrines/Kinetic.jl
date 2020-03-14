# ------------------------------------------------------------
# Reconstruction and Slope Limiters
# ------------------------------------------------------------


export vanleer,
       minmod,
       reconstruct2,
       reconstruct3


vanleer(sL, sR) = 
(fortsign(1., sL) + fortsign(1., sR)) * abs(sL) * abs(sR) / (abs(sL) + abs(sR) + 1.e-7)


minmod(sL, sR) = 
0.5 * (fortsign(1., sL) + fortsign(1., sR)) * min(abs(sR), abs(sL))


reconstruct2(wL::Float64, wR::Float64, Δx::Float64) =
(wR - wL) / Δx

reconstruct2(wL::AbstractArray{Float64,1}, wR::AbstractArray{Float64,1}, Δx::Float64) =
(wR .- wL) ./ Δx


function reconstruct3( wL::Float64, wN::Float64, wR::Float64, ΔxL::Float64, ΔxR::Float64, 
                       limiter="vanleer"::AbstractString )

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

function reconstruct3( wL::AbstractArray{Float64,1}, wN::AbstractArray{Float64,1}, wR::AbstractArray{Float64,1}, 
                       ΔxL::Float64, ΔxR::Float64, limiter="vanleer"::AbstractString )

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