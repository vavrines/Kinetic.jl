# ============================================================
# Thermodynamics
# ============================================================

"""
Calculate heat capacity ratio

"""
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


"""
Calculate speed of sound

"""
sound_speed(λ::Real, γ::Real) = (0.5 * γ / λ)^0.5

sound_speed(prim::AbstractArray{<:Real,1}, γ::Real) = sound_speed(prim[end], γ)


function sound_speed(prim::AbstractArray{<:Real,2}, γ::Real)
    c = similar(prim, axes(prim, 2))
    for j in eachindex(c)
        c[j] = sound_speed(prim[end, j], γ)
    end

    return maximum(c)
end