# ============================================================
# Theory
# ============================================================


export velocity_moments, 
	   maxwellian, 
	   conserve_prim, 
	   prim_conserve, 
	   heat_capacity_ratio, 
	   ref_vis, 
	   sos, 
	   collision_time


# ------------------------------------------------------------
# Calculate conservative moments from distribution function
# ------------------------------------------------------------
velocity_moments(f::AbstractArray{Float64,1}, u::AbstractArray{Float64,1}, ω::AbstractArray{Float64,1}, n::Int) =
sum(@. ω * u^n * f)


velocity_moments(f::AbstractArray{Float64,2}, u::AbstractArray{Float64,2}, ω::AbstractArray{Float64,2}, n::Int) =
sum(@. ω * u^n * f)


# ------------------------------------------------------------
# Calculate equilibrium distribution function
# ------------------------------------------------------------
maxwellian(u::AbstractArray{Float64,1}, ρ::Union{Int,Float64}, U::Union{Int,Float64}, λ::Union{Int,Float64}) =
@. ρ * (λ / π)^0.5 * exp(-λ * (u - U)^2)

maxwellian(u::AbstractArray{Float64,1}, prim::Array{Float64,1}) =
maxwellian(u, prim[1], prim[2], prim[3])

maxwellian(u::AbstractArray{Float64,1}, prim::Array{Int,1}) =
maxwellian(u, Float64.(prim))


maxwellian(u::AbstractArray{Float64,2}, v::AbstractArray{Float64,2}, ρ::Union{Int,Float64}, U::Union{Int,Float64}, V::Union{Int,Float64}, λ::Union{Int,Float64}) =
@. ρ * (λ / π) * exp(-λ * ((u - U)^2 + (v - V)^2))

maxwellian(u::AbstractArray{Float64,2}, v::AbstractArray{Float64,2}, prim::Array{Float64,1}) =
maxwellian(u, v, prim[1], prim[2], prim[3], prim[4])

maxwellian(u::AbstractArray{Float64,2}, v::AbstractArray{Float64,2}, prim::Array{Int,1}) =
maxwellian(u, v, Float64.(prim))


# ------------------------------------------------------------
# Calculate conservative/primitive variables
# ------------------------------------------------------------
function prim_conserve(prim::Array{Float64,1}, γ::Union{Int,Float64})

	W = similar(prim)

	if length(prim) == 3 # 1D
		W[1] = prim[1]
		W[2] = prim[1] * prim[2]
		W[3] = 0.5 * prim[1] / prim[3] / (γ - 1.0) + 0.5 * prim[1] * prim[2]^2
	elseif length(prim) == 4 # 2D
		W[1] = prim[1]
		W[2] = prim[1] * prim[2]
		W[3] = prim[1] * prim[3]
		W[4] = 0.5 * prim[1] / prim[4] / (γ - 1.0) + 0.5 * prim[1] * (prim[2]^2 + prim[3]^2)
	elseif length(prim) == 5 # 3D
		W[1] = prim[1]
		W[2] = prim[1] * prim[2]
		W[3] = prim[1] * prim[3]
		W[4] = prim[1] * prim[4]
		W[5] = 0.5 * prim[1] / prim[5] / (γ - 1.0) + 0.5 * prim[1] * (prim[2]^2 + prim[3]^2 + prim[4]^2)
	else
		println("prim -> w : dimension error")
	end

	return W

end

prim_conserve(prim::Array{Int,1}, γ::Union{Int,Float64}) = prim_conserve(Float64.(prim), γ)

prim_conserve(ρ::Union{Int,Float64}, U::Union{Int,Float64}, λ::Union{Int,Float64}, γ::Union{Int,Float64}) = 
prim_conserve([ρ, U, λ], γ)

prim_conserve(ρ::Union{Int,Float64}, U::Union{Int,Float64}, V::Union{Int,Float64}, λ::Union{Int,Float64}, γ::Union{Int,Float64}) = 
prim_conserve([ρ, U, V, λ], γ)


function conserve_prim(W::Array{Float64,1}, γ::Union{Int,Float64})

	prim = similar(W)

	if length(W) == 3 # 1D
		prim[1] = W[1]
		prim[2] = W[2] / W[1]
		prim[3] = 0.5 * W[1] / (γ - 1.0) / (W[3] - 0.5 * W[2]^2 / W[1])
	elseif length(W) == 4 # 2D
		prim[1] = W[1]
		prim[2] = W[2] / W[1]
		prim[3] = W[3] / W[1]
		prim[4] = 0.5 * W[1] / (γ - 1.0) / (W[4] - 0.5 * (W[2]^2 + W[3]^2) / W[1])
	elseif length(W) == 5 # 3D
		prim[1] = W[1]
		prim[2] = W[2] / W[1]
		prim[3] = W[3] / W[1]
		prim[4] = W[4] / W[1]
		prim[5] = 0.5 * W[1] / (γ - 1.0) / (W[5] - 0.5 * (W[2]^2 + W[3]^2 + W[4]^2) / W[1])
	else
		println("w -> prim : dimension error")
	end

	return prim

end

conserve_prim(W::Array{Int,1}, γ::Union{Int,Float64}) = conserve_prim(Float64.(W), γ)

conserve_prim(ρ::Union{Int,Float64}, M::Union{Int,Float64}, E::Union{Int,Float64}, gamma::Union{Int,Float64}) = 
conserve_prim([ρ, M, E], gamma)

conserve_prim(ρ::Union{Int,Float64}, MX::Union{Int,Float64}, MY::Union{Int,Float64}, E::Union{Int,Float64}, gamma::Union{Int,Float64}) = 
conserve_prim([ρ, MX, MY, E], gamma)


# ------------------------------------------------------------
# Calculate heat capacity ratio
# ------------------------------------------------------------
function heat_capacity_ratio(K::Union{Int, Float64}, D::Int)
	
	if D == 1
		γ = (K + 3.) / (K + 1.)
	elseif D == 2
		γ = (K + 4.) / (K + 2.)
	elseif D ==3
		γ = (K + 5.) / (K + 3.)
	end

	return γ

end


# ------------------------------------------------------------
# Calculate speed of sound
# ------------------------------------------------------------
sos(λ::Union{Int,Float64}, γ::Union{Int,Float64}) = (0.5 * γ / λ)^0.5

sos(prim::Array{Float64,1}, γ::Union{Int,Float64}) = sos(prim[end], γ)

sos(prim::Array{Int,1}, γ::Union{Int,Float64}) = sos(Float64.(prim), γ)


# ------------------------------------------------------------
# Calculate reference viscosity
# 1. variable hard sphere (VHS) model
# ------------------------------------------------------------
ref_vis(Kn::Union{Int,Float64}, alpha::Union{Int,Float64}, omega::Union{Int,Float64}) = 
5.0 * (alpha + 1.) * (alpha + 2.) * √π / (4. * alpha * (5. - 2. * omega) * (7. - 2. * omega)) * Kn


# ------------------------------------------------------------
# Calculate collision time
# 1. variable hard sphere (VHS) model
# ------------------------------------------------------------
collision_time(prim::Array{Float64,1}, muRef::Union{Int,Float64}, omega::Union{Int,Float64}) = 
muRef * 2. * prim[end]^(1. - omega) / prim[1]

collision_time(prim::Array{Int,1}, muRef::Union{Int,Float64}, omega::Union{Int,Float64}) = 
collision_time(Float64.(prim), muRef, omega)