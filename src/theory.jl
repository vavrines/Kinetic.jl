# ============================================================
# Theory
# ============================================================


export velocity_moments, 
	   maxwellian, 
	   conserve_prim, 
	   prim_conserve, 
	   gas_constant, 
	   ref_vis, 
	   sos, 
	   collision_time


# ------------------------------------------------------------
# Calculate conservative moments from distribution function
# ------------------------------------------------------------
velocity_moments(f::Array{Float64,1}, u::Array{Float64,1}, ω::Array{Float64,1}, N::Int64) =
sum(@. ω * u^N * f)


# ------------------------------------------------------------
# Calculate equilibrium distribution function
# ------------------------------------------------------------
maxwellian(u::Array{Float64,1}, ρ::Float64, U::Float64, λ::Float64) =
@. ρ * (λ / π)^0.5 * exp(-λ * (u - U)^2)


maxwellian(u::Array{Float64,2}, v::Array{Float64,2}, ρ::Float64, U::Float64, V::Float64, λ::Float64) =
@. ρ * (λ / π) * exp(-λ * ((u - U)^2 + (v - V)^2))


# ------------------------------------------------------------
# Calculate conservative/primitive variables
# ------------------------------------------------------------
function prim_conserve(prim::Array{Float64,1}, gamma::Float64)

	W = similar(prim)

	if length(prim) == 3 # 1D
		W[1] = prim[1]
		W[2] = prim[1] * prim[2]
		W[3] = 0.5 * prim[1] / prim[3] / (gamma - 1.0) + 0.5 * prim[1] * prim[2]^2
	elseif length(prim) == 4 # 2D
		W[1] = prim[1]
		W[2] = prim[1] * prim[2]
		W[3] = prim[1] * prim[3]
		W[4] = 0.5 * prim[1] / prim[4] / (gamma - 1.0) + 0.5 * prim[1] * (prim[2]^2 + prim[3]^2)
	else
		println("prim -> w : dimension error")
	end

	return W

end


function conserve_prim(W::Array{Float64,1}, gamma::Float64)

	prim = similar(W)

	if length(W) == 3 # 1D
		prim[1] = W[1]
		prim[2] = W[2] / W[1]
		prim[3] = 0.5 * W[1] / (gamma - 1.0) / (W[3] - 0.5 * W[2]^2 / W[1])
	elseif length(W) == 4 # 2D
		prim[1] = W[1]
		prim[2] = W[2] / W[1]
		prim[3] = W[3] / W[1]
		prim[4] = 0.5 * W[1] / (gamma - 1.0) / (W[4] - 0.5 * (W[2]^2 + W[3]^2) / W[1])
	else
		println("w -> prim : dimension error")
	end

	return prim

end


prim_conserve(ρ::Float64, U::Float64, λ::Float64, gamma::Float64) = 
prim_conserve([ρ, U, λ], gamma)


# ------------------------------------------------------------
# Calculate physical property parameters
# ------------------------------------------------------------
function gas_constant(inK, dim::Int64) 
	
	if dim == 1
		gam = (inK + 3.) / (inK + 1.)
	elseif dim == 2
		gam = (inK + 4.) / (inK + 2.)
	end

	return gam

end


ref_vis(Kn::Float64, alpha::Float64, omega::Float64) = 
5.0 * (alpha + 1.) * (alpha + 2.) * √π / (4. * alpha * (5. - 2. * omega) * (7. - 2. * omega)) * Kn


sos(λ::Float64, γ::Float64) = (0.5 * γ / λ)^0.5

sos(prim::Array{Float64,1}, γ::Float64) = sos(prim[end], γ)


collision_time(prim::Array{Float64,1}, muRef::Float64, omega::Float64) = muRef * 2. * prim[end]^(1. - omega) / prim[1]
