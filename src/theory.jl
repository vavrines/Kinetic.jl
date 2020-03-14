# ============================================================
# Theory
# ============================================================


export gauss_moments,
	   moments_conserve,
	   moments_conserve_slope,
	   discrete_moments, 
	   maxwellian, 
	   conserve_prim, 
	   prim_conserve, 
	   heat_capacity_ratio, 
	   ref_vhs_vis, 
	   sound_speed, 
	   vhs_collision_time,
	   aap_hs_collision_time,
	   aap_hs_prim


"""
Velocity moments of particle distribution function
1. theoretical form
"""

# ------------------------------------------------------------
# Calculate directional velocity moments of Gaussian
# G = (λ / π)^(D / 2) * exp[-λ(c^2 + ξ^2)]
# ------------------------------------------------------------
function gauss_moments(prim::Array{<:Real,1}, inK::Real)

	MuL = OffsetArray{Float64}(undef, 0:6); MuR = similar(MuL); Mu = similar(MuL)

	MuL[0] = 0.5 * SpecialFunctions.erfc(-sqrt(prim[end]) * prim[2])
	MuL[1] = prim[2] * MuL[0] + 0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])
	MuR[0] = 0.5 * SpecialFunctions.erfc(sqrt(prim[end]) * prim[2])
	MuR[1] = prim[2] * MuR[0] - 0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])

	Threads.@threads for i=2:6
		MuL[i] = prim[2] * MuL[i-1] + 0.5 * (i-1) * MuL[i-2] / prim[end]
		MuR[i] = prim[2] * MuR[i-1] + 0.5 * (i-1) * MuR[i-2] / prim[end]
	end

	@. Mu = MuL + MuR

	if length(prim) == 3
		Mxi = OffsetArray{Float64}(undef, 0:2)
		Mxi[0] = 1.0
		Mxi[1] = 0.5 * inK / prim[end]
		Mxi[2] = (inK^2 + 2.0 * inK) / (4.0 * prim[end]^2)

		return Mu, Mxi, MuL, MuR
	elseif length(prim) == 4
		Mv = OffsetArray{Float64}(undef, 0:6)
		Mv[0] = 1.0
    	Mv[1] = prim[3]
		Threads.@threads for i=2:6
			Mv[i] = prim[3] * Mv[i-1] + 0.5 * (i-1) * Mv[i-2] / prim[end]
		end

		Mxi = OffsetArray{Float64}(undef, 0:2)
		Mxi[0] = 1.0
		Mxi[1] = 0.5 * inK / prim[end]
		Mxi[2] = (inK^2 + 2.0 * inK) / (4.0 * prim[end]^2)

		return Mu, Mv, Mxi, MuL, MuR
	elseif length(prim) == 5
		Mv = OffsetArray{Float64}(undef, 0:6)
		Mv[0] = 1.
		Mv[1] = prim[3]
		Threads.@threads for i=2:6
			Mv[i] = prim[3] * Mv[i-1] + 0.5 * (i-1) * Mv[i-2] / prim[end]
		end

		Mw = OffsetArray{Float64}(undef, 0:6)
		Mw[0] = 1.
		Mw[1] = prim[4]
		Threads.@threads for i=2:6
			Mw[i] = prim[4] * Mw[i-1] + 0.5 * (i-1) * Mw[i-2] / prim[end]
		end

		return Mu, Mv, Mw, MuL, MuR
	end

end


# ------------------------------------------------------------
# Calculate conservative moments
# ------------------------------------------------------------
function moments_conserve(Mu::OffsetArray{Float64,1}, Mxi::OffsetArray{Float64,1}, alpha::Int, delta::Int)

    uv = zeros(3)

    uv[1] = Mu[alpha] * Mxi[delta÷2]
    uv[2] = Mu[alpha+1] * Mxi[delta÷2]
    uv[3] = 0.5 * (Mu[alpha+2] * Mxi[delta÷2] + Mu[alpha] * Mxi[(delta+2)÷2])

    return uv

end


function moments_conserve( Mu::OffsetArray{Float64,1}, Mv::OffsetArray{Float64,1}, Mw::OffsetArray{Float64,1}, 
						   alpha::Int, beta::Int, delta::Int )

	if length(Mw) == 3
		uv = zeros(4)

		uv[1] = Mu[alpha] * Mv[beta] * Mw[delta÷2]
		uv[2] = Mu[alpha+1] * Mv[beta] * Mw[delta÷2]
		uv[3] = Mu[alpha] * Mv[beta+1] * Mw[delta÷2]
		uv[4] = 0.5 * (Mu[alpha+2] * Mv[beta] * Mw[delta÷2] + Mu[alpha] * Mv[beta+2] * Mw[delta÷2] + 
				Mu[alpha] * Mv[beta] * Mw[(delta+2)÷2])	
	else
		uv = zeros(5)

		uv[1] = Mu[alpha] * Mv[beta] * Mw[delta]
		uv[2] = Mu[alpha+1] * Mv[beta] * Mw[delta]
		uv[3] = Mu[alpha] * Mv[beta+1] * Mw[delta]
		uv[4] = Mu[alpha] * Mv[beta] * Mw[delta+1]
		uv[5] = 0.5 * (Mu[alpha+2] * Mv[beta] * Mw[delta] + Mu[alpha] * Mv[beta+2] * Mw[delta] + 
				Mu[alpha] * Mv[beta] * Mw[delta+2])
	end

	return uv

end


# ------------------------------------------------------------
# Calculate slope-related conservative moments
# a = a1 + u * a2 + 0.5 * u^2 * a3
# ------------------------------------------------------------
function moments_conserve_slope(a::Array{Float64,1}, Mu::OffsetArray{Float64,1}, Mxi::OffsetArray{Float64,1}, alpha::Int)

	au = @. a[1] * moments_conserve(Mu, Mxi, alpha+0, 0) +
         a[2] * moments_conserve(Mu, Mxi, alpha+1, 0) +
         0.5 * a[3] * moments_conserve(Mu, Mxi, alpha+2, 0) +
         0.5 * a[3] * moments_conserve(Mu, Mxi, alpha+0, 2)

    return au

end


function moments_conserve_slope( a::Array{Float64,1}, Mu::OffsetArray{Float64,1}, Mv::OffsetArray{Float64,1}, Mxi::OffsetArray{Float64,1}, 
								 alpha::Int, beta::Int )

    au = @. a[1] * moments_conserve(Mu, Mv, Mxi, alpha+0, beta+0, 0) +
         a[2] * moments_conserve(Mu, Mv, Mxi, alpha+1, beta+0, 0) +
         a[3] * moments_conserve(Mu, Mv, Mxi, alpha+0, beta+1, 0) +
         0.5 * a[4] * moments_conserve(Mu, Mv, Mxi, alpha+2, beta+0, 0) +
         0.5 * a[4] * moments_conserve(Mu, Mv, Mxi, alpha+0, beta+2, 0) +
         0.5 * a[4] * moments_conserve(Mu, Mv, Mxi, alpha+0, beta+0, 2)

    return au

end

function moments_conserve_slope( a::Array{Float64,1}, Mu::OffsetArray{Float64,1}, Mv::OffsetArray{Float64,1}, 
								 Mw::OffsetArray{Float64,1}, alpha::Int, beta::Int, delta::Int )

	au = @. a[1] * moments_conserve(Mu, Mv, Mw, alpha+0, beta+0, delta+0) +
			a[2] * moments_conserve(Mu, Mv, Mw, alpha+1, beta+0, delta+0) +
			a[3] * moments_conserve(Mu, Mv, Mw, alpha+0, beta+1, delta+0) +
			a[4] * moments_conserve(Mu, Mv, Mw, alpha+0, beta+0, delta+1) +
			0.5 * a[5] * moments_conserve(Mu, Mv, Mw, alpha+2, beta+0, delta+0) +
			0.5 * a[5] * moments_conserve(Mu, Mv, Mw, alpha+0, beta+2, delta+0) +
			0.5 * a[5] * moments_conserve(Mu, Mv, Mw, alpha+0, beta+0, delta+2)

	return au

end


"""
Velocity moments of particle distribution function
2. discrete form
"""

#--- 1D ---#
discrete_moments(f::AbstractArray{Float64,1}, u::AbstractArray{Float64,1}, ω::AbstractArray{Float64,1}, n::Int) =
sum(@. ω * u^n * f)


#--- 2D ---#
discrete_moments(f::AbstractArray{Float64,2}, u::AbstractArray{Float64,2}, ω::AbstractArray{Float64,2}, n::Int) =
sum(@. ω * u^n * f)


#--- 3D ---#
discrete_moments(f::AbstractArray{Float64,3}, u::AbstractArray{Float64,3}, ω::AbstractArray{Float64,3}, n::Int) =
sum(@. ω * u^n * f)


"""
Equilibrium in discrete form 
1. Gas: Maxwellian

# >@param[in] : particle velocity quadrature points
# >@param[in] : density, velocity and inverse of temperature

# >@return : Maxwellian distribution function
"""

#--- 1D ---#
maxwellian(u::AbstractArray{Float64,1}, ρ::Real, U::Real, λ::Real) =
@. ρ * (λ / π)^0.5 * exp(-λ * (u - U)^2)

maxwellian(u::AbstractArray{Float64,1}, prim::Array{<:Real,1}) =
maxwellian(u, prim[1], prim[2], prim[end]) # in case of input with length 4/5


#--- 2D ---#
maxwellian(u::AbstractArray{Float64,2}, v::AbstractArray{Float64,2}, ρ::Real, U::Real, V::Real, λ::Real) =
@. ρ * (λ / π) * exp(-λ * ((u - U)^2 + (v - V)^2))

maxwellian(u::AbstractArray{Float64,2}, v::AbstractArray{Float64,2}, prim::Array{<:Real,1}) =
maxwellian(u, v, prim[1], prim[2], prim[3], prim[end]) # in case of input with length 5 


#--- 3D ---#
maxwellian(u::AbstractArray{Float64,3}, v::AbstractArray{Float64,3}, w::AbstractArray{Float64,3}, 
		   ρ::Real, U::Real, V::Real, W::Real, λ::Real) =
@. ρ * (λ / π)^1.5 * exp(-λ * ((u - U)^2 + (v - V)^2+ (w - W)^2))

maxwellian(u::AbstractArray{Float64,3}, v::AbstractArray{Float64,3}, w::AbstractArray{Float64,3}, prim::Array{<:Real,1}) =
maxwellian(u, v, w, prim[1], prim[2], prim[3], prim[4], prim[5])


"""
Flow variables with conservative and primitive forms
"""

# ------------------------------------------------------------
# primitive -> conservative
# ------------------------------------------------------------
function prim_conserve(prim::Array{<:Real,1}, γ::Real)

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

prim_conserve(ρ::Real, U::Real, λ::Real, γ::Real) = 
prim_conserve([ρ, U, λ], γ)

prim_conserve(ρ::Real, U::Real, V::Real, λ::Real, γ::Real) = 
prim_conserve([ρ, U, V, λ], γ)


# ------------------------------------------------------------
# conservative -> primitive
# ------------------------------------------------------------
function conserve_prim(W::Array{<:Real,1}, γ::Real)

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

conserve_prim(ρ::Real, M::Real, E::Real, gamma::Real) = 
conserve_prim([ρ, M, E], gamma)

conserve_prim(ρ::Real, MX::Real, MY::Real, E::Real, gamma::Real) = 
conserve_prim([ρ, MX, MY, E], gamma)


"""
Thermodynamical properties
"""

# ------------------------------------------------------------
# Calculate heat capacity ratio
# ------------------------------------------------------------
function heat_capacity_ratio(K::Real, D::Int)
	
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
sound_speed(λ::Real, γ::Real) = (0.5 * γ / λ)^0.5

sound_speed(prim::Array{<:Real,1}, γ::Real) = sound_speed(prim[end], γ)


"""
Single component gas models
"""

# ------------------------------------------------------------
# Calculate reference viscosity
# 1. variable hard sphere (VHS) model
# ------------------------------------------------------------
ref_vhs_vis(Kn::Real, alpha::Real, omega::Real) = 
5.0 * (alpha + 1.) * (alpha + 2.) * √π / (4. * alpha * (5. - 2. * omega) * (7. - 2. * omega)) * Kn


# ------------------------------------------------------------
# Calculate collision time
# 1. variable hard sphere (VHS) model
# ------------------------------------------------------------
vhs_collision_time(prim::Array{<:Real,1}, muRef::Real, omega::Real) = 
muRef * 2. * prim[end]^(1. - omega) / prim[1]


"""
Multiple component gas models
"""

# ------------------------------------------------------------
# Calculate mixture collision time from AAP model
# ------------------------------------------------------------
function aap_hs_collision_time(prim::Array{<:Real,2}, mi::Real, ni::Real, me::Real, ne::Real, kn::Real)

	ν = zeros(axes(prim, 2))

	ν[1] = prim[1,1] / (mi * (ni + ne)) * 4. * sqrt(π) / 3. * sqrt(1. / prim[end,1] + 1. / prim[end,1]) / (sqrt(2.) * π * kn) + 
		   prim[1,2] / (me * (ni + ne)) * 4. * sqrt(π) / 3. * sqrt(1. / prim[end,1] + 1. / prim[end,2]) / (sqrt(2.) * π * kn)
	ν[2] = prim[1,1] / (mi * (ni + ne)) * 4. * sqrt(π) / 3. * sqrt(1. / prim[end,1] + 1. / prim[end,2]) / (sqrt(2.) * π * kn) + 
		   prim[1,2] / (me * (ni + ne)) * 4. * sqrt(π) / 3. * sqrt(1. / prim[end,2] + 1. / prim[end,2]) / (sqrt(2.) * π * kn)

	return 1. ./ ν

end


# ------------------------------------------------------------
# Calculate mixture primitive variables from AAP model
# ------------------------------------------------------------
function aap_hs_prim(prim::Array{<:Real,2}, tau::Array{<:Real,1}, mi::Real, ni::Real, me::Real, ne::Real, kn::Real)

	mixprim = similar(prim)

	mixprim[1,:] = deepcopy(prim[1,:])
	mixprim[2,1] = prim[2,1] + tau[1] / kn * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,2] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (prim[2,2] - prim[2,1])
	mixprim[2,2] = prim[2,2] + tau[2] / kn * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,1] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (prim[2,1] - prim[2,2])
	mixprim[3,1] = prim[3,1] + tau[1] / kn * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,2] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (prim[3,2] - prim[3,1])
	mixprim[3,2] = prim[3,2] + tau[2] / kn * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,1] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (prim[3,1] - prim[3,2])
	mixprim[4,1] = prim[4,1] + tau[1] / kn * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,2] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (prim[4,2] - prim[4,1])
	mixprim[4,2] = prim[4,2] + tau[2] / kn * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,1] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (prim[4,1] - prim[4,2])
	
	mixprim[5,1] = 1. / (1. / prim[end,1] - 2. / 3. * (mixprim[2,1] - prim[2,1])^2 - 2. / 3. * (mixprim[3,1] - prim[3,1])^2 - 2. / 3. * (mixprim[4,1] - prim[4,1])^2 + 
				   tau[1] / kn * 2. * mi / (mi + me) * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,2] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (1. / prim[end,2] * me / mi - 1. / prim[end,1] + 
				   2. / 3. * me / mi * (prim[2,2] - prim[2,1])^2 + 2. / 3. * me / mi * (prim[3,2] - prim[3,1])^2 + 2. / 3. * me / mi * (prim[4,2] - prim[4,1])^2))
	mixprim[5,2] = 1. / (1. / prim[end,2] - 2. / 3. * (mixprim[2,2] - prim[2,2])^2 - 2. / 3. * (mixprim[3,2] - prim[3,2])^2 - 2. / 3. * (mixprim[4,2] - prim[4,2])^2 + 
				   tau[2] / kn * 2. * me / (mi + me) * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,1] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (1. / prim[end,1] * mi / me - 1. / prim[end,2] + 
				   2. / 3. * mi / me * (prim[2,1] - prim[2,2])^2 + 2. / 3. * mi / me * (prim[3,1] - prim[3,2])^2 + 2. / 3. * mi / me * (prim[4,1] - prim[4,2])^2))

	return mixprim

end