# ============================================================
# Theory
# ============================================================


export gauss_moments,
	   mixture_gauss_moments,
	   moments_conserve,
	   mixture_moments_conserve,
	   moments_conserve_slope,
	   mixture_moments_conserve_slope,
	   discrete_moments, 
	   maxwellian, 
	   mixture_maxwellian,
	   conserve_prim, 
	   mixture_conserve_prim,
	   prim_conserve, 
	   mixture_prim_conserve,
	   heat_capacity_ratio, 
	   ref_vhs_vis, 
	   sound_speed, 
	   vhs_collision_time,
	   aap_hs_collision_time,
	   aap_hs_prim,
	   aap_hs_diffeq,
	   shift_pdf!,
	   em_coefficients


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

	Threads.@threads for i = 2:6
		MuL[i] = prim[2] * MuL[i - 1] + 0.5 * (i - 1) * MuL[i - 2] / prim[end]
		MuR[i] = prim[2] * MuR[i - 1] + 0.5 * (i - 1) * MuR[i - 2] / prim[end]
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
		Threads.@threads for i = 2:6
			Mv[i] = prim[3] * Mv[i - 1] + 0.5 * (i - 1) * Mv[i - 2] / prim[end]
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
		Threads.@threads for i = 2:6
			Mv[i] = prim[3] * Mv[i - 1] + 0.5 * (i - 1) * Mv[i - 2] / prim[end]
		end

		Mw = OffsetArray{Float64}(undef, 0:6)
		Mw[0] = 1.
		Mw[1] = prim[4]
		Threads.@threads for i = 2:6
			Mw[i] = prim[4] * Mw[i - 1] + 0.5 * (i - 1) * Mw[i - 2] / prim[end]
		end

		return Mu, Mv, Mw, MuL, MuR
	end

end


function mixture_gauss_moments(prim::Array{<:Real,2}, inK::Real)

    Mu = OffsetArray{Float64}(undef, 0:6, axes(prim, 2)); MuL = similar(Mu); MuR = similar(Mu)
	
	if size(prim, 1) == 3
		Mxi = OffsetArray{Float64}(undef, 0:2, axes(prim, 2))
		for j in axes(prim, 2)
			_tu, _txi, _tuL, _tuR = Kinetic.gauss_moments(prim[:,j], inK)

			Mu[:,j] .= _tu
			Mxi[:,j] .= _txi
			MuL[:,j] .= _tuL
			MuR[:,j] .= _tuR
		end

		return Mu, Mxi, MuL, MuR
	elseif size(prim, 1) == 4
		Mv = OffsetArray{Float64}(undef, 0:6, axes(prim, 2))
		Mxi = OffsetArray{Float64}(undef, 0:2, axes(prim, 2))
		for j in axes(prim, 2)
			_tu, _tv, _txi, _tuL, _tuR = Kinetic.gauss_moments(prim[:,j], inK)

			Mu[:,j] .= _tu
			Mv[:,j] .= _tv
			Mxi[:,j] .= _txi
			MuL[:,j] .= _tuL
			MuR[:,j] .= _tuR
		end

		return Mu, Mv, Mxi, MuL, MuR
	elseif size(prim, 1) == 5
		Mv = OffsetArray{Float64}(undef, 0:6, axes(prim, 2))
		Mw = OffsetArray{Float64}(undef, 0:6, axes(prim, 2))
		
		for j in axes(prim, 2)
			_tu, _tv, _tw, _tuL, _tuR = Kinetic.gauss_moments(prim[:,j], inK)

			Mu[:,j] .= _tu
			Mv[:,j] .= _tv
			Mw[:,j] .= _tw
			MuL[:,j] .= _tuL
			MuR[:,j] .= _tuR
		end

		return Mu, Mv, Mw, MuL, MuR
	end

end


# ------------------------------------------------------------
# Calculate conservative moments
# ------------------------------------------------------------
function moments_conserve(Mu::OffsetArray{Float64,1}, Mxi::OffsetArray{Float64,1}, alpha::Int, delta::Int)

    uv = zeros(3)

    uv[1] = Mu[alpha] * Mxi[delta ÷ 2]
    uv[2] = Mu[alpha + 1] * Mxi[delta ÷ 2]
    uv[3] = 0.5 * (Mu[alpha + 2] * Mxi[delta ÷ 2] + Mu[alpha] * Mxi[(delta + 2) ÷ 2])

    return uv

end


function moments_conserve( Mu::OffsetArray{Float64,1}, Mv::OffsetArray{Float64,1}, Mw::OffsetArray{Float64,1}, 
						   alpha::Int, beta::Int, delta::Int )

	if length(Mw) == 3
		uv = zeros(4)

		uv[1] = Mu[alpha] * Mv[beta] * Mw[delta ÷ 2]
		uv[2] = Mu[alpha + 1] * Mv[beta] * Mw[delta ÷ 2]
		uv[3] = Mu[alpha] * Mv[beta + 1] * Mw[delta ÷ 2]
		uv[4] = 0.5 * (Mu[alpha + 2] * Mv[beta] * Mw[delta ÷ 2] + Mu[alpha] * Mv[beta + 2] * Mw[delta ÷ 2] + 
				Mu[alpha] * Mv[beta] * Mw[(delta + 2) ÷ 2])	
	else
		uv = zeros(5)

		uv[1] = Mu[alpha] * Mv[beta] * Mw[delta]
		uv[2] = Mu[alpha + 1] * Mv[beta] * Mw[delta]
		uv[3] = Mu[alpha] * Mv[beta + 1] * Mw[delta]
		uv[4] = Mu[alpha] * Mv[beta] * Mw[delta + 1]
		uv[5] = 0.5 * (Mu[alpha + 2] * Mv[beta] * Mw[delta] + Mu[alpha] * Mv[beta + 2] * Mw[delta] + 
				Mu[alpha] * Mv[beta] * Mw[delta + 2])
	end

	return uv

end


function mixture_moments_conserve( Mu::OffsetArray{Float64,2}, Mxi::OffsetArray{Float64,2},  
								   alpha::Int, beta::Int, delta::Int )

	Muv = zeros(3, axes(Mu, 2))
	for j in axes(Muv, 2)
		Muv[:,j] .= moments_conserve(Mu[:,j], Mxi[:,j], alpha, delta)
	end

	return Muv

end


function mixture_moments_conserve( Mu::OffsetArray{Float64,2}, Mv::OffsetArray{Float64,2}, Mw::OffsetArray{Float64,2}, 
								   alpha::Int, delta::Int )
	
	Muv = ifelse(length(Mw) == 3, zeros(4, axes(Mu, 2)), zeros(5, axes(Mu, 2)))
	for j in axes(Muv, 2)
		Muv[:,j] .= moments_conserve(Mu[:,j], Mv[:,j], Mw[:,j], alpha, beta, delta)
	end
	
	return Muv

end

# ------------------------------------------------------------
# Calculate slope-related conservative moments
# a = a1 + u * a2 + 0.5 * u^2 * a3
# ------------------------------------------------------------
function moments_conserve_slope(a::Array{Float64,1}, Mu::OffsetArray{Float64,1}, Mxi::OffsetArray{Float64,1}, alpha::Int)

	au = @. a[1] * moments_conserve(Mu, Mxi, alpha + 0, 0) +
         a[2] * moments_conserve(Mu, Mxi, alpha + 1, 0) +
         0.5 * a[3] * moments_conserve(Mu, Mxi, alpha + 2, 0) +
         0.5 * a[3] * moments_conserve(Mu, Mxi, alpha + 0, 2)

    return au

end


function moments_conserve_slope( a::Array{Float64,1}, Mu::OffsetArray{Float64,1}, Mv::OffsetArray{Float64,1}, Mxi::OffsetArray{Float64,1}, 
								 alpha::Int, beta::Int )

    au = @. a[1] * moments_conserve(Mu, Mv, Mxi, alpha + 0, beta + 0, 0) +
         a[2] * moments_conserve(Mu, Mv, Mxi, alpha + 1, beta + 0, 0) +
         a[3] * moments_conserve(Mu, Mv, Mxi, alpha + 0, beta + 1, 0) +
         0.5 * a[4] * moments_conserve(Mu, Mv, Mxi, alpha + 2, beta + 0, 0) +
         0.5 * a[4] * moments_conserve(Mu, Mv, Mxi, alpha + 0, beta + 2, 0) +
         0.5 * a[4] * moments_conserve(Mu, Mv, Mxi, alpha + 0, beta + 0, 2)

    return au

end

function moments_conserve_slope( a::Array{Float64,1}, Mu::OffsetArray{Float64,1}, Mv::OffsetArray{Float64,1}, 
								 Mw::OffsetArray{Float64,1}, alpha::Int, beta::Int, delta::Int )

	au = @. a[1] * moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 0, delta + 0) +
			a[2] * moments_conserve(Mu, Mv, Mw, alpha + 1, beta + 0, delta + 0) +
			a[3] * moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 1, delta + 0) +
			a[4] * moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 0, delta + 1) +
			0.5 * a[5] * moments_conserve(Mu, Mv, Mw, alpha + 2, beta + 0, delta + 0) +
			0.5 * a[5] * moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 2, delta + 0) +
			0.5 * a[5] * moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 0, delta + 2)

	return au

end


function mixture_moments_conserve_slope(a::Array{Float64,2}, Mu::OffsetArray{Float64,2}, Mxi::OffsetArray{Float64,2}, alpha::Int)

	au = zeros(3, axes(a, 2))
	for j in axes(au, 2)
		au[:,j] .= moments_conserve_slope(a[:,j], Mu[:,j], Mxi[:,j], alpha)
	end

    return au

end


function mixture_moments_conserve_slope( a::Array{Float64,2}, Mu::OffsetArray{Float64,2}, Mv::OffsetArray{Float64,2}, Mxi::OffsetArray{Float64,2}, 
										 alpha::Int, beta::Int)

	au = zeros(4, axes(a, 2))
	for j in axes(au, 2)
		au[:,j] .= moments_conserve_slope(a[:,j], Mu[:,j], Mv[:,j], Mxi[:,j], alpha, beta)
	end

    return au

end


function mixture_moments_conserve_slope( a::Array{Float64,2}, Mu::OffsetArray{Float64,2}, Mv::OffsetArray{Float64,2}, Mw::OffsetArray{Float64,2}, 
										 alpha::Int, beta::Int, delta::Int)

	au = zeros(5, axes(a, 2))
	for j in axes(au, 2)
		au[:,j] .= moments_conserve_slope(a[:,j], Mu[:,j], Mv[:,j], Mw[:,j], alpha, beta, delta)
	end

	return au

end


"""
Velocity moments of particle distribution function
2. discrete form
"""

# --- 1D ---#
discrete_moments(f::AbstractArray{Float64,1}, u::AbstractArray{Float64,1}, ω::AbstractArray{Float64,1}, n::Int) =
sum(@. ω * u^n * f)


# --- 2D ---#
discrete_moments(f::AbstractArray{Float64,2}, u::AbstractArray{Float64,2}, ω::AbstractArray{Float64,2}, n::Int) =
sum(@. ω * u^n * f)


# --- 3D ---#
discrete_moments(f::AbstractArray{Float64,3}, u::AbstractArray{Float64,3}, ω::AbstractArray{Float64,3}, n::Int) =
sum(@. ω * u^n * f)


"""
Equilibrium in discrete form 
1. Gas: Maxwellian

# >@param[in] : particle velocity quadrature points
# >@param[in] : density, velocity and inverse of temperature

# >@return : Maxwellian distribution function
"""

# --- 1D ---#
maxwellian(u::AbstractArray{Float64,1}, ρ::Real, U::Real, λ::Real) =
@. ρ * (λ / π)^0.5 * exp(-λ * (u - U)^2)

maxwellian(u::AbstractArray{Float64,1}, prim::Array{<:Real,1}) =
maxwellian(u, prim[1], prim[2], prim[end]) # in case of input with length 4/5


function mixture_maxwellian(u::AbstractArray{Float64,2}, prim::Array{<:Real,2})
	mixM = zeros(axes(u))
	for j in axes(mixM, 2)
		mixM[:,j] .= maxwellian(u[:,j], prim[:,j])
	end

	return mixM
end


# --- 2D ---#
maxwellian(u::AbstractArray{Float64,2}, v::AbstractArray{Float64,2}, ρ::Real, U::Real, V::Real, λ::Real) =
@. ρ * (λ / π) * exp(-λ * ((u - U)^2 + (v - V)^2))

maxwellian(u::AbstractArray{Float64,2}, v::AbstractArray{Float64,2}, prim::Array{<:Real,1}) =
maxwellian(u, v, prim[1], prim[2], prim[3], prim[end]) # in case of input with length 5 


function mixture_maxwellian(u::AbstractArray{Float64,3}, v::AbstractArray{Float64,3}, prim::Array{<:Real,2})
	mixM = zeros(axes(u))
	for k in axes(mixM, 3)
		mixM[:,:,k] .= maxwellian(u[:,:,k], v[:,:,k], prim[:,k])
	end

	return mixM
end


# --- 3D ---#
maxwellian(u::AbstractArray{Float64,3}, v::AbstractArray{Float64,3}, w::AbstractArray{Float64,3}, 
		   ρ::Real, U::Real, V::Real, W::Real, λ::Real) =
@. ρ * (λ / π)^1.5 * exp(-λ * ((u - U)^2 + (v - V)^2 + (w - W)^2))

maxwellian(u::AbstractArray{Float64,3}, v::AbstractArray{Float64,3}, w::AbstractArray{Float64,3}, prim::Array{<:Real,1}) =
maxwellian(u, v, w, prim[1], prim[2], prim[3], prim[4], prim[5])


function mixture_maxwellian(u::AbstractArray{Float64,4}, v::AbstractArray{Float64,4}, w::AbstractArray{Float64,4}, prim::Array{<:Real,2})
	mixM = zeros(axes(u))
	for l in axes(mixM, 4)
		mixM[:,:,:,l] .= maxwellian(u[:,:,:,l], v[:,:,:,l], w[:,:,:,l], prim[:,l])
	end

	return mixM
end


"""
Flow variables with conservative and primitive forms
"""

# ------------------------------------------------------------
# primitive -> conservative
# ------------------------------------------------------------
function prim_conserve(prim::Array{<:Real,1}, γ::Real)

	W = zeros(axes(prim))

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


function mixture_prim_conserve(prim::Array{<:Real,2}, γ::Real)
	w = zeros(axes(prim))
	for j in axes(w, 2)
		w[:,j] .= prim_conserve(prim[:,j], γ)
	end

	return w
end


# ------------------------------------------------------------
# conservative -> primitive
# ------------------------------------------------------------
function conserve_prim(W::Array{<:Real,1}, γ::Real)

	prim = zeros(axes(W))

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


function mixture_conserve_prim(w::Array{<:Real,2}, γ::Real)
	prim = zeros(axes(w))
	for j in axes(prim, 2)
		prim[:,j] .= conserve_prim(w[:,j], γ)
	end

	return prim
end


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
	elseif D == 3
		γ = (K + 5.) / (K + 3.)
	end

	return γ

end


# ------------------------------------------------------------
# Calculate speed of sound
# ------------------------------------------------------------
sound_speed(λ::Real, γ::Real) = (0.5 * γ / λ)^0.5

sound_speed(prim::Array{<:Real,1}, γ::Real) = sound_speed(prim[end], γ)

function sound_speed(prim::Array{<:Real,2}, γ::Real)
	c = zeros(axes(prim, 2))
	for j in eachindex(c)
		c[j] = sound_speed(prim[end,j], γ)
	end

	return maximum(c)
end


"""
Single component gas models
"""

# ------------------------------------------------------------
# Calculate reference viscosity
# 1. variable hard sphere (VHS) model
# ------------------------------------------------------------
ref_vhs_vis(Kn::Real, alpha::Real, omega::Real) = 
5. * (alpha + 1.) * (alpha + 2.) * √π / (4. * alpha * (5. - 2. * omega) * (7. - 2. * omega)) * Kn


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

	if size(prim, 1) == 3
		mixprim[1,:] = deepcopy(prim[1,:])
		mixprim[2,1] = prim[2,1] + tau[1] / kn * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,2] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (prim[2,2] - prim[2,1])
		mixprim[2,2] = prim[2,2] + tau[2] / kn * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,1] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (prim[2,1] - prim[2,2])
		mixprim[3,1] = 1. / (1. / prim[end,1] - 2. / 3. * (mixprim[2,1] - prim[2,1])^2 + 
					   tau[1] / kn * 2. * mi / (mi + me) * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,2] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (1. / prim[end,2] * me / mi - 1. / prim[end,1] + 
					   2. / 3. * me / mi * (prim[2,2] - prim[2,1])^2))
		mixprim[3,2] = 1. / (1. / prim[end,2] - 2. / 3. * (mixprim[2,2] - prim[2,2])^2 + 
					   tau[2] / kn * 2. * me / (mi + me) * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,1] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (1. / prim[end,1] * mi / me - 1. / prim[end,2] + 
					   2. / 3. * mi / me * (prim[2,1] - prim[2,2])^2))
	elseif size(prim, 1) == 4
		mixprim[1,:] = deepcopy(prim[1,:])
		mixprim[2,1] = prim[2,1] + tau[1] / kn * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,2] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (prim[2,2] - prim[2,1])
		mixprim[2,2] = prim[2,2] + tau[2] / kn * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,1] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (prim[2,1] - prim[2,2])
		mixprim[3,1] = prim[3,1] + tau[1] / kn * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,2] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (prim[3,2] - prim[3,1])
		mixprim[3,2] = prim[3,2] + tau[2] / kn * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,1] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (prim[3,1] - prim[3,2])
		mixprim[4,1] = 1. / (1. / prim[end,1] - 2. / 3. * (mixprim[2,1] - prim[2,1])^2 - 2. / 3. * (mixprim[3,1] - prim[3,1])^2 + 
					   tau[1] / kn * 2. * mi / (mi + me) * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,2] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (1. / prim[end,2] * me / mi - 1. / prim[end,1] + 
					   2. / 3. * me / mi * (prim[2,2] - prim[2,1])^2 + 2. / 3. * me / mi * (prim[3,2] - prim[3,1])^2))
		mixprim[4,2] = 1. / (1. / prim[end,2] - 2. / 3. * (mixprim[2,2] - prim[2,2])^2 - 2. / 3. * (mixprim[3,2] - prim[3,2])^2 + 
					   tau[2] / kn * 2. * me / (mi + me) * (4. * sqrt(2.) / (3. * sqrt(π)) * prim[1,1] / (ni + ne) / (mi + me) * sqrt(1. / prim[end,1] + 1. / prim[end,2])) * (1. / prim[end,1] * mi / me - 1. / prim[end,2] + 
					   2. / 3. * mi / me * (prim[2,1] - prim[2,2])^2 + 2. / 3. * mi / me * (prim[3,1] - prim[3,2])^2))
	elseif size(prim, 1) == 5
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
	else
		println("AAP mixture : dimension error")
	end

	return mixprim

end


# ------------------------------------------------------------
# Mixture source term function for DifferentialEquations.jl
# ------------------------------------------------------------
function aap_hs_diffeq(du, u, p, t)
	
	I₁, I₂, I₃, I₄, I₅, E₁, E₂, E₃, E₄, E₅ = u
    τᵢ, τₑ, mi, ni, me, ne, kn, γ = p

    τ = [ τᵢ, τₑ ]
    w = [ I₁ E₁;
          I₂ E₂;
          I₃ E₃;
          I₄ E₄;
          I₅ E₅ ]
    
    # modified variables
	prim = mixture_conserve_prim(w, γ)
	mixprim = aap_hs_prim(prim, τ, mi, ni, me, ne, kn)
	mixw = mixture_conserve_prim(prim, γ)

    du[1] = (mixw[1,1] - I₁) / τᵢ
    du[2] = (mixw[2,1] - I₂) / τᵢ
    du[3] = (mixw[3,1] - I₃) / τᵢ
    du[4] = (mixw[4,1] - I₄) / τᵢ
    du[5] = (mixw[5,1] - I₅) / τᵢ
    du[6] = (mixw[1,2] - E₁) / τₑ
    du[7] = (mixw[2,2] - E₂) / τₑ
    du[8] = (mixw[3,2] - E₃) / τₑ
    du[9] = (mixw[4,2] - E₄) / τₑ
    du[10] = (mixw[5,2] - E₅) / τₑ

	nothing
	
end


"""
Shift distribution function by external force
"""

# ------------------------------------------------------------
# Single-component gas
# ------------------------------------------------------------
function shift_pdf!(f::AbstractArray{Float64,1}, a::Float64, du::Float64, dt::Float64)

	q0 = eachindex(f) |> first # for OffsetArray
	q1 = eachindex(f) |> last

	if a > 0
		shift = Int(floor(a * dt / du)) # only for uniform velocity grid
		for k=q1:-1:q0+shift
			f[k] = f[k-shift]
		end
		for k=q0:shift+q0-1
			f[k] = 0.
		end

		for k=q0+1:q1
			f[k] += (dt * a - du * shift) * (f[k-1] - f[k]) / du
		end
	else
		shift = Int(floor(-a * dt / du))
		for k=q0:q1-shift
			f[k] = f[k+shift]
		end
		for k=q1-shift+1:q1
			f[k] = 0.
		end

		for k=q0:q1-1
			f[k] += (dt * a + du * shift) * (f[k] - f[k+1]) / du
		end
	end

	f[q0] = f[q0+1]
	f[q1] = f[q1-1]

end

# ------------------------------------------------------------
# Multi-component gas
# ------------------------------------------------------------
function shift_pdf!(f::AbstractArray{Float64,2}, a::Array{Float64,1}, du::Array{Float64,1}, dt::Float64)
	for j in axes(f, 2)
		shift_pdf!(f[:,j], a[j], du[j], dt)
	end
end


function em_coefficients( prim::Array{Float64,2}, E::Array{Float64,1}, B::Array{Float64,1}, mr::Float64, 
						  lD::Float64, rL::Float64, dt::Float64 )

	A = zeros(9, 9)
	A[1,1] = -1. / (2. * rL)
	A[2,2] = -1. / (2. * rL)
	A[3,3] = -1. / (2. * rL)
	A[4,1] = mr / (2. * rL)
	A[5,2] = mr / (2. * rL)
	A[6,3] = mr / (2. * rL)
	A[7,1] = 1. / (dt)
	A[8,2] = 1. / (dt)
	A[9,3] = 1. / (dt)

	A[1,4] = 1. / (dt)
	A[1,5] = -B[3] / (2. * rL)
	A[1,6] = B[2] / (2. * rL)
	A[2,4] = B[3] / (2. * rL)
	A[2,5] = 1. / (dt)
	A[2,6] = -B[1] / (2. * rL)
	A[3,4] = -B[2] / (2. * rL)
	A[3,5] = B[1] / (2. * rL)
	A[3,6] = 1. / (dt)

	A[4,7] = 1. / (dt)
	A[4,8] = mr * B[3] / (2. * rL)
	A[4,9] = -mr * B[2] / (2. * rL)
	A[5,7] = -mr * B[3] / (2. * rL)
	A[5,8] = 1. / (dt)
	A[5,9] = mr * B[1] / (2. * rL)
	A[6,7] = mr * B[2] / (2. * rL)
	A[6,8] = -mr * B[1] / (2. * rL)
	A[6,9] = 1. / (dt)

	A[7,4] = prim[1,1] / (2. * rL * lD^2)
	A[8,5] = prim[1,1] / (2. * rL * lD^2)
	A[9,6] = prim[1,1] / (2. * rL * lD^2)
	A[7,7] = -(prim[1,2] * mr) / (2. * rL * lD^2)
	A[8,8] = -(prim[1,2] * mr) / (2. * rL * lD^2)
	A[9,9] = -(prim[1,2] * mr) / (2. * rL * lD^2)

	b = zeros(9)
	b[1] = prim[2,1] / (dt) + E[1] / (2. * rL) - B[2] * prim[4,1] / (2. * rL) + B[3] * prim[3,1] / (2. * rL)
	b[2] = prim[3,1] / (dt) + E[2] / (2. * rL) - B[3] * prim[2,1] / (2. * rL) + B[1] * prim[4,1] / (2. * rL)
	b[3] = prim[4,1] / (dt) + E[3] / (2. * rL) - B[1] * prim[3,1] / (2. * rL) + B[2] * prim[2,1] / (2. * rL)
	b[4] = prim[2,2] / (dt) - mr * E[1] / (2. * rL) + mr * B[2] * prim[4,2] / (2. * rL) - mr * B[3] * prim[3,2] / (2. * rL)
	b[5] = prim[3,2] / (dt) - mr * E[2] / (2. * rL) + mr * B[3] * prim[2,2] / (2. * rL) - mr * B[1] * prim[4,2] / (2. * rL)
	b[6] = prim[4,2] / (dt) - mr * E[3] / (2. * rL) + mr * B[1] * prim[3,2] / (2. * rL) - mr * B[2] * prim[2,2] / (2. * rL)
	b[7] = E[1] / (dt) - prim[1,1] * prim[2,1] / (2. * rL * lD^2) + prim[1,2] * prim[2,2] * mr / (2. * rL * lD^2)
	b[8] = E[2] / (dt) - prim[1,1] * prim[3,1] / (2. * rL * lD^2) + prim[1,2] * prim[3,2] * mr / (2. * rL * lD^2)
	b[9] = E[3] / (dt) - prim[1,1] * prim[4,1] / (2. * rL * lD^2) + prim[1,2] * prim[4,2] * mr / (2. * rL * lD^2)

	return A, b

end