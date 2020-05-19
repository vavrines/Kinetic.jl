# ============================================================
# Data Structures
# ============================================================


export Setup,
	   GasProperty,
	   PlasmaProperty,
       IB1D1F,
	   IB1D2F,
	   MIB1D1F,
	   MIB1D2F,
	   MIB1D4F,
	   ControlVolume1D1F,
	   ControlVolume1D2F,
	   MControlVolume1D1F,
	   MControlVolume1D2F,
	   MControlVolume1D4F,
	   Interface1D1F,
	   Interface1D2F,
	   MInterface1D1F,
	   MInterface1D2F,
	   MInterface1D4F,
	   Solution1D1F,
	   Solution1D2F


# ------------------------------------------------------------
# Structure of computational setup
# ------------------------------------------------------------
struct Setup{S, I, E, F} <: AbstractSetup

    case :: S
	space :: S
	nSpecies :: I
	interpOrder :: I
	limiter :: S
	cfl :: E
	maxTime :: F

	function Setup( case::AbstractString, space::AbstractString, nSpecies::Int, interpOrder::Int, 
					limiter::AbstractString, cfl::Real, maxTime::Real )

		new{typeof(case), typeof(nSpecies), typeof(cfl), typeof(maxTime)}(
			case, space, nSpecies, interpOrder, limiter, cfl, maxTime)
    
    end

end


# ------------------------------------------------------------
# Structure of property
# ------------------------------------------------------------
struct GasProperty <: AbstractProperty

	Kn :: Float64
	Ma :: Float64
	Pr :: Float64
	K :: Float64
	γ :: Float64
	ω :: Float64
	αᵣ :: Float64
	ωᵣ :: Float64
	μᵣ :: Float64

    function GasProperty( KN::Real, MA::Real, PR::Real, INK::Real, GAMMA::Real, OMEGA::Real,
    			 		  ALPHAREF::Real, OMEGAREF::Real, MUREF::Real )

    	Kn = Float64(KN)
    	Ma = Float64(MA)
    	Pr = Float64(PR)

		K = Float64(INK)
		γ = Float64(GAMMA)
		ω = Float64(OMEGA)

		αᵣ = Float64(ALPHAREF)
		ωᵣ = Float64(OMEGAREF)
		μᵣ = Float64(MUREF)

		# inner constructor method
		new(Kn, Ma, Pr, K, γ, ω, αᵣ, ωᵣ, μᵣ)
    
    end

end


struct PlasmaProperty <: AbstractProperty

	Kn :: Array{Float64,1}
	Ma :: Float64
	Pr :: Float64
	K :: Float64
	γ :: Float64
	
	mi :: Float64
	ni::Float64
	me :: Float64
	ne::Float64
	lD :: Float64
	rL::Float64 

	sol :: Float64
	χ :: Float64
	ν :: Float64
	A1p :: Array{Float64,2}
	A1n :: Array{Float64,2}
	D1 :: Array{Float64,1}

    function PlasmaProperty( KN::Array{<:Real,1}, MA::Real, PR::Real, 
    			 		  	 INK::Real, GAMMA::Real, 
    			 		  	 MI::Real, NI::Real, ME::Real, NE::Real,
    			 		  	 LD::Real, RL::Real, SOL::Real, CHI::Real, NU::Real )

    	Kn = Float64.(KN)
    	Ma = Float64(MA)
    	Pr = Float64(PR)
		K = Float64(INK)
		γ = Float64(GAMMA)

		mi = Float64(MI); ni = Float64(NI)
		me = Float64(ME); ne = Float64(NE)
		lD = Float64(LD); rL = Float64(RL)

		sol = Float64(SOL)
		χ = Float64(CHI)
		ν = Float64(NU)
		
		# A^+
		A1p = Array{Float64}(undef, 8, 8)
		A1p[1,1] = (sol * χ) / 2.
		A1p[7,1] = χ / 2.
		A1p[2,2] = sol / 2.
		A1p[6,2] = 0.5
		A1p[3,3] = sol / 2.
		A1p[5,3] = -1. / 2.
		A1p[4,4] = (sol * ν) / 2.
		A1p[8,4] = (sol^2 * ν) / 2.
		A1p[3,5] = -sol^2 / 2.
		A1p[5,5] = sol / 2.
		A1p[2,6] = sol^2 / 2.
		A1p[6,6] = sol / 2.
		A1p[1,7] = (sol^2 * χ) / 2.
		A1p[7,7] = (sol * χ) / 2.
		A1p[4,8] = ν / 2.
		A1p[8,8] = (sol * ν) / 2.

		# A^-
		A1n = Array{Float64}(undef, 8, 8)
		A1n[1,1] = -(sol * χ) / 2.
		A1n[7,1] = χ / 2.
		A1n[2,2] = -sol / 2.
		A1n[6,2] = 1. / 2.
		A1n[3,3] = -sol / 2.
		A1n[5,3] = -1. / 2.
		A1n[4,4] = -(sol * ν) / 2.
		A1n[8,4] = (sol^2 * ν) / 2.
		A1n[3,5] = -sol^2 / 2.
		A1n[5,5] = -sol / 2.
		A1n[2,6] = sol^2 / 2.
		A1n[6,6] = -sol / 2.
		A1n[1,7] = (sol^2 * χ) / 2.
		A1n[7,7] = -(sol * χ) / 2.
		A1n[4,8] = ν / 2.
		A1n[8,8] = -(sol * ν) / 2.

		D1 = [sol, sol, sol * χ, sol * ν, -sol, -sol, -sol * χ, -sol * ν]

		# inner constructor method
		new(Kn, Ma, Pr, K, γ, mi, ni, me, ne, lD, rL, sol, χ, ν, A1p, A1n, D1)
    
    end

end


# ------------------------------------------------------------
# Structure of initial and boundary conditions
# ------------------------------------------------------------
struct IB1D1F{A,B} <: AbstractCondition

	wL :: A
	primL :: A
	fL :: B
	bcL :: A

	wR :: A
	primR :: A
	fR :: B
	bcR :: A

	# 1D1F1V
    function IB1D1F( wL::Array{<:Real,1}, primL::Array{<:Real,1}, 
    			     fL::AbstractArray{<:Real,1}, bcL::Array{<:Real,1}, 
    			     wR::Array{<:Real,1}, primR::Array{<:Real,1}, 
    			     fR::AbstractArray{<:Real,1}, bcR::Array{<:Real,1} )

		new{typeof(wL), typeof(fL)}(wL, primL, fL, bcL, wR, primR, fR, bcR)
    
	end
	
	# 1D1F3V
	function IB1D1F( wL::Array{<:Real,1}, primL::Array{<:Real,1}, 
					 fL::AbstractArray{Float64,3}, bcL::Array{<:Real,1}, 
					 wR::Array{<:Real,1}, primR::Array{<:Real,1}, 
					 fR::AbstractArray{Float64,3}, bcR::Array{<:Real,1} )

		new{typeof(wL), typeof(fL)}(wL, primL, fL, bcL, wR, primR, fR, bcR)

	end

end

struct IB1D2F <: AbstractCondition

	# initial condition
	wL :: Array{Float64,1}
	primL :: Array{Float64,1}
    hL :: AbstractArray{Float64,1}
    bL :: AbstractArray{Float64,1}
	bcL :: Array{Float64,1}

	wR :: Array{Float64,1}
	primR :: Array{Float64,1}
    hR :: AbstractArray{Float64,1}
    bR :: AbstractArray{Float64,1}
	bcR :: Array{Float64,1}

    function IB1D2F( WL::Array{<:Real,1}, PRIML::Array{<:Real,1}, 
    			     HL::AbstractArray{Float64,1}, BL::AbstractArray{<:Real,1}, BCL::Array{<:Real,1}, 
    			     WR::Array{<:Real,1}, PRIMR::Array{<:Real,1}, 
    			     HR::AbstractArray{Float64,1}, BR::AbstractArray{<:Real,1}, BCR::Array{<:Real,1} )

    	wL = Float64.(WL); primL = Float64.(PRIML); hL = deepcopy(HL); bL = Float64.(BL); bcL = Float64.(BCL)
    	wR = Float64.(WR); primR = Float64.(PRIMR); hR = deepcopy(HR); bR = Float64.(BR); bcR = Float64.(BCR)

		# inner constructor
		new(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR)
    
    end

end


struct MIB1D1F <: AbstractCondition

	# initial/boundary condition
	wL :: Array{Float64,2}
	primL :: Array{Float64,2}
    fL :: AbstractArray{Float64,2}
	bcL :: Array{Float64,2}

	wR :: Array{Float64,2}
	primR :: Array{Float64,2}
    fR :: AbstractArray{Float64,2}
	bcR :: Array{Float64,2}

	function MIB1D1F( WL::Array{<:Real,2}, PRIML::Array{<:Real,2}, FL::AbstractArray{Float64,2}, BCL::Array{<:Real,2}, 
					  WR::Array{<:Real,2}, PRIMR::Array{<:Real,2}, FR::AbstractArray{Float64,2}, BCR::Array{<:Real,2} )

		wL = Float64.(WL); primL = Float64.(PRIML); fL = deepcopy(FL); bcL = Float64.(BCL); 
		wR = Float64.(WR); primR = Float64.(PRIMR); fR = deepcopy(FR); bcR = Float64.(BCR); 

		# inner constructor
		new(wL, primL, fL, bcL, wR, primR, fR, bcR)
    
    end

end


struct MIB1D2F <: AbstractCondition

	# initial/boundary condition
	wL :: Array{Float64,2}
	primL :: Array{Float64,2}
    hL :: AbstractArray{Float64,2}
	bL :: AbstractArray{Float64,2}
	bcL :: Array{Float64,2}

	wR :: Array{Float64,2}
	primR :: Array{Float64,2}
    hR :: AbstractArray{Float64,2}
	bR :: AbstractArray{Float64,2}
	bcR :: Array{Float64,2}

	function MIB1D2F( WL::Array{<:Real,2}, PRIML::Array{<:Real,2}, HL::AbstractArray{Float64,2}, 
					  BL::AbstractArray{Float64,2}, BCL::Array{<:Real,2}, 
					  WR::Array{<:Real,2}, PRIMR::Array{<:Real,2}, HR::AbstractArray{Float64,2}, 
					  BR::AbstractArray{Float64,2}, BCR::Array{<:Real,2} )

						
		wL = Float64.(WL); primL = Float64.(PRIML); hL = deepcopy(HL); bL = Float64.(BL); bcL = Float64.(BCL); 
		wR = Float64.(WR); primR = Float64.(PRIMR); hR = deepcopy(HR); bR = Float64.(BR); bcR = Float64.(BCR); 

		# inner constructor
		new(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR)
    
    end

end


struct MIB1D4F <: AbstractCondition

	# initial/boundary condition
	wL :: Array{Float64,2}
	primL :: Array{Float64,2}
    h0L :: AbstractArray{Float64,2}
	h1L :: AbstractArray{Float64,2}
	h2L :: AbstractArray{Float64,2}
	h3L :: AbstractArray{Float64,2}
	bcL :: Array{Float64,2}
	EL :: Array{Float64,1}
	BL :: Array{Float64,1}
	lorenzL :: Array{Float64,2}

	wR :: Array{Float64,2}
	primR :: Array{Float64,2}
    h0R :: AbstractArray{Float64,2}
	h1R :: AbstractArray{Float64,2}
	h2R :: AbstractArray{Float64,2}
	h3R :: AbstractArray{Float64,2}
	bcR :: Array{Float64,2}
	ER :: Array{Float64,1}
	BR :: Array{Float64,1}
	lorenzR :: Array{Float64,2}

	function MIB1D4F( WL::Array{<:Real,2}, PRIML::Array{<:Real,2}, H0L::AbstractArray{Float64,2}, 
					  H1L::AbstractArray{Float64,2}, H2L::AbstractArray{Float64,2}, H3L::AbstractArray{Float64,2}, 
					  BCL::Array{<:Real,2}, EFIELDL::Array{<:Real,1}, BFIELDL::Array{<:Real,1}, LL::Array{<:Real,2}, 
					  WR::Array{<:Real,2}, PRIMR::Array{<:Real,2}, H0R::AbstractArray{Float64,2}, 
					  H1R::AbstractArray{Float64,2}, H2R::AbstractArray{Float64,2}, H3R::AbstractArray{Float64,2}, 
					  BCR::Array{<:Real,2}, EFIELDR::Array{<:Real,1}, BFIELDR::Array{<:Real,1}, LR::Array{<:Real,2} )

						
		wL = Float64.(WL); primL = Float64.(PRIML); h0L = deepcopy(H0L); h1L = Float64.(H1L); h2L = Float64.(H2L); h3L = Float64.(H3L)
		bcL = Float64.(BCL); EL = Float64.(EFIELDL); BL = Float64.(BFIELDL); lorenzL = Float64.(LL)
		wR = Float64.(WR); primR = Float64.(PRIMR); h0R = deepcopy(H0R); h1R = Float64.(H1R); h2R = Float64.(H2R); h3R = Float64.(H3R); 
		bcR = Float64.(BCR); ER = Float64.(EFIELDR); BR = Float64.(BFIELDR); lorenzR = Float64.(LR)

		# inner constructor
		new( wL, primL, h0L, h1L, h2L, h3L, bcL, EL, BL, lorenzL, 
			 wR, primR, h0R, h1R, h2R, h3R, bcR, ER, BR, lorenzR )
    
    end

end


# ------------------------------------------------------------
# Structure of control volume
# ------------------------------------------------------------
mutable struct ControlVolume1D1F{F,A,B} <: AbstractControlVolume1D

	x :: F
	dx :: F

	w :: A
	prim :: A
	sw :: A

	f :: B
	sf :: B

	function ControlVolume1D1F(X::Real, DX::Real, W::Array, PRIM::Array, F::AbstractArray)

		x = deepcopy(X)
		dx = deepcopy(DX)

		w = deepcopy(W)
		prim = deepcopy(PRIM)
		sw = zeros(typeof(W[1]), axes(w))

		f = deepcopy(F)
		sf = zeros(typeof(F[1]), axes(f))

		new{typeof(x),typeof(w),typeof(f)}(x, dx, w, prim, sw, f, sf)

	end

end


mutable struct ControlVolume1D2F <: AbstractControlVolume1D

	x :: Float64
	dx :: Float64

	w :: Array{Float64,1}
	prim :: Array{Float64,1}
	sw :: Array{Float64,1}

	h :: AbstractArray{Float64,1}
	b :: AbstractArray{Float64,1}
	sh :: AbstractArray{Float64,1}
	sb :: AbstractArray{Float64,1}

	function ControlVolume1D2F( X::Real, DX::Real, 
								w0::Array{<:Real,1}, prim0::Array{<:Real,1}, 
								h0::AbstractArray{Float64,1}, b0::AbstractArray{Float64,1} )

		x = Float64(X)
		dx = Float64(DX)

		w = Float64.(w0)
		prim = Float64.(prim0)
		sw = zeros(axes(w))

		h = deepcopy(h0)
		b = deepcopy(b0)
		sh = zeros(axes(h))
		sb = zeros(axes(b))

		new(x, dx, w, prim, sw, h, b, sh, sb)

	end

end


mutable struct MControlVolume1D1F <: AbstractControlVolume1D

	x :: Float64
	dx :: Float64

	w :: Array{Float64,2}
	prim :: Array{Float64,2}
	sw :: Array{Float64,2}

	f :: AbstractArray{Float64,2}
	sf :: AbstractArray{Float64,2}

	function MControlVolume1D1F( X::Real, DX::Real, 
								 w0::Array{<:Real,2}, prim0::Array{<:Real,2}, 
								 F0::AbstractArray{Float64,2} )

		x = Float64(X)
		dx = Float64(DX)

		w = Float64.(w0)
		prim = Float64.(prim0)
		sw = zeros(axes(w))

		f = deepcopy(F0)
		sf = zeros(axes(F0))

		new(x, dx, w, prim, sw, f, sf)

	end

end


mutable struct MControlVolume1D2F <: AbstractControlVolume1D

	x :: Float64
	dx :: Float64

	w :: Array{Float64,2}
	prim :: Array{Float64,2}
	sw :: Array{Float64,2}

	h :: AbstractArray{Float64,2}
	b :: AbstractArray{Float64,2}
	sh :: AbstractArray{Float64,2}
	sb :: AbstractArray{Float64,2}

	function MControlVolume1D2F( X::Real, DX::Real, 
								 w0::Array{<:Real,2}, prim0::Array{<:Real,2}, 
								 H0::AbstractArray{Float64,2}, B0::AbstractArray{Float64,2} )

		x = Float64(X)
		dx = Float64(DX)

		w = Float64.(w0)
		prim = Float64.(prim0)
		sw = zeros(axes(w))

		h = deepcopy(H0)
		b = deepcopy(B0)
		sh = zeros(axes(H0))
		sb = zeros(axes(B0))

		new(x, dx, w, prim, sw, h, b, sh, sb)

	end

end


mutable struct MControlVolume1D4F <: AbstractControlVolume1D

	x :: Float64
	dx :: Float64

	w :: Array{Float64,2}
	prim :: Array{Float64,2}
	sw :: Array{Float64,2}

	h0 :: AbstractArray{Float64,2}; h1 :: AbstractArray{Float64,2}
	h2 :: AbstractArray{Float64,2}; h3 :: AbstractArray{Float64,2}
	sh0 :: AbstractArray{Float64,2}; sh1 :: AbstractArray{Float64,2}
	sh2 :: AbstractArray{Float64,2}; sh3 :: AbstractArray{Float64,2}

	E :: Array{Float64,1}; B :: Array{Float64,1}
	ϕ :: Float64; ψ :: Float64
	lorenz :: Array{Float64,2}

	function MControlVolume1D4F( X::Real, DX::Real, 
								 w0::Array{<:Real,2}, prim0::Array{<:Real,2}, 
								 H0::AbstractArray{Float64,2}, H1::AbstractArray{Float64,2}, 
								 H2::AbstractArray{Float64,2}, H3::AbstractArray{Float64,2},
							 	 E0::Array{Float64,1}, B0::Array{Float64,1}, L::Array{Float64,2} )

		x = Float64(X)
		dx = Float64(DX)

		w = Float64.(w0)
		prim = Float64.(prim0)
		sw = zeros(axes(w))

		h0 = deepcopy(H0); h1 = deepcopy(H1); h2 = deepcopy(H2); h3 = deepcopy(H3)
		sh0 = zeros(axes(H0)); sh1 = zeros(axes(H1)); sh2 = zeros(axes(H2)); sh3 = zeros(axes(H3))

		E = deepcopy(E0); B = deepcopy(B0); 
		ϕ = 0.; ψ = 0.
		lorenz = deepcopy(L)

		new(x, dx, w, prim, sw, h0, h1, h2, h3, sh0, sh1, sh2, sh3, E, B, ϕ, ψ, lorenz)

	end

end


# ------------------------------------------------------------
# Structure of cell interface
# ------------------------------------------------------------
mutable struct Interface1D1F{A, B} <: AbstractInterface1D

	fw :: A
	ff :: B

	function Interface1D1F(f::AbstractArray{<:Real,1}) # for ghost cell

		fw = zeros(typeof(f[1]), 3)
		ff = zeros(typeof(f[1]), axes(f))

		new{typeof(fw), typeof(ff)}(fw, ff)

	end

	function Interface1D1F(f::AbstractArray{<:Real,3})

		fw = zeros(typeof(f[1]), 5)
		ff = zeros(typeof(f[1]), axes(f))

		new{typeof(fw), typeof(ff)}(fw, ff)

	end

end


mutable struct Interface1D2F <: AbstractInterface1D

	fw :: Array{Float64,1}
	fh :: AbstractArray{Float64,1}
	fb :: AbstractArray{Float64,1}

	function Interface1D2F(f::AbstractArray{Float64,1}) # for ghost cell

		fw = zeros(3)
		fh = zeros(axes(f))
		fb = zeros(axes(f))

		new(fw, fh, fb)

	end

end


mutable struct MInterface1D1F <: AbstractInterface1D

	fw :: Array{Float64,2}
	ff :: AbstractArray{Float64,2}

	function MInterface1D1F(f::AbstractArray{Float64,2})

		fw = zeros(5, axes(f, 2))
		ff = zeros(axes(f))

		new(fw, ff)

	end

end


mutable struct MInterface1D2F <: AbstractInterface1D

	fw :: Array{Float64,2}
	fh :: AbstractArray{Float64,2}
	fb :: AbstractArray{Float64,2}

	function MInterface1D2F(f::AbstractArray{Float64,2})

		fw = zeros(5, axes(f, 2))
		fh = zeros(axes(f))
		fb = zeros(axes(f))

		new(fw, fh, fb)

	end

end


mutable struct MInterface1D4F <: AbstractInterface1D

	fw :: Array{Float64,2}
	fh0 :: AbstractArray{Float64,2}
	fh1 :: AbstractArray{Float64,2}
	fh2 :: AbstractArray{Float64,2}
	fh3 :: AbstractArray{Float64,2}
	femL :: Array{Float64,1}
	femR :: Array{Float64,1}

	function MInterface1D4F(f::AbstractArray{Float64,2})

		fw = zeros(5, axes(f, 2))
		fh0 = zeros(axes(f))
		fh1 = zeros(axes(f))
		fh2 = zeros(axes(f))
		fh3 = zeros(axes(f))
		femL = zeros(8)
		femR = zeros(8)

		new(fw, fh0, fh1, fh2, fh3, femL, femR)

	end

end


mutable struct Solution1D1F{A,B} <: AbstractSolution

	w :: A
	prim :: A
	sw :: A
	f :: B
	sf :: B

	function Solution1D1F(w::Array, prim::Array, f::AbstractArray)
		sw = zeros(typeof(w[1]), axes(w))
		sf = zeros(typeof(f[1]), axes(f))

		new{typeof(w), typeof(f)}(w, prim, sw, f, sf)
	end

	function Solution1D1F(w::Array, prim::Array, sw::Array, f::AbstractArray, sf::AbstractArray)
		new{typeof(w), typeof(f)}(w, prim, sw, f, sf)
	end

end


mutable struct Solution1D2F{A,B} <: AbstractSolution

	w :: A
	prim :: A
	sw :: A
	h :: B
	b :: B
	sh :: B
	sb :: B

	function Solution1D2F(w::Array, prim::Array, h::AbstractArray, b::AbstractArray)
		sw = zeros(typeof(w[1]), axes(w))
		sh = zeros(typeof(h[1]), axes(h))
		sb = zeros(typeof(b[1]), axes(b))

		new{typeof(w), typeof(h)}(w, prim, sw, h, b, sh, sb)
	end

	function Solution1D2F(w::Array, prim::Array, sw::Array, h::AbstractArray, b::AbstractArray, sh::AbstractArray, sb::AbstractArray)
		new{typeof(w), typeof(h)}(w, prim, sw, h, b, sh, sb)
	end

end