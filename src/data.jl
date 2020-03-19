# ============================================================
# Data Structures
# ============================================================


export Setup,
	   GasProperty,
	   PlasmaProperty,
       IB1D1F,
	   IB1D2F,
	   MIB1D4F,
	   ControlVolume1D1F,
	   MControlVolume1D4F,
	   Interface1D1F,
	   MInterface1D4F


# ------------------------------------------------------------
# Structure of computational setup
# ------------------------------------------------------------
struct Setup <: AbstractSetup

    case :: AbstractString
	space :: AbstractString
	nSpecies :: Int64
	interpOrder :: Int64
	limiter :: AbstractString
	cfl :: Float64
	maxTime :: Float64

	function Setup( CASE::AbstractString, SPACE::AbstractString, SPECIES::Int, ORDER::Int, 
					LM::AbstractString, CFL::Real, TIME::Real )

    	case = CASE
		space = SPACE
		nSpecies = SPECIES
		interpOrder = Int64(ORDER)
		limiter = LM
    	cfl = Float64(CFL)
    	maxTime = Float64(TIME)

		# inner constructor method
		new(case, space, nSpecies, interpOrder, limiter, cfl, maxTime)
    
    end

end


# ------------------------------------------------------------
# Structure of property
# ------------------------------------------------------------
struct GasProperty <: AbstractProperty

	Kn :: Float64; Ma :: Float64; Pr :: Float64
	K :: Float64; γ :: Float64; ω :: Float64
    αᵣ :: Float64; ωᵣ :: Float64; μᵣ :: Float64

    function GasProperty( KN::Real, MA::Real, PR::Real, 
    			 		  INK::Real, GAMMA::Real, OMEGA::Real,
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

	Kn :: Array{Float64,1}; Ma :: Float64; Pr :: Float64
	K :: Float64; γ :: Float64;
	
	mi :: Float64; ni::Float64
	me :: Float64; ne::Float64
	lD :: Float64; rL::Float64 

	sol :: Float64; χ :: Float64; ν :: Float64
	A1p :: Array{Float64,2}; A1n :: Array{Float64,2}; D1 :: Array{Float64,1}

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
struct IB1D1F <: AbstractCondition

	# initial condition
	wL :: Array{Float64,1}
	primL :: Array{Float64,1}
	hL :: AbstractArray{Float64,1}
	bcL :: Array{Float64,1}

	wR :: Array{Float64,1}
	primR :: Array{Float64,1}
	hR :: AbstractArray{Float64,1}
	bcR :: Array{Float64,1}

    function IB1D1F( WL::Array{<:Real,1}, PRIML::Array{<:Real,1}, 
    			     HL::AbstractArray{Float64,1}, BCL::Array{<:Real,1}, 
    			     WR::Array{<:Real,1}, PRIMR::Array{<:Real,1}, 
    			     HR::AbstractArray{Float64,1}, BCR::Array{<:Real,1} )

    	wL = Float64.(WL); primL = Float64.(PRIML); hL = deepcopy(HL); bcL = Float64.(BCL)
    	wR = Float64.(WR); primR = Float64.(PRIMR); hR = deepcopy(HR); bcR = Float64.(BCR)

		# inner constructor
		new(wL, primL, hL, bcL, wR, primR, hR, bcR)
    
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
mutable struct ControlVolume1D1F <: AbstractControlVolume1D

	x :: Float64
	dx :: Float64

	w :: Array{Float64,1}
	prim :: Array{Float64,1}
	sw :: Array{Float64,1}

	f :: AbstractArray{Float64,1}
	sf :: AbstractArray{Float64,1}

	function ControlVolume1D1F( X::Real, DX::Real, 
							    w0::Array{<:Real,1}, prim0::Array{<:Real,1}, f0::AbstractArray{Float64,1} )

		x = Float64(X)
		dx = Float64(DX)

		w = Float64.(w0)
		prim = Float64.(prim0)
		sw = zeros(axes(w))

		f = deepcopy(f0)
		sf = zeros(axes(f))

		new(x, dx, w, prim, sw, f, sf)

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

	function ControlVolume1D1F( X::Real, DX::Real, 
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
mutable struct Interface1D1F <: AbstractInterface1D

	fw :: Array{Float64,1}
	ff :: AbstractArray{Float64,1}

	function Interface1D1F(f::AbstractArray{Float64,1}) # for ghost cell

		fw = zeros(3)
		ff = zeros(axes(f))

		new(fw, ff)

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