# ============================================================
# Data Structures
# ============================================================


export Setup,
	   GasProperty,
	   PlasmaProperty,
       IB1D1F,
	   IB1D2F,
	   ControlVolume1D1F,
	   Interface1D1F


# ------------------------------------------------------------
# Structure of computational setup
# ------------------------------------------------------------
struct Setup <: AbstractSetup

    case :: AbstractString
    space :: AbstractString
	interpOrder :: Int64
	limiter :: AbstractString
	cfl :: Float64
	maxTime :: Float64

    function Setup(CASE::AbstractString, SPACE::AbstractString, ORDER::Int, LM::AbstractString, CFL::Real, TIME::Real)

    	case = CASE
        space = SPACE
		interpOrder = Int64(ORDER)
		limiter = LM
    	cfl = Float64(CFL)
    	maxTime = Float64(TIME)

		# inner constructor method
		new(case, space, interpOrder, limiter, cfl, maxTime)
    
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

	Kn :: Float64; Ma :: Float64; Pr :: Float64
	K :: Float64; γ :: Float64; 

	sol :: Float64; χ :: Float64; ν :: Float64
	A1p :: Array{Float64,2}; A1n :: Array{Float64,2}; D1 :: Array{Float64,1}

    function GasProperty( KN::Array{<:Real,1}, MA::Real, PR::Real, 
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


# ------------------------------------------------------------
# Structure of control volume
# ------------------------------------------------------------
mutable struct ControlVolume1D1F <: AbstractControlVolume1D

	x :: Float64
	dx :: Float64

	w :: Array{Float64,1}
	prim :: Array{Float64,1}
	sw :: Array{Float64,1}

	h :: AbstractArray{Float64,1}
	sh :: AbstractArray{Float64,1}

	function ControlVolume1D1F( X::Real, DX::Real, 
							    w0::Array{<:Real,1}, prim0::Array{<:Real,1}, h0::AbstractArray{Float64,1} )

		x = Float64(X)
		dx = Float64(DX)

		w = Float64.(w0)
		prim = Float64.(prim0)
		sw = zeros(axes(w))

		h = deepcopy(h0)
		sh = zeros(axes(h))

		new(x, dx, w, prim, sw, h, sh)

	end

end


# ------------------------------------------------------------
# Structure of cell interface
# ------------------------------------------------------------
mutable struct Interface1D1F <: AbstractInterface1D

	fw :: Array{Float64,1}
	fh :: AbstractArray{Float64,1}

	function Interface1D1F(f::AbstractArray{Float64,1})

		fw = zeros(3)
		fh = zeros(axes(f))

		new(fw, fh)

	end

end