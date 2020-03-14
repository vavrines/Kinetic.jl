# ============================================================
# Data Structures
# ============================================================


export Setup,
       GasProperty,
       IB1D1F,
	   IB1D2F,
	   ControlVolume1D1F,
	   Interface1D1F


# ------------------------------------------------------------
# Structure of computational setup
# ------------------------------------------------------------
struct Setup <: AbstractSetup

    isNewStart :: Bool
    case :: AbstractString
    space :: AbstractString
	interpOrder :: Int64
	cfl :: Float64
	maxTime :: Float64

    function Setup(NS::Bool, CASE::AbstractString, SPACE::AbstractString, ORDER::Int, CFL::Number, TIME::Number)

        isNewStart = NS
    	case = CASE
        space = SPACE
        interpOrder = Int64(ORDER)
    	cfl = Float64(CFL)
    	maxTime = Float64(TIME)

		# inner constructor method
		new(isNewStart, case, space, interpOrder, cfl, maxTime)
    
    end

end


# ------------------------------------------------------------
# Structure of gas property
# ------------------------------------------------------------
struct GasProperty <: AbstractProperty

	Kn :: Float64; Ma :: Float64; Pr :: Float64
	K :: Float64; γ :: Float64; ω :: Float64
    αᵣ :: Float64; ωᵣ :: Float64; μᵣ :: Float64

    function GasProperty( KN::Number, MA::Number, PR::Number, 
    			 		  INK::Number, GAMMA::Number, OMEGA::Number,
    			 		  ALPHAREF::Number, OMEGAREF::Number, MUREF::Number )

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

    function IB1D1F( WL::Array{Float64,1}, PRIML::Array{Float64,1}, 
    			     HL::AbstractArray{Float64,1}, BCL::Array{Float64,1}, 
    			     WR::Array{Float64,1}, PRIMR::Array{Float64,1}, 
    			     HR::AbstractArray{Float64,1}, BCR::Array{Float64,1} )

    	wL = WL; primL = PRIML; hL = HL; bcL = BCL
    	wR = WR; primR = PRIMR; hR = HR; bcR = BCR

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

    function IB1D2F( WL::Array{Float64,1}, PRIML::Array{Float64,1}, 
    			     HL::AbstractArray{Float64,1}, BL::AbstractArray{Float64,1}, BCL::Array{Float64,1}, 
    			     WR::Array{Float64,1}, PRIMR::Array{Float64,1}, 
    			     HR::AbstractArray{Float64,1}, BR::AbstractArray{Float64,1}, BCR::Array{Float64,1} )

    	wL = WL; primL = PRIML; hL = HL; bL = BL; bcL = BCL
    	wR = WR; primR = PRIMR; hR = HR; bR = BR; bcR = BCR

		# inner constructor
		new(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR)
    
    end

end


# ------------------------------------------------------------
# Structure of control volume
# ------------------------------------------------------------
mutable struct ControlVolume1D1F <: AbstractControlVolume

	x :: Float64
	dx :: Float64

	w :: Array{Float64,1}
	prim :: Array{Float64,1}
	sw :: Array{Float64,1}

	h :: Array{Float64,1}
	sh :: Array{Float64,1}

	function ControlVolume1D1F( X::Number, DX::Number, 
							    w0::Array{<:Number,1}, prim0::Array{<:Number,1}, h0::AbstractArray{Float64,1} )

		x = X
		dx = DX

		w = w0
		prim = prim0
		sw = zeros(axes(w))

		h = h0
		sh = zeros(axes(h))

		new(x, dx, w, prim, sw, h, sh)

	end

end


# ------------------------------------------------------------
# Structure of cell interface
# ------------------------------------------------------------
mutable struct Interface1D1F <: AbstractInterface

	fw :: Array{Float64,1}
	fh :: Array{Float64,1}

	function Interface1D1F(nu::Int)

		fw = zeros(3)
		fh = zeros(nu)

		new(fw, fh)

	end

end