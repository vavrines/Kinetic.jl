# ============================================================
# Methods of Particle Velocity Space
# ============================================================


export VSpace1D,
       VSpace2D,
       newton_cotes


# ------------------------------------------------------------
# Structure of velocity space
# ------------------------------------------------------------
mutable struct VSpace1D <: AbstractVelocitySpace

	u0 :: Float64
	u1 :: Float64
	nu :: Int64
	u :: Array{Float64,1}
    weights :: Array{Float64,1}

	VSpace1D() = VSpace1D(-5.0, 5.0, 50, "rectangle")

    function VSpace1D(U0::AbstractFloat, U1::AbstractFloat, UNUM::Int, TYPE="rectangle"::String)

		u0 = U0
		u1 = U1
		nu = UNUM
        du = (U1 - U0) / UNUM

		u = zeros(UNUM)
		weights = similar(u)

		if TYPE == "rectangle" #// rectangular formula
			for i=1:UNUM
				u[i] = U0 + (i - 0.5) * du
			end
			weights = ones(UNUM) .* du

		elseif TYPE == "newton" #// newton-cotes formula
			for i=1:UNUM
				u[i] = U0 + (i - 0.5) * du
				weights[i] = newton_cotes(i, UNUM) * du
            end
            
		elseif TYPE == "gauss" #// gaussian integration
            println("Gaussian integration coming soon")
            
		else
            println("error: no velocity quadrature rule")
		end

		# inner constructor method
		new(u0, u1, nu, u, weights)
    
    end # constructor

end # struct


mutable struct VSpace2D <: AbstractVelocitySpace

	u0 :: Float64
	u1 :: Float64
	nu :: Int64
	v0 :: Float64
	v1 :: Float64
	nv :: Int64
	u :: Array{Float64,2}
	v :: Array{Float64,2}
    weights :: Array{Float64,2}

	VSpace2D() = VSpace2D(-5.0, 5.0, 28, -5.0, 5.0, 28, 1)

    function VSpace2D(U0::AbstractFloat, U1::AbstractFloat, UNUM::Int, 
                          V0::AbstractFloat, V1::AbstractFloat, VNUM::Int, TYPE=1::Int)

		u0 = U0
		u1 = U1
		nu = UNUM
		v0 = V0
		v1 = V1
		nv = VNUM
        du = (U1 - U0) / UNUM
        dv = (V1 - V0) / VNUM
            
		if TYPE == "rectangle" #// rectangular formula
			u = zeros(UNUM, VNUM)
			v = zeros(UNUM, VNUM)
			for i=1:UNUM
				for j=1:VNUM
					u[i,j] = U0 + (i - 0.5) * du
					v[i,j] = V0 + (j - 0.5) * dv
				end
			end
			weights = ones(UNUM, VNUM) .* du .* dv

		elseif TYPE == "newton" #// newton-cotes formula
			u = zeros(UNUM, VNUM)
			v = zeros(UNUM, VNUM)
			weights = zeros(UNUM, VNUM)
			for i=1:UNUM
				for j=1:VNUM
					u[i,j] = U0 + (i - 0.5) * du
					v[i,j] = V0 + (j - 0.5) * dv

					weights[i,j] = (newton_cotes(i, UNUM) * du) * (newton_cotes(j, VNUM) * dv)
				end
			end

		elseif TYPE == "gauss" #// gaussian integration
			println("Gaussian integration coming soon")
		else
			println("error: no velocity quadrature rule")
		end

		# inner constructor method
		new(u0, u1, nu, v0, v1, nv, u, v, weights)
    
    end # constructor

end # struct


# ------------------------------------------------------------
# Newton-Cotes rule
# ------------------------------------------------------------
function newton_cotes(idx::Int, num::Int)

	if idx == 1 || idx == num
	    nc_coeff = 14.0 / 45.0
	elseif (idx - 5) % 4 == 0
		nc_coeff = 28.0 / 45.0
	elseif (idx - 3) % 4 == 0
	    nc_coeff = 24.0 / 45.0
	else
	    nc_coeff = 64.0 / 45.0
	end

	return nc_coeff

end