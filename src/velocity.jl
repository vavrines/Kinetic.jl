# ============================================================
# Methods of Particle Velocity Space
# ============================================================


export VMesh1D,
       VMesh2D,
       newton_cotes


# ------------------------------------------------------------
# Structure of velocity space
# ------------------------------------------------------------
mutable struct VMesh1D <: AbstractVelocityMesh

	u0 :: Float64; u1 :: Float64; nu :: Int64
	u :: AbstractArray{Float64,1}; du :: AbstractArray{Float64,1}; weights :: AbstractArray{Float64,1}

	VMesh1D() = VMesh1D(-5, 5, 50)
	VMesh1D(U0::Union{Int, AbstractFloat}, U1::Union{Int, AbstractFloat}) = VMesh1D(U0, U1, 50)

	function VMesh1D( U0::Union{Int, AbstractFloat}, U1::Union{Int, AbstractFloat}, UNUM::Int, 
					  TYPE="rectangle"::String, NG=0::Int)

		u0 = Float64(U0); u1 = Float64(U1); nu = UNUM; δ = (u1 - u0) / nu
		u = OffsetArray{Float64}(undef, 1-NG:nu+NG); du = similar(u); weights = similar(u)

		if TYPE == "rectangle" #// rectangular formula
			for i in eachindex(u)
				u[i] = u0 + (i - 0.5) * δ
				du[i] = δ
				weights[i] = δ
			end

		elseif TYPE == "newton" #// newton-cotes formula
			for i in eachindex(u)
				u[i] = u0 + (i - 0.5) * δ
				du[i] = δ
				weights[i] = newton_cotes(i+NG, UNUM+NG*2) * δ
            end
            
		elseif TYPE == "gauss" #// gaussian integration
            println("Gaussian integration coming soon")
            
		else
            println("error: no velocity quadrature rule")
		end

		# inner constructor method
		new(u0, u1, nu, u, du, weights)
    
    end # constructor

end # struct


mutable struct VMesh2D <: AbstractVelocityMesh

	u0 :: Float64; u1 :: Float64; nu :: Int64
	v0 :: Float64; v1 :: Float64; nv :: Int64
	u :: Array{Float64,2}; v :: Array{Float64,2}
	du :: Array{Float64,2}; dv :: Array{Float64,2}
    weights :: Array{Float64,2}

	VMesh2D() = VMesh2D(-5, 5, 28, -5, 5, 28)
	VMesh2D(U0::Union{Int, AbstractFloat}, U1::Union{Int, AbstractFloat}, 
			V0::Union{Int, AbstractFloat}, V1::Union{Int, AbstractFloat}) = 
	VMesh2D(U0, U1, 28, V0, V1, 28)

    function VMesh2D( U0::Union{Int, AbstractFloat}, U1::Union{Int, AbstractFloat}, UNUM::Int, 
					  V0::Union{Int, AbstractFloat}, V1::Union{Int, AbstractFloat}, VNUM::Int, 
					  TYPE="rectangle"::String, NGU=0::Int, NGV=0::Int )

		u0 = Float64(U0); u1 = Float64(U1); nu = UNUM; δu = (u1 - u0) / nu
		v0 = Float64(V0); v1 = Float64(V1); nv = VNUM; δv = (v1 - v0) / nv
		u = OffsetArray{Float64}(undef, 1-NGU:nu+NGU, 1-NGV:nv+NGV)
        v = similar(u); du = similar(u); dv = similar(u); weights = similar(u)

		if TYPE == "rectangle" #// rectangular formula
			for j in axes(u, 2)
				for i in axes(u, 1)
					u[i,j] = u0 + (i - 0.5) * δu
					v[i,j] = v0 + (j - 0.5) * δv
					du[i,j] = δu
					dv[i,j] = δv
					weights[i,j] = δu * δv
				end
			end

		elseif TYPE == "newton" #// newton-cotes formula
			for j in axes(u, 2)
				for i in axes(u, 1)
					u[i,j] = u0 + (i - 0.5) * δu
					v[i,j] = v0 + (j - 0.5) * δv
					du[i,j] = δu
					dv[i,j] = δv
					weights[i,j] = newton_cotes(i+NGU, UNUM+NGU*2) * δu * newton_cotes(j+NGV, VNUM+NGV*2) * δv
				end
			end

		elseif TYPE == "gauss" #// gaussian integration
			println("Gaussian integration coming soon")

		else
			println("error: no velocity quadrature rule")
		end

		# inner constructor method
		new(u0, u1, nu, v0, v1, nv, u, v, du, dv, weights)
    
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