# ============================================================
# Methods of Particle Velocity Space
# ============================================================


export VSpace1D,
	   VSpace2D,
	   MVSpace1D,
	   MVSpace2D,
       newton_cotes


# ------------------------------------------------------------
# Structure of velocity space
# ------------------------------------------------------------
mutable struct VSpace1D <: AbstractVelocitySpace

	u0 :: Float64; u1 :: Float64; nu :: Int64
	u :: AbstractArray{Float64,1}; du :: AbstractArray{Float64,1}; weights :: AbstractArray{Float64,1}

	VSpace1D() = VSpace1D(-5, 5, 50)
	VSpace1D(U0::Real, U1::Real) = VSpace1D(U0, U1, 50)

	function VSpace1D( U0::Real, U1::Real, UNUM::Int, 
					   TYPE="rectangle"::String, NG=0::Int )

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


mutable struct VSpace2D <: AbstractVelocitySpace

	u0 :: Float64; u1 :: Float64; nu :: Int64
	v0 :: Float64; v1 :: Float64; nv :: Int64
	u :: Array{Float64,2}; v :: Array{Float64,2}
	du :: Array{Float64,2}; dv :: Array{Float64,2}
    weights :: Array{Float64,2}

	VSpace2D() = VSpace2D(-5, 5, 28, -5, 5, 28)
	VSpace2D(U0::Real, U1::Real, V0::Real, V1::Real) = VSpace2D(U0, U1, 28, V0, V1, 28)

    function VSpace2D( U0::Real, U1::Real, UNUM::Int, 
					   V0::Real, V1::Real, VNUM::Int, 
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
# Structure of multi-component velocity space
# ------------------------------------------------------------
mutable struct MVSpace1D <: AbstractVelocitySpace

	u0 :: Array{Float64,1}; u1 :: Array{Float64,1}; nu :: Int64
	u :: AbstractArray{Float64,2}; du :: AbstractArray{Float64,2}; weights :: AbstractArray{Float64,2}

	MVSpace1D() = MVSpace1D(-5, 5, -10, 10, 28)
	MVSpace1D(U0::Real, U1::Real, V0::Real, V1::Real) = MVSpace1D(U0, U1, V0, V1, 28)

	function MVSpace1D( Ui0::Real, Ui1::Real, Ue0::Real, Ue1::Real, 
						UNUM::Int, TYPE="rectangle"::String, NG=0::Int )

		u0 = [Ui0, Ue0]
		u1 = [Ui1, Ue1]
		nu = UNUM
		δ = (u1 .- u0) ./ nu
		u = OffsetArray{Float64}(undef, 1-NG:nu+NG, 1:2)
		du = similar(u); weights = similar(u)

		if TYPE == "rectangle" #// rectangular formula
			for j in axes(u, 2), i in axes(u, 1)
				u[i,j] = u0[j] + (i - 0.5) * δ[j]
				du[i,j] = δ[j]
				weights[i,j] = δ[j]
			end
		elseif TYPE == "newton" #// newton-cotes formula
			for j in axes(u, 2), i in axes(u, 1)
				u[i,j] = u0[j] + (i - 0.5) * δ[j]
				du[i,j] = δ[j]
				weights[i,j] = newton_cotes(i+NG, UNUM+NG*2) * δ[j]
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


mutable struct MVSpace2D <: AbstractVelocitySpace

	u0 :: Array{Float64,1}; u1 :: Array{Float64,1}; nu :: Int64
	v0 :: Array{Float64,1}; v1 :: Array{Float64,1}; nv :: Int64
	u :: AbstractArray{Float64,3}; v :: AbstractArray{Float64,3}
	du :: AbstractArray{Float64,3}; dv :: AbstractArray{Float64,3}
    weights :: AbstractArray{Float64,3}

	MVSpace2D() = MVSpace2D(-5, 5, -10, 10, 28, -5, 5, -10, 10, 28)
	MVSpace2D(U0::Real, U1::Real, V0::Real, V1::Real) = MVSpace2D(U0, U1, 28, V0, V1, 28)

    function MVSpace2D( Ui0::Real, Ui1::Real, Ue0::Real, Ue1::Real, UNUM::Int, 
					    Vi0::Real, Vi1::Real, Ve0::Real, Ve1::Real, VNUM::Int, 
					    TYPE="rectangle"::String, NGU=0::Int, NGV=0::Int )

		u0 = Float64.([Ui0, Ue0]); u1 = Float64.([Ui1, Ue1]); nu = UNUM; δu = (u1 .- u0) ./ nu
		v0 = Float64.([Vi0, Ve0]); v1 = Float64.([Vi1, Ve1]); nv = VNUM; δv = (v1 .- v0) ./ nv
		u = OffsetArray{Float64}(undef, 1-NGU:nu+NGU, 1-NGV:nv+NGV, 1:2)
		v = similar(u); du = similar(u); dv = similar(u); weights = similar(u)

		if TYPE == "rectangle" #// rectangular formula
			for k in axes(u, 3), j in axes(u, 2), i in axes(u, 1)
				u[i,j,k] = u0[k] + (i - 0.5) * δu[k]
				v[i,j,k] = v0[k] + (j - 0.5) * δv[k]
				du[i,j,k] = δu[k]
				dv[i,j,k] = δv[k]
				weights[i,j,k] = δu[k] * δv[k]
			end
		elseif TYPE == "newton" #// newton-cotes formula
			for k in axes(u, 3), j in axes(u, 2), i in axes(u, 1)
				u[i,j,k] = u0[k] + (i - 0.5) * δu[k]
				v[i,j,k] = v0[k] + (j - 0.5) * δv[k]
				du[i,j,k] = δu[k]
				dv[i,j,k] = δv[k]
				weights[i,j,k] = newton_cotes(i+NGU, UNUM+NGU*2) * δu[k] * newton_cotes(j+NGV, VNUM+NGV*2) * δv[k]
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
	elseif (idx-5) % 4 == 0
		nc_coeff = 28.0 / 45.0
	elseif (idx-3) % 4 == 0
	    nc_coeff = 24.0 / 45.0
	else
	    nc_coeff = 64.0 / 45.0
	end

	return nc_coeff

end