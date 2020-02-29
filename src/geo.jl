# ============================================================
# Geometrical Methods
# ============================================================


export PMesh1D,
       PMesh2D,
       uniform_mesh,
       global_frame, 
       local_frame


# ------------------------------------------------------------
# Structure of mesh
# ------------------------------------------------------------
mutable struct PMesh1D <: AbstractPhysicalMesh

	x0 :: Float64; x1 :: Float64; nx :: Int64
	x :: AbstractArray{Float64,1}; dx :: AbstractArray{Float64,1}

    PMesh1D() = PMesh1D(0, 1, 100)
    PMesh1D(X0::Union{Int, AbstractFloat}, X1::Union{Int, AbstractFloat}) = PMesh1D(X0, X1, 100)
    
    function PMesh1D( X0::Union{Int, AbstractFloat}, X1::Union{Int, AbstractFloat}, XNUM::Int, 
                      TYPE="uniform"::String, NG=0::Int)

        x0 = Float64(X0); x1 = Float64(X1); nx = XNUM; δ = (x1 - x0) / nx
        x = OffsetArray{Float64}(undef, 1-NG:nx+NG); dx = similar(x)

		if TYPE == "uniform" #// uniform mesh
            for i in eachindex(x)
                x[i] = x0 + (i - 0.5) * δ
                dx[i] = δ
            end
		end

		# inner constructor method
		new(x0, x1, nx, x, dx)
    
    end

end # struct


mutable struct PMesh2D <: AbstractPhysicalMesh

	x0 :: Float64; x1 :: Float64; nx :: Int64
	y0 :: Float64; y1 :: Float64; ny :: Int64
	x :: Array{Float64,2}; y :: Array{Float64,2}
    dx :: Array{Float64,2}; dy :: Array{Float64,2}

    PMesh2D() = PMesh2D(0, 1, 45, 0, 1, 45)
	PMesh2D(X0::Union{Int, AbstractFloat}, X1::Union{Int, AbstractFloat}, 
			Y0::Union{Int, AbstractFloat}, Y1::Union{Int, AbstractFloat}) = 
	PMesh2D(X0, X1, 45, Y0, Y1, 45)

    function PMesh2D( X0::Union{Int, AbstractFloat}, X1::Union{Int, AbstractFloat}, XNUM::Int, 
    				  Y0::Union{Int, AbstractFloat}, Y1::Union{Int, AbstractFloat}, YNUM::Int, 
    				  TYPE="uniform"::String, NGX=0::Int, NGY=0::Int)

		x0 = Float64(X0); x1 = Float64(X1); nx = XNUM; δx = (x1 - x0) / nx
        y0 = Float64(Y0); y1 = Float64(Y1); ny = YNUM; δy = (y1 - y0) / ny
        x = OffsetArray{Float64}(undef, 1-NGX:nx+NGX, 1-NGY:ny+NGY)
        y = OffsetArray{Float64}(undef, 1-NGX:nx+NGX, 1-NGY:ny+NGY)
        dx = OffsetArray{Float64}(undef, 1-NGX:nx+NGX, 1-NGY:ny+NGY)
        dy = OffsetArray{Float64}(undef, 1-NGX:nx+NGX, 1-NGY:ny+NGY)

		if TYPE == "uniform" # rectangular formula
            for j in axes(x, 2)
                for i in axes(x, 1)
                    x[i,j] = x0 + (i - 0.5) * δx
                    y[i,j] = y0 + (j - 0.5) * δy
                    dx[i,j] = δx
                    dy[i,j] = δy
                end
            end
		end

		# inner constructor method
		new(x0, x1, nx, y0, y1, ny, x, y, dx, dy)
    
    end

end # struct


function uniform_mesh(x0::AbstractFloat, xnum::Int, dx::AbstractFloat)

    points = zeros(xnum)
    for i=1:xnum
        points[i] = x0 + (i - 0.5) * dx
    end

    return points

end


function global_frame(w::Array{Float64,1}, cosa::AbstractFloat, sina::AbstractFloat) 
    
    if length(w) == 2
        G = [ w[1] * cosa - w[2] * sina, w[1] * sina + w[2] * cosa ]
    elseif length(w) == 4
        G = [ w[1], w[2] * cosa - w[3] * sina, w[2] * sina + w[3] * cosa, w[4] ]
    else
        println("error: local -> global")
    end

    return G

end


function local_frame(w::Array{Float64,1}, cosa::AbstractFloat, sina::AbstractFloat)
    
    if length(w) == 2
        L = [ w[1] * cosa + w[2] * sina, w[2] * cosa - w[1] * sina]
    elseif length(w) == 4
        L = [ w[1], w[2] * cosa + w[3] * sina, w[3] * cosa - w[2] * sina, w[4] ]
    else
        println("error: global -> local")
    end

    return L

end