# ============================================================
# Geometrical Methods
# ============================================================


export PSpace1D, PSpace2D, uniform_mesh, global_frame, local_frame, meshgrid


```
Structured mesh

```
struct PSpace1D <: AbstractPhysicalSpace

    x0::Float64
    x1::Float64
    nx::Int64
    x::AbstractArray{Float64,1}
    dx::AbstractArray{Float64,1}

    PSpace1D() = PSpace1D(0, 1, 100)
    PSpace1D(X0::Real, X1::Real) = PSpace1D(X0, X1, 100)

    function PSpace1D(
        X0::Real,
        X1::Real,
        XNUM::Int,
        TYPE = "uniform"::String,
        NG = 0::Int,
    )

        x0 = Float64(X0)
        x1 = Float64(X1)
        nx = XNUM
        δ = (x1 - x0) / nx
        x = OffsetArray{Float64}(undef, 1-NG:nx+NG)
        dx = similar(x)

        if TYPE == "uniform" # // uniform mesh
            for i in eachindex(x)
                x[i] = x0 + (i - 0.5) * δ
                dx[i] = δ
            end
        end

        # inner constructor method
        new(x0, x1, nx, x, dx)

    end

end # struct


struct PSpace2D <: AbstractPhysicalSpace

    x0::Float64
    x1::Float64
    nx::Int64
    y0::Float64
    y1::Float64
    ny::Int64
    x::AbstractArray{Float64,2}
    y::AbstractArray{Float64,2}
    dx::AbstractArray{Float64,2}
    dy::AbstractArray{Float64,2}

    PSpace2D() = PSpace2D(0, 1, 45, 0, 1, 45)
    PSpace2D(X0::Real, X1::Real, Y0::Real, Y1::Real) =
        PSpace2D(X0, X1, 45, Y0, Y1, 45)

    function PSpace2D(
        X0::Real,
        X1::Real,
        XNUM::Int,
        Y0::Real,
        Y1::Real,
        YNUM::Int,
        TYPE = "uniform"::String,
        NGX = 0::Int,
        NGY = 0::Int,
    )

        x0 = Float64(X0)
        x1 = Float64(X1)
        nx = XNUM
        δx = (x1 - x0) / nx
        y0 = Float64(Y0)
        y1 = Float64(Y1)
        ny = YNUM
        δy = (y1 - y0) / ny
        x = OffsetArray{Float64}(undef, 1-NGX:nx+NGX, 1-NGY:ny+NGY)
        y = OffsetArray{Float64}(undef, 1-NGX:nx+NGX, 1-NGY:ny+NGY)
        dx = OffsetArray{Float64}(undef, 1-NGX:nx+NGX, 1-NGY:ny+NGY)
        dy = OffsetArray{Float64}(undef, 1-NGX:nx+NGX, 1-NGY:ny+NGY)

        if TYPE == "uniform" # rectangular formula
            for j in axes(x, 2)
                for i in axes(x, 1)
                    x[i, j] = x0 + (i - 0.5) * δx
                    y[i, j] = y0 + (j - 0.5) * δy
                    dx[i, j] = δx
                    dy[i, j] = δy
                end
            end
        end

        # inner constructor method
        new(x0, x1, nx, y0, y1, ny, x, y, dx, dy)

    end

end # struct


function uniform_mesh(x0::Real, xnum::Int, dx::Real)

    points = zeros(xnum)
    for i = 1:xnum
        points[i] = x0 + (i - 0.5) * dx
    end

    return points

end


function global_frame(w::AbstractArray{<:Real,1}, cosa::Real, sina::Real)

    if length(w) == 2
        G = [w[1] * cosa - w[2] * sina, w[1] * sina + w[2] * cosa]
    elseif length(w) == 4
        G = [w[1], w[2] * cosa - w[3] * sina, w[2] * sina + w[3] * cosa, w[4]]
    else
        throw("local -> global: dimension dismatch")
    end

    return G

end


function global_frame(
    w::AbstractArray{<:Real,1},
    dirccos::AbstractArray{<:Real,2},
)

    if length(w) == 3
        G = [
            w[1] * dirccos[1, 1] + w[2] * dirccos[2, 1] + w[3] * dirccos[3, 1],
            w[1] * dirccos[1, 2] + w[2] * dirccos[2, 2] + w[3] * dirccos[3, 2],
            w[1] * dirccos[1, 3] + w[2] * dirccos[2, 3] + w[3] * dirccos[3, 3],
        ]
    elseif length(w) == 5
        G = [
            w[1],
            w[2] * dirccos[1, 1] + w[3] * dirccos[2, 1] + w[4] * dirccos[3, 1],
            w[2] * dirccos[1, 2] + w[3] * dirccos[2, 2] + w[4] * dirccos[3, 2],
            w[2] * dirccos[1, 3] + w[3] * dirccos[2, 3] + w[4] * dirccos[3, 3],
            w[5],
        ]
    else
        throw("local -> global: dimension dismatch")
    end

    return G

end


function local_frame(w::AbstractArray{<:Real,1}, cosa::Real, sina::Real)

    if length(w) == 2
        L = [w[1] * cosa + w[2] * sina, w[2] * cosa - w[1] * sina]
    elseif length(w) == 4
        L = [w[1], w[2] * cosa + w[3] * sina, w[3] * cosa - w[2] * sina, w[4]]
    else
        throw("global -> local: dimension dismatch")
    end

    return L

end


function local_frame(
    w::AbstractArray{<:Real,1},
    dirccos::AbstractArray{<:Real,2},
)

    if length(w) == 3
        L = [
            w[1] * dirccos[1, 1] + w[2] * dirccos[1, 2] + w[3] * dirccos[1, 3],
            w[1] * dirccos[2, 1] + w[2] * dirccos[2, 2] + w[3] * dirccos[2, 3],
            w[1] * dirccos[3, 1] + w[2] * dirccos[3, 2] + w[3] * dirccos[3, 3],
        ]
    elseif length(w) == 5
        L = [
            w[1],
            w[2] * dirccos[1, 1] + w[3] * dirccos[1, 2] + w[4] * dirccos[1, 3],
            w[2] * dirccos[2, 1] + w[3] * dirccos[2, 2] + w[4] * dirccos[2, 3],
            w[2] * dirccos[3, 1] + w[3] * dirccos[3, 2] + w[4] * dirccos[3, 3],
            w[5],
        ]
    else
        throw("global -> local: dimension dismatch")
    end

    return L

end


function meshgrid(x::AbstractArray{<:Real,1}, y::AbstractArray{<:Real,1})
    @assert ndims(x) == ndims(y) == 1

    X = [i for j in y, i in x]
    Y = [j for j in y, i in x]

    return X, Y
end


function meshgrid(
    x::AbstractArray{<:Real,1},
    y::AbstractArray{<:Real,1},
    z::AbstractArray{<:Real,1},
)
    @assert ndims(x) == ndims(y) == ndims(z) == 1

    X = [i for k in z, j in y, i in x]
    Y = [j for k in z, j in y, i in x]
    Z = [k for k in z, j in y, i in x]

    return X, Y, Z
end
