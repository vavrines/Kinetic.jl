# ============================================================
# Geometrical Methods
# ============================================================

export global_frame, local_frame
export PSpace1D, PSpace2D, uniform_mesh, meshgrid
export UnstructMesh 
export read_mesh, mesh_connectivity_2D, mesh_center_2D, mesh_area_2D


"""
Transform local flow variables to global frame

"""
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


"""
Transform global flow variables to local frame

"""
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


# ------------------------------------------------------------
# Structured Mesh
# ------------------------------------------------------------

"""
1D physical space with structured mesh

"""
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


"""
2D Physical space with structured mesh

"""
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


"""
Generate uniform mesh

"""
function uniform_mesh(x0::Real, xnum::Int, dx::Real)

    points = zeros(xnum)
    for i = 1:xnum
        points[i] = x0 + (i - 0.5) * dx
    end

    return points

end


"""
Equivalent structured mesh generator as matlab
* 2D
* 3D

"""
function meshgrid(x::AbstractArray{<:Real,1}, y::AbstractArray{<:Real,1})
    X = [i for j in y, i in x]
    Y = [j for j in y, i in x]

    return X, Y
end


function meshgrid(
    x::AbstractArray{<:Real,1},
    y::AbstractArray{<:Real,1},
    z::AbstractArray{<:Real,1},
)
    X = [i for k in z, j in y, i in x]
    Y = [j for k in z, j in y, i in x]
    Z = [k for k in z, j in y, i in x]

    return X, Y, Z
end


# ------------------------------------------------------------
# Unstructured Mesh
# ------------------------------------------------------------

"""
Physical space with unstructured mesh

"""
struct UnstructMesh{A,B} <: AbstractPhysicalSpace

    nodes::A # locations of vertex points
    cells::B # node indices of elements

    function UnstructMesh(nodes, cells)
        new{typeof(nodes),typeof(cells)}(nodes, cells)
    end

end # struct


"""
Read mesh file
* nodes : are saved with 3D coordinates (z=0 for 2D case)
* cells : node ids inside cells
"""
function read_mesh(file)
    meshio = pyimport("meshio")
    m0 = meshio.read(file)
    nodes = m0.points
    cells = m0.cells[end][2] .+ 1 # python data is zero-indexed

    return nodes, cells
end


"""
Compute connectivity of 2D unstructured mesh

"""
function mesh_connectivity_2D(cells::AbstractArray{<:Int,2})

    nNodesPerCell = size(cells, 2)
    nCells = size(cells, 1)
    nEdgesMax = nNodesPerCell * nCells

    tmpEdgeNodes = -ones(Int, nEdgesMax, 2)
    tmpEdgeCells = -ones(Int, nEdgesMax, 2)

    counter = 0
    for i = 1:nCells, k = 1:nNodesPerCell
        isNewEdge = true
        for j = 1:counter
            if tmpEdgeNodes[j, :] == [cells[i, k], cells[i, k%nNodesPerCell+1]] ||
               tmpEdgeNodes[j, :] == [cells[i, k%nNodesPerCell+1], cells[i, k]]
                isNewEdge = false
                tmpEdgeCells[j, 2] = i
            end
        end
        if isNewEdge
            counter += 1
            tmpEdgeNodes[counter, 1] = cells[i, k]
            tmpEdgeNodes[counter, 2] = cells[i, k%nNodesPerCell+1]
            tmpEdgeCells[counter, 1] = i
        end
    end

    nEdges = counter
    edgeNodes = tmpEdgeNodes[1:nEdges, :]
    edgeCells = tmpEdgeCells[1:nEdges, :]

    cellNeighbors = -ones(Int, nCells, nNodesPerCell)
    for i = 1:nCells, k = 1:nNodesPerCell, j = 1:nEdges
        if edgeNodes[j, 1] == cells[i, k] &&
           edgeNodes[j, 2] == cells[i, k%nNodesPerCell+1] ||
           edgeNodes[j, 1] == cells[i, k%nNodesPerCell+1] && edgeNodes[j, 2] == cells[i, k]
            if edgeCells[j, 1] != i && edgeCells[j, 2] == i
                cellNeighbors[i, k] = edgeCells[j, 1]
            elseif edgeCells[j, 1] == i && edgeCells[j, 2] != i
                cellNeighbors[i, k] = edgeCells[j, 2]
            else
                throw("wrong info in neighboring cells of edge")
            end
        end
    end

    return edgeNodes, edgeCells, cellNeighbors

end


"""
Compute areas of 2D elements

"""
function mesh_area_2D(
    nodes::AbstractArray{<:AbstractFloat,2},
    cells::AbstractArray{<:Int,2},
)

    ΔS = zeros(size(cells, 1))

    if size(cells, 2) == 3 # triangular mesh
        for i in eachindex(ΔS)
            ΔS[i] = abs(
                (
                    nodes[cells[i, 1], 1] *
                    (nodes[cells[i, 2], 2] - nodes[cells[i, 3], 2]) +
                    nodes[cells[i, 2], 1] *
                    (nodes[cells[i, 3], 2] - nodes[cells[i, 1], 2]) +
                    nodes[cells[i, 3], 1] * (nodes[cells[i, 1], 2] - nodes[cells[i, 2], 2])
                ) / 2,
            )
        end
    elseif size(cells, 2) == 4 # quadrilateral mesh
        for i in eachindex(ΔS)
            d1 = [
                nodes[cells[i][1]][1] - nodes[cells[i][2]][1],
                nodes[cells[i][1]][2] - nodes[cells[i][2]][2],
            ]
            d2 = [
                nodes[cells[i][2]][1] - nodes[cells[i][3]][1],
                nodes[cells[i][2]][2] - nodes[cells[i][3]][2],
            ]
            d3 = [
                nodes[cells[i][3]][1] - nodes[cells[i][4]][1],
                nodes[cells[i][3]][2] - nodes[cells[i][4]][2],
            ]
            d4 = [
                nodes[cells[i][4]][1] - nodes[cells[i][1]][1],
                nodes[cells[i][4]][2] - nodes[cells[i][1]][2],
            ]

            a = sqrt(d1[1]^2 + d1[2]^2)
            b = sqrt(d2[1]^2 + d2[2]^2)
            c = sqrt(d3[1]^2 + d3[2]^2)
            d = sqrt(d4[1]^2 + d4[2]^2)
            T = 0.5 * (a + b + c + d)

            alpha = acos((d4[1] * d1[1] + d4[2] * d1[2]) / (a * d))
            beta = acos((d2[1] * d3[1] + d2[2] * d3[2]) / (b * c))

            ΔS[i] = sqrt(
                (T - a) * (T - b) * (T - c) * (T - d) -
                a * b * c * d * cos(0.5 * (alpha + beta)) * cos(0.5 * (alpha + beta)),
            )
        end
    end

    return ΔS

end


"""
Compute central points of 2D elements

"""
function mesh_center_2D(nodes, cells)

    cellMidPoints = zeros(size(cells, 1), 2)
    for i in axes(cellMidPoints, 1) # nCells
        for j in axes(cells, 2) # nNodesPerCell
            cellMidPoints[i, :] .+= nodes[cells[i, j], :]
        end
    end
    cellMidPoints ./= size(cells, 2)

    return cellMidPoints

end

