# ============================================================
# Coffee Creamer of Quadrature
# ============================================================


"""
Spherical linear interpolation

`slerp(pt1::AbstractArray{<:Real,1}, pt2::AbstractArray{<:Real,1}, n::Int)`

"""
function slerp(pt1::AbstractArray{<:Real,1}, pt2::AbstractArray{<:Real,1}, n::Int)
    if norm(pt1 - pt2) < 1e-10 # same points
        return repeat(pt1, 1, n) # return n copies of that point
    end

    omega = acos(dot(pt1, pt2))
    t = collect(range(0, stop = 1, length = n))

    return (sin.((1 .- t) * omega) / sin.(omega))' .* pt1 +
           (sin.(t * omega) / sin.(omega))' .* pt2
end


"""
Cleaner for all duplicate (non unique) entries of quadrature points and triangles

`unique(Points::Array{Float64,2}, Triangles::Array{Int64,2})`

* @arg Points : quadrature points
* @arg Triangles : triangulation
* @return xyz & triangulation : new quadrature points and triangulation

"""
function unique(Points::Array{Float64,2}, Triangles::Array{Int64,2})

    nPoints = size(Points)[2]
    nTriangles = size(Triangles)[2]

    map = collect(1:nPoints)
    unique = ones(Int64, nPoints)

    """
    map is a vector, where each entry contains the mapping to the ID in the unique vector we want to have.
    We start with map(i) = i and then iterate through all the points
    if it turns out, that point 10 == point 4, we will set map(10) = 4 (and map(4)=4)
    if then also point 20 == point 10 == point 4, we will set map(20) == 4 as well.
    unique simply stores a boolean as for whether or not the value is unique or a duplicate.
    This does however keep the first occurrence as unique.
    Meaning for a vector
    [1 2 3 3 4 5 3]
    unique would be
    [1 1 1 0 1 1 0]
    """

    for i = 1:nPoints
        for j = i+1:nPoints
            dist = sqrt(
                (Points[1, i] - Points[1, j])^2 +
                (Points[2, i] - Points[2, j])^2 +
                (Points[3, i] - Points[3, j])^2,
            )
            if (dist < 1e-6) # equal points
                map[j] = min(map[j], i) # take the minimal ID
                unique[j] = 0 # this point is no longer unique
            end
        end
    end

    """
    We are still not done
    Consider again the vector from above
    [1 2 3 3 4 5 3]
    map was at first
    [0 1 2 3 4 5 6]
    and is now
    [0 1 2 2 4 5 2]
    meaning that we have to renumber certain entries because
    we want to get
    [0 1 2 2 3 4 2] , i.e. the "real count
    """

    uniqueCounter = 1 # the "real" count
    for i = 1:nPoints
        if unique[i] == 1
            map[i] = uniqueCounter
            uniqueCounter += 1
        else
            map[i] = map[map[i]]
        end
    end

    # now we go through the triangles and update the IDs
    triangulation = zeros(Int64, 3, nTriangles)
    for i = 1:nTriangles
        for j = 1:3
            idx = Triangles[j, i]
            triangulation[j, i] = min(idx, map[idx])
        end
    end

    # count unique points
    nQuadPoints = 0
    for i = 1:nPoints
        nQuadPoints += (unique[i] == 1)
    end

    # write unique points to XYZ
    xyz = zeros(Float64, 3, nQuadPoints)
    counter = 1
    for i = 1:nPoints
        if unique[i] == 1
            xyz[:, counter] = Points[:, i]
            counter += 1
        end
    end

    return xyz, triangulation
    
end
