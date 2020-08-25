# ============================================================
# Quadrature Methods
# ============================================================


"""
Gauss-Legendre quadrature

`legendre_quadrature(n::Int)`

* @arg n : quadrature order (MUST be even)
* @return points : quadrature points in 3D coordinate
* @return weights : quadrature weights

"""
function legendre_quadrature(n::Int)

    pointsxyz = zeros(n * n, 3)
    weights = zeros(n * n)

    # construct Gauss quadrature
    mu, gaussweights = gausslegendre(n)

    # transform between (mu,phi) and (x,y,z)
    phi = [(k + 0.5) * pi / n for k = 0:2*n-1] # equidistance in z axis
    range = 1:n÷2 # only use upper half of the sphere as quadrature point due to pseudo 3D

    x = sqrt.(1.0 .- mu[range] .^ 2) .* cos.(phi)'
    y = sqrt.(1.0 .- mu[range] .^ 2) .* sin.(phi)'
    z = mu[range] .* ones(size(phi))'
    weights = 2.0 * pi / n * repeat(gaussweights[range], 1, 2 * n)

    # assign 
    pointsxyz[:, 1] .= x[:]
    pointsxyz[:, 2] .= y[:]
    pointsxyz[:, 3] .= z[:]
    weights = weights[:]

    return pointsxyz, weights

end


"""
Octaeder quadrature

`octa_quadrature(n::Int, slerpflag = true::Bool)`

* @arg n : quadrature order
* @arg slerpflag : flag of spherical linear interpolation
* @return points
* @return triangulation

"""
function octa_quadrature(n::Int, slerpflag = true::Bool)

    # integral range
    pt0 = [0.0, 0.0, 1.0]
    pt1 = [0.0, 1.0, 0.0]
    pt2 = [1.0, 0.0, 0.0]

    # slerp / linspace
    if slerpflag
        pts01 = slerp(pt0, pt1, n)
        pts02 = slerp(pt0, pt2, n)
    else
        pts01 = linspace(pt0, pt1, n)
        pts02 = linspace(pt0, pt2, n)
    end

    nptsoctant = Int64(n * (n + 1) / 2)
    pts = zeros(3, nptsoctant)

    # generate points in planar geometry
    counter = 0
    for i = 1:n
        if slerpflag
            if i == 1
                tmp = pts01[:, 1]
            else
                tmp = slerp(pts01[:, i], pts02[:, i], i)
            end
        else
            tmp = linspace(pts01[i], pts02[i], i)
        end
        for j = 1:i
            counter += 1
            if slerpflag
                pts[:, counter] = tmp[:, j]
            else
                pts[:, counter] = tmp[j]
            end
        end

    end

    # project points onto sphere
    for i = 1:nptsoctant
        pts[:, i] = pts[:, i] / norm(pts[:, i])
    end

    # enumerate over points and write their connectivity
    ids = zeros(Int64, n, n) # matrix that assigns an ID to the points
    nTrianglesOctant = Int64(n * n - 2 * n + 1)
    triangles = zeros(Int64, 3, nTrianglesOctant) # matrix that contains all triangles

    counter = 0
    for i = 1:n
        for j = 1:i
            counter += 1
            ids[i, j] = counter
        end
    end

    # create triangles
    counter = 0
    tmp = zeros(Int64, 1)
    for i = 1:n
        for j = 1:i-1
            tmp = [ids[i, j], ids[i, j+1], ids[i-1, j]]
            counter += 1
            triangles[:, counter] = tmp
        end
        if i < n
            for j = 1:i-1
                tmp = [ids[i, j], ids[i, j+1], ids[i+1, j+1]]
                counter += 1
                triangles[:, counter] = tmp
            end
        end
    end

    # now we have the quadrature points and triangles for a single octant
    ptsAll = deepcopy(pts)

    tmp = deepcopy(pts)
    tmp[1, :] *= -1.0
    ptsAll = deepcopy(hcat(ptsAll, tmp))

    tmp = deepcopy(pts)
    tmp[2, :] *= -1.0
    ptsAll = deepcopy(hcat(ptsAll, tmp))

    tmp = deepcopy(pts)
    tmp[1, :] *= -1.0
    tmp[2, :] *= -1.0
    ptsAll = deepcopy(hcat(ptsAll, tmp))

    tmp = deepcopy(ptsAll)
    tmp[3, :] *= -1.0
    ptsAll = deepcopy(hcat(ptsAll, tmp))

    trianglesAll = deepcopy(triangles)
    for i = 1:7
        trianglesAll = (hcat(trianglesAll, triangles .+ i .* nptsoctant))
    end
    ptsAll, triangulation = unique(ptsAll, trianglesAll)

    points = permutedims(ptsAll)
    triangulation = permutedims(triangulation)

    return points, triangulation

end
