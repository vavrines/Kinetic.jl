using KitBase, LinearAlgebra
using KitBase.ProgressMeter: @showprogress

cd(@__DIR__)
ps = KitBase.UnstructPSpace("../../assets/mesh/linesource.su2")

begin
    # quadrature
    quadratureorder = 10
    points, weights = KitBase.legendre_quadrature(quadratureorder)
    nq = size(points, 1)
    vs = KitBase.UnstructVSpace(-1.0, 1.0, nq, points, weights)

    # IC
    s2 = 0.03^2
    init_field(x, y) = max(1e-4, 1.0 / (4.0 * π * s2) * exp(-(x^2 + y^2) / 4.0 / s2))

    # particle
    SigmaS = ones(size(ps.cellid, 1))
    SigmaA = zeros(size(ps.cellid, 1))
    SigmaT = SigmaS + SigmaA

    # time
    tspan = (0.0, 0.3)
    cfl = 0.7
end

ctr = Array{KitBase.ControlVolumeUS1F}(undef, size(ps.cellid, 1))
for i in eachindex(ctr)
    n = Vector{Float64}[]
    for j = 1:3
        push!(
            n,
            KitBase.unit_normal(
                ps.points[ps.facePoints[ps.cellFaces[i, j], 1], :],
                ps.points[ps.facePoints[ps.cellFaces[i, j], 2], :],
            ),
        )

        if dot(ps.faceCenter[ps.cellFaces[i, j], :] .- ps.cellCenter[i, :], n[j]) < 0
            n[j] .= -n[j]
        end
    end

    phi = zeros(nq)
    phi .= init_field(ps.cellCenter[i, 1], ps.cellCenter[i, 2])

    w = sum(weights .* phi)
    dx = [
        KitBase.point_distance(
            ps.cellCenter[i, :],
            ps.points[ps.cellid[i, 1], :],
            ps.points[ps.cellid[i, 2], :],
        ),
        KitBase.point_distance(
            ps.cellCenter[i, :],
            ps.points[ps.cellid[i, 2], :],
            ps.points[ps.cellid[i, 3], :],
        ),
        KitBase.point_distance(
            ps.cellCenter[i, :],
            ps.points[ps.cellid[i, 3], :],
            ps.points[ps.cellid[i, 1], :],
        ),
    ]

    ctr[i] = KitBase.ControlVolumeUS1F(n, ps.cellCenter[i, :], dx, w, w, phi)
end

face = Array{KitBase.Interface2D1F}(undef, size(ps.facePoints, 1))
for i in eachindex(face)
    len = norm(ps.points[ps.facePoints[i, 1], :] .- ps.points[ps.facePoints[i, 2], :])
    n = KitBase.unit_normal(ps.points[ps.facePoints[i, 1], :], ps.points[ps.facePoints[i, 2], :])

    if !(-1 in ps.faceCells[i, :])
        n0 = ps.cellCenter[ps.faceCells[i, 2], :] .- ps.cellCenter[ps.faceCells[i, 1], :]
    else
        n0 = zero(n)
    end
    if dot(n, n0) < 0
        n .= -n
    end

    fw = 0.0
    ff = zeros(nq)

    face[i] = KitBase.Interface2D1F(len, n[1], n[2], fw, ff)
end

dt = 1.2 / 150 * cfl
nt = tspan[2] ÷ dt |> Int
@showprogress for iter = 1:nt
    @inbounds Threads.@threads for i in eachindex(face)
        velo = vs.u[:, 1] .* face[i].n[1] + vs.u[:, 2] .* face[i].n[2]
        if !(-1 in ps.faceCells[i, :])
            KitBase.flux_kfvs!(
                face[i].ff,
                ctr[ps.faceCells[i, 1]].f,
                ctr[ps.faceCells[i, 2]].f,
                velo,
                dt,
            )
        end
    end

    @inbounds Threads.@threads for i in eachindex(ctr)
        if ps.cellType[i] == 0
            for j = 1:3
                dirc = sign(dot(ctr[i].n[j], face[ps.cellFaces[i, j]].n))
                @. ctr[i].f -=
                    dirc * face[ps.cellFaces[i, j]].ff * face[ps.cellFaces[i, j]].len /
                    ps.cellArea[i]
            end

            integral = KitBase.discrete_moments(ctr[i].f, vs.weights)
            integral /= 4.0 * π
            @. ctr[i].f += (integral - ctr[i].f) * dt
        end
    end
end

cdata = zeros(length(ctr), 1)
for i in eachindex(ctr)
    cdata[i, 1] = sum(ctr[i].f .* vs.weights)
end
write_vtk(ps.points, ps.cellid, cdata)
