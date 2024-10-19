using KitBase, LinearAlgebra
using KitML.Solaris.Optim
using KitBase.ProgressMeter: @showprogress

# one-cell simplification
begin
    quadratureorder = 5
    points, weights = KitBase.octa_quadrature(quadratureorder)
    nq = size(points, 1)
    L = 1
    ne = (L + 1)^2

    α = zeros(ne)
    u0 = [2.0, 0.0, 0.0, 0.0]
    m = KitBase.eval_spherharmonic(points, L)

    res = KitBase.optimize_closure(α, m, weights, u0, KitBase.maxwell_boltzmann_dual)
    u = KitBase.realizable_reconstruct(
        res.minimizer,
        m,
        weights,
        KitBase.maxwell_boltzmann_dual_prime,
    )
end

# multi-cell case
begin
    # space
    x0 = -1.5
    x1 = 1.5
    y0 = -1.5
    y1 = 1.5
    nx = 50#100
    ny = 50#100
    dx = (x1 - x0) / nx
    dy = (y1 - y0) / ny
    pspace = KitBase.PSpace2D(x0, x1, nx, y0, y1, ny)

    # time
    tEnd = 1.0
    cfl = 0.5
    dt = cfl / 2 * (dx * dy) / (dx + dy)

    # quadrature
    quadratureorder = 5
    points, weights = KitBase.octa_quadrature(quadratureorder)
    nq = size(points, 1)
    vspace = KitBase.VSpace1D{Float64,Int64,typeof(points),typeof(weights)}(
        -1.0,
        1.0,
        nq,
        points,
        zero(points),
        weights,
    )

    # particle
    SigmaS = 1.0 .* ones(ny, nx)
    SigmaA = 0.0 .* ones(ny, nx)
    SigmaT = SigmaS + SigmaA

    # moments
    L = 1
    ne = (L + 1)^2
    phi = zeros(ne, nx, ny)
    α = zeros(ne, nx, ny)
    m = KitBase.eval_spherharmonic(points, L)
end

# initial distribution
begin
    s2 = 0.03^2
    flr = 1e-4
    init_field(x, y) = max(flr, 1.0 / (4.0 * pi * s2) * exp(-(x^2 + y^2) / 4.0 / s2))
    for j in 1:nx
        for i in 1:ny
            y = y0 + (i - 0.5) * dy
            x = x0 + (j - 0.5) * dx
            # only zeroth order moment is non-zero
            phi[1, i, j] = init_field(x, y)
        end
    end
end

global t = 0.0
flux1 = zeros(ne, nx + 1, ny)
flux2 = zeros(ne, nx, ny + 1)

@showprogress for iter in 1:20
    # regularization
    Threads.@threads for j in 1:ny
        for i in 1:nx
            res = KitBase.optimize_closure(
                α[:, i, j],
                m,
                vspace.weights,
                phi[:, i, j],
                KitBase.maxwell_boltzmann_dual,
            )
            α[:, i, j] .= res.minimizer

            phi[:, i, j] .= KitBase.realizable_reconstruct(
                res.minimizer,
                m,
                vspace.weights,
                KitBase.maxwell_boltzmann_dual_prime,
            )
        end
    end

    # flux
    fη1 = zeros(nq)
    for j in 1:ny
        for i in 2:nx
            KitBase.flux_kfvs!(
                fη1,
                KitBase.maxwell_boltzmann_dual.(α[:, i-1, j]' * m)[:],
                KitBase.maxwell_boltzmann_dual.(α[:, i, j]' * m)[:],
                vspace.u[:, 1],
                dt,
            )

            for k in axes(flux1, 1)
                flux1[k, i, j] = sum(m[k, :] .* vspace.weights .* fη1)
            end
        end
    end

    fη2 = zeros(nq)
    for i in 1:nx
        for j in 2:ny
            KitBase.flux_kfvs!(
                fη2,
                KitBase.maxwell_boltzmann_dual.(α[:, i, j-1]' * m)[:],
                KitBase.maxwell_boltzmann_dual.(α[:, i, j]' * m)[:],
                vspace.u[:, 2],
                dt,
            )

            for k in axes(flux2, 1)
                flux2[k, i, j] = sum(m[k, :] .* (vspace.weights .* fη2))
            end
        end
    end

    # update
    for j in 2:ny-1
        for i in 2:nx-1
            for q in 1:1
                phi[q, i, j] =
                    phi[q, i, j] +
                    (flux1[q, i, j] - flux1[q, i+1, j]) / dx +
                    (flux2[q, i, j] - flux2[q, i, j+1]) / dy +
                    (SigmaS[i, j] * phi[q, i, j] - SigmaT[i, j] * phi[q, i, j]) * dt
            end

            for q in 2:ne
                phi[q, i, j] =
                    phi[q, i, j] +
                    (flux1[q, i, j] - flux1[q, i+1, j]) / dx +
                    (flux2[q, i, j] - flux2[q, i, j+1]) / dy +
                    (-SigmaT[i, j] * phi[q, i, j]) * dt
            end
        end
    end

    global t += dt
end

using Plots
contourf(pspace.x[1:nx, 1], pspace.y[1, 1:ny], α[1, :, :])
contourf(pspace.x[1:nx, 1], pspace.y[1, 1:ny], phi[1, :, :])
