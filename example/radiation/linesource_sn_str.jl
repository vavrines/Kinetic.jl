using KitBase
using KitBase.ProgressMeter: @showprogress

begin
    # space
    x0 = -1.5
    x1 = 1.5
    y0 = -1.5
    y1 = 1.5
    nx = 100
    ny = 100
    dx = (x1 - x0) / nx
    dy = (y1 - y0) / ny

    pspace = PSpace2D(x0, x1, nx, y0, y1, ny)

    # time
    tEnd = 1.0
    cfl = 0.95

    # quadrature
    quadratureorder = 6
    points, weights = octa_quadrature(quadratureorder)
    nq = size(points, 1)

    # particle
    SigmaS = 1.0 .* ones(ny + 4, nx + 4)
    SigmaA = 0.0 .* ones(ny + 4, nx + 4)
    SigmaT = SigmaS + SigmaA
end

# initial distribution
phi = zeros(nq, nx, ny)
s2 = 0.03^2
flr = 1e-4
init_field(x, y) = max(flr, 1.0 / (4.0 * pi * s2) * exp(-(x^2 + y^2) / 4.0 / s2))
for j in 1:nx
    for i in 1:ny
        y = y0 + dy / 2 + (i - 3) * dy
        x = x0 + dx / 2 + (j - 3) * dx
        for q in 1:nq
            phi[q, i, j] = init_field(x, y) / 4.0 / π
        end
    end
end

dt = cfl / 2 * (dx * dy) / (dx + dy)
global t = 0.0

flux1 = zeros(nq, nx + 1, ny)
flux2 = zeros(nq, nx, ny + 1)

@showprogress for iter in 1:50
    for i in 2:nx, j in 1:ny
        tmp = @view flux1[:, i, j]
        flux_kfvs!(tmp, phi[:, i-1, j], phi[:, i, j], points[:, 1], dt)
    end
    for i in 1:nx, j in 2:ny
        tmp = @view flux2[:, i, j]
        flux_kfvs!(tmp, phi[:, i, j-1], phi[:, i, j], points[:, 2], dt)
    end

    for j in 1:ny, i in 1:nx
        integral = discrete_moments(phi[:, i, j], weights)
        integral *= 1.0 / 4.0 / pi

        for q in 1:nq
            phi[q, i, j] =
                phi[q, i, j] +
                (flux1[q, i, j] - flux1[q, i+1, j]) / dx +
                (flux2[q, i, j] - flux2[q, i, j+1]) / dy +
                (integral - phi[q, i, j]) * dt
        end
    end

    global t += dt
end

ρ = zeros(nx, ny)
for i in 1:nx, j in 1:ny
    ρ[i, j] = discrete_moments(phi[:, i, j], weights)
end

using Plots
contourf(pspace.x[1:nx, 1], pspace.y[1, 1:ny], ρ[:, :])
