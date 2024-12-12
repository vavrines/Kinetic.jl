using KitBase, Plots
using KitBase.ProgressMeter: @showprogress
using Base.Threads: @threads

###
# 2F2V
###

function step!(
    w::T1,
    prim::T1,
    h::T2,
    b::T2,
    fwL::T1,
    fhL::T2,
    fbL::T2,
    fwR::T1,
    fhR::T2,
    fbR::T2,
    fwD::T1,
    fhD::T2,
    fbD::T2,
    fwU::T1,
    fhU::T2,
    fbU::T2,
    u::T3,
    v::T3,
    weights::T3,
    K,
    γ,
    μᵣ,
    ω,
    a,
    Δs,
    dt,
    RES,
    AVG,
) where {T1,T2,T3}
    #--- store W^n ---#
    w_old = deepcopy(w)

    #--- calculate force ---#
    # hermite
    no = 4
    hsource = hermite_force(h, u, v, weights, prim, no, a) * dt
    bsource = hermite_force(b, u, v, weights, prim, no, a) * dt

    # difference
    #=hsource = zero(h)
    bsource = zero(b)
    for j in 2:size(u, 2)-1
        for i in 2:size(u, 1)-1
            hsource[i, j] = (h[i, j+1] - h[i, j-1]) / (2 * vs.dv[1]) * a[2] * dt
            bsource[i, j] = (b[i, j+1] - b[i, j-1]) / (2 * vs.dv[1]) * a[2] * dt
        end
    end=#

    #--- update W^{n+1} ---#
    @. w += (fwL - fwR + fwD - fwU) / Δs
    prim .= conserve_prim(w, γ)
    prim[2] = prim[2] + a[1] * dt
    prim[3] = prim[3] + a[2] * dt
    w .= prim_conserve(prim, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    #--- calculate M^{n+1} and tau^{n+1} ---#
    MH = maxwellian(u, v, prim)
    MB = MH .* K ./ (2.0 * prim[end])
    τ = vhs_collision_time(prim, μᵣ, ω)

    #--- update distribution function ---#
    for i in eachindex(u)
        h[i] +=
            (fhL[i] - fhR[i] + fhD[i] - fhU[i]) / Δs + dt / τ * (MH[i] - h[i]) - hsource[i]
        b[i] +=
            (fbL[i] - fbR[i] + fbD[i] - fbU[i]) / Δs + dt / τ * (MB[i] - b[i]) - bsource[i]
    end
end

function up!(KS, ctr, a1face, a2face, dt, residual;)
    nx, ny, dx, dy = KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy
    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for j in 1:ny
        for i in 1:nx
            step!(
                ctr[i, j].w,
                ctr[i, j].prim,
                ctr[i, j].h,
                ctr[i, j].b,
                a1face[i, j].fw,
                a1face[i, j].fh,
                a1face[i, j].fb,
                a1face[i+1, j].fw,
                a1face[i+1, j].fh,
                a1face[i+1, j].fb,
                a2face[i, j].fw,
                a2face[i, j].fh,
                a2face[i, j].fb,
                a2face[i, j+1].fw,
                a2face[i, j+1].fh,
                a2face[i, j+1].fb,
                KS.vs.u,
                KS.vs.v,
                KS.vs.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                [0.0, ϕ],
                dx[i, j] * dy[i, j],
                dt,
                sumRes,
                sumAvg,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * nx * ny) / (sumAvg[i] + 1.e-7)
    end

    return nothing
end

set = Setup(;
    case="cacity",
    space="2d2f2v",
    boundary=["maxwell", "maxwell", "maxwell", "maxwell"],
    limiter="minmod",
    cfl=0.8,
    maxTime=2,
)
ps = PSpace2D(0, 1, 16, 0, 1, 16)
vs = VSpace2D(-5, 5, 28, -5, 5, 28)
gas = Gas(; Kn=0.075, K=1.0)

ϕ = -0.3
fw = function (x, y, p)
    ρ = exp(2.0 * y * ϕ)
    prim = [ρ, 0.0, 0.0, 1.0]
    return prim_conserve(prim, gas.γ)
end
ff = function (x, y, p)
    w = fw(x, y, p)
    prim = conserve_prim(w, gas.γ)
    h = maxwellian(vs.u, vs.v, prim)
    b = energy_maxwellian(h, prim, gas.K)
    return h, b
end
bc = function (x, y, p)
    if y == 1.0
        return [1.0, 0.15, 0.0, 1.0]
    else
        return [1.0, 0.0, 0.0, 1.0]
    end
end
ib = IB2F(fw, ff, bc, ())

ks = SolverSet(set, ps, vs, gas, ib)
ctr, a1face, a2face = init_fvm(ks)

t = 0.0
dt = timestep(ks, ctr, 0.0)
nt = ks.set.maxTime ÷ dt |> Int
res = zeros(4)

@showprogress for iter in 1:nt
    evolve!(ks, ctr, a1face, a2face, dt)
    up!(ks, ctr, a1face, a2face, dt, res)

    if maximum(res) < 1e-6
        break
    end
end

plot(ks, ctr)

# heat flux
hf = zeros(ps.nx, ps.ny, 2)
for i in 1:ps.nx, j in 1:ps.ny
    hf[i, j, :] .=
        heat_flux(ctr[i, j].h, ctr[i, j].b, ctr[i, j].prim, vs.u, vs.v, vs.weights)
end

# tecplot
#sol = extract_sol(ks, ctr)
#sol = cat(sol, hf, dims=3)
#KB.write_tec(ps.x, ps.y, sol)

###
# 1F2V
###

function step!(
    w::T1,
    prim::T1,
    h::T2,
    fwL::T1,
    fhL::T2,
    fwR::T1,
    fhR::T2,
    fwD::T1,
    fhD::T2,
    fwU::T1,
    fhU::T2,
    u::T3,
    v::T3,
    weights::T3,
    K,
    γ,
    μᵣ,
    ω,
    a,
    Δs,
    dt,
    RES,
    AVG,
) where {T1,T2,T3}
    #--- store W^n ---#
    w_old = deepcopy(w)

    #--- calculate force ---#
    # hermite
    #no = 4
    #hsource = hermite_force(h, u, v, weights, prim, no, a) * dt

    # difference
    hsource = zero(h)
    for j in 2:size(u, 2)-1
        for i in 2:size(u, 1)-1
            hsource[i, j] = (h[i, j+1] - h[i, j-1]) / (2 * vs.dv[1]) * a[2] * dt
        end
    end

    #--- update W^{n+1} ---#
    @. w += (fwL - fwR + fwD - fwU) / Δs
    prim .= conserve_prim(w, γ)
    prim[2] = prim[2] + a[1] * dt
    prim[3] = prim[3] + a[2] * dt
    w .= prim_conserve(prim, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    #--- calculate M^{n+1} and tau^{n+1} ---#
    MH = maxwellian(u, v, prim)
    τ = vhs_collision_time(prim, μᵣ, ω)

    #--- update distribution function ---#
    for i in eachindex(u)
        h[i] +=
            (fhL[i] - fhR[i] + fhD[i] - fhU[i]) / Δs + dt / τ * (MH[i] - h[i]) - hsource[i]
    end
end

function up!(KS, ctr, a1face, a2face, dt, residual;)
    nx, ny, dx, dy = KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy
    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for j in 1:ny
        for i in 1:nx
            step!(
                ctr[i, j].w,
                ctr[i, j].prim,
                ctr[i, j].f,
                a1face[i, j].fw,
                a1face[i, j].ff,
                a1face[i+1, j].fw,
                a1face[i+1, j].ff,
                a2face[i, j].fw,
                a2face[i, j].ff,
                a2face[i, j+1].fw,
                a2face[i, j+1].ff,
                KS.vs.u,
                KS.vs.v,
                KS.vs.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                [0.0, ϕ],
                dx[i, j] * dy[i, j],
                dt,
                sumRes,
                sumAvg,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * nx * ny) / (sumAvg[i] + 1.e-7)
    end

    return nothing
end

set = Setup(;
    case="cacity",
    space="2d1f2v",
    boundary=["maxwell", "maxwell", "maxwell", "maxwell"],
    limiter="minmod",
    cfl=0.8,
    maxTime=2,
)
ps = PSpace2D(0, 1, 16, 0, 1, 16)
vs = VSpace2D(-5, 5, 28, -5, 5, 28)
gas = Gas(; Kn=0.075, K=0.0, γ=3.0)

ϕ = -0.5
fw = function (x, y, p)
    ρ = exp(2.0 * y * ϕ)
    prim = [ρ, 0.0, 0.0, 1.0]
    return prim_conserve(prim, gas.γ)
end
ff = function (x, y, p)
    w = fw(x, y, p)
    prim = conserve_prim(w, gas.γ)
    h = maxwellian(vs.u, vs.v, prim)
    return h
end
bc = function (x, y, p)
    if y == 1.0
        return [1.0, 0.15, 0.0, 1.0]
    else
        return [1.0, 0.0, 0.0, 1.0]
    end
end
ib = IB2F(fw, ff, bc, ())

ks = SolverSet(set, ps, vs, gas, ib)
ctr, a1face, a2face = init_fvm(ks)

t = 0.0
dt = timestep(ks, ctr, 0.0)
nt = ks.set.maxTime ÷ dt |> Int
res = zeros(4)

@showprogress for iter in 1:nt
    evolve!(ks, ctr, a1face, a2face, dt)
    up!(ks, ctr, a1face, a2face, dt, res)

    if maximum(res) < 1e-6
        break
    end
end

plot(ks, ctr)

# heat flux
hf = zeros(ps.nx, ps.ny, 2)
for i in 1:ps.nx, j in 1:ps.ny
    hf[i, j, :] .=
        heat_flux(ctr[i, j].f, ctr[i, j].prim, vs.u, vs.v, vs.weights)
end

# tecplot
#sol = extract_sol(ks, ctr)
#@. sol[:, :, 4] = 1 / sol[:, :, 4]
#sol = cat(sol, hf, dims=3)
#KB.write_tec(ps.x, ps.y, sol)
