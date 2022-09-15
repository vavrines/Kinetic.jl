using KitBase
using KitBase.OffsetArrays
using KitBase.ProgressMeter: @showprogress

cd(@__DIR__)
D = read_dict("shear.txt")
for key in keys(D)
    s = Symbol(key)
    @eval $s = $(D[key])
end

begin
    γ = heat_capacity_ratio(inK, 2)
    set = set_setup(D)
    pSpace = PSpace1D(x0, x1, nx, nxg)
    vSpace =
        VSpace2D(umin, umax, nu, vmin, vmax, nv; type = vMeshType, ngu = nug, ngv = nvg)

    μᵣ = ref_vhs_vis(knudsen, alphaRef, omegaRef)
    gas = set_property(D)

    primL = [1.0, 0.0, 1.0, 1.0]
    primR = [1.0, 0.0, -1.0, 2.0]

    wL = prim_conserve(primL, γ)
    wR = prim_conserve(primR, γ)

    HL = maxwellian(vSpace.u, vSpace.v, primL)
    HR = maxwellian(vSpace.u, vSpace.v, primR)

    BL = HL .* inK ./ (2.0 * primL[end])
    BR = HR .* inK ./ (2.0 * primR[end])

    bc = zeros(4)

    p = (
        x0 = x0,
        x1 = x1,
        wL = wL,
        wR = wR,
        primL = primL,
        primR = primR,
        HL = HL,
        HR = HR,
        BL = BL,
        BR = BR,
    )

    fw = function (x, p)
        if x <= (p.x0 + p.x1) / 2
            return p.wL
        else
            return p.wR
        end
    end
    ff = function (x, p)
        if x <= (p.x0 + p.x1) / 2
            return p.HL, p.BL
        else
            return p.HR, p.BR
        end
    end
    bc = function (x, p)
        if x <= (p.x0 + p.x1) / 2
            return p.primL
        else
            return p.primR
        end
    end

    ib = IB2F(fw, ff, bc, p)

    ks = SolverSet(set, pSpace, vSpace, gas, ib)
end

begin
    ctr = OffsetArray{ControlVolume2F}(undef, axes(ks.pSpace.x, 1))
    face = Array{Interface2F}(undef, ks.pSpace.nx + 1)

    idx0 = (eachindex(pSpace.x)|>collect)[1]
    idx1 = (eachindex(pSpace.x)|>collect)[end]

    for i in eachindex(ctr)
        w = ks.ib.fw(ks.ps.x[i], ks.ib.p)
        prim = conserve_prim(w, γ)
        h, b = ks.ib.ff(ks.ps.x[i], ks.ib.p)
        ctr[i] = ControlVolume(w, prim, h, b, 2)
    end

    for i = 1:ks.pSpace.nx+1
        fw = deepcopy(ks.ib.fw(ks.ps.x[1], ks.ib.p))
        ff = deepcopy(ks.ib.ff(ks.ps.x[1], ks.ib.p)[1])
        face[i] = Interface(fw, ff, ff, 2)
    end
end

begin
    iter = 0
    res = zeros(4)
    simTime = 0.0
    dt = Kinetic.timestep(ks, ctr, simTime)
    maxTime = vhs_collision_time(ks.ib.bc(ks.ps.x0, ks.ib.p), μᵣ, omega)
    nt = Int(floor(maxTime / dt))
end

# There're no default solver for 1D simulation with 2D setting
# Let's do it manually
@showprogress for iter = 1:nt
    #Kinetic.reconstruct!(ks, ctr)

    @inbounds Threads.@threads for i in eachindex(face)
        flux_kfvs!(
            face[i].fw,
            face[i].fh,
            face[i].fb,
            ctr[i-1].h,
            ctr[i-1].b,
            ctr[i].h,
            ctr[i].b,
            ks.vSpace.u,
            ks.vSpace.v,
            ks.vSpace.weights,
            dt,
            1.0,
        )
    end

    @inbounds Threads.@threads for i = 1:ks.pSpace.nx
        #--- store W^n and calculate shakhov term ---#
        w_old = deepcopy(ctr[i].w)

        #--- update W^{n+1} ---#
        @. ctr[i].w += (face[i].fw - face[i+1].fw) / ks.ps.dx[i]
        ctr[i].prim .= conserve_prim(ctr[i].w, ks.gas.γ)

        #--- calculate M^{n+1} and tau^{n+1} ---#
        MH = maxwellian(ks.vSpace.u, ks.vSpace.v, ctr[i].prim)
        MB = MH .* ks.gas.K ./ (2.0 * ctr[i].prim[end])
        τ = vhs_collision_time(ctr[i].prim, ks.gas.μᵣ, ks.gas.ω)

        #--- update distribution function ---#
        for q in axes(MH, 2), p in axes(MH, 1)
            ctr[i].h[p, q] =
                (
                    ctr[i].h[p, q] +
                    (face[i].fh[p, q] - face[i+1].fh[p, q]) / ks.ps.dx[i] +
                    dt / τ * MH[p, q]
                ) / (1.0 + dt / τ)
            ctr[i].b[p, q] =
                (
                    ctr[i].b[p, q] +
                    (face[i].fb[p, q] - face[i+1].fb[p, q]) / ks.ps.dx[i] +
                    dt / τ * MB[p, q]
                ) / (1.0 + dt / τ)
        end
    end
end

sol = zeros(ks.pSpace.nx, 10)
for i = 1:ks.pSpace.nx
    sol[i, 1:3] = ctr[i].prim[1:3]
    sol[i, 4] = 1.0 / ctr[i].prim[4]
end

using Plots
plot(ks.pSpace.x[1:ks.pSpace.nx], sol[:, 1])
plot(ks.pSpace.x[1:ks.pSpace.nx], sol[:, 2])
plot(ks.pSpace.x[1:ks.pSpace.nx], sol[:, 3])
plot(ks.pSpace.x[1:ks.pSpace.nx], sol[:, 4])

u1d = VSpace1D(umin, umax, nu)
f = zeros(ks.vSpace.nv)
for j in axes(ctr[1].h, 2), i in axes(ctr[1].h, 1)
    f[j] =
        0.5 * (
            sum(@. u1d.weights * ctr[ks.pSpace.nx÷2].h[:, j]) +
            sum(@. u1d.weights * ctr[ks.pSpace.nx÷2+1].h[:, j])
        )
end
plot(ks.vSpace.v[end÷2, :, 1], f)
