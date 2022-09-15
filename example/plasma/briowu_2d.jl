using KitBase, KitBase.OffsetArrays
using KitBase.ProgressMeter: @showprogress

cd(@__DIR__)
D = read_dict("briowu_2d.txt")
for key in keys(D)
    s = Symbol(key)
    @eval $s = $(D[key])
end

begin
    γ = heat_capacity_ratio(D[:inK], 3)
    set = set_setup(D)
    ps = set_geometry(D)

    ue0 = umin * sqrt(mi / me)
    ue1 = umax * sqrt(mi / me)
    ve0 = vmin * sqrt(mi / me)
    ve1 = vmax * sqrt(mi / me)
    kne = knudsen * (me / mi)

    vs = MVSpace2D(umin, umax, ue0, ue1, nu, vmin, vmax, ve0, ve1, nv; type = vMeshType, ngu = nug, ngv = nvg)
    plasma = Plasma1D([knudsen,kne], mach, prandtl, inK, γ, mi, ni, me, ne, lD, rL, sol, echi, bnu)

    begin
        # upstream
        primL = zeros(5, 2)
        primL[1, 1] = 1.0 * mi
        primL[2, 1] = 0.0
        primL[3, 1] = 0.0
        primL[4, 1] = 0.0
        primL[5, 1] = mi / 1.0
        primL[1, 2] = 1.0 * me
        primL[2, 2] = 0.0
        primL[3, 2] = 0.0
        primL[4, 2] = 0.0
        primL[5, 2] = me / 1.0

        wL = mixture_prim_conserve(primL, γ)
        h0L = mixture_maxwellian(vs.u, vs.v, primL)

        h1L = similar(h0L)
        h2L = similar(h0L)
        for j in axes(h0L, 3)
            h1L[:, :, j] .= primL[4, j] .* h0L[:, :, j]
            h2L[:, :, j] .= (primL[4, j]^2 + 1.0 / (2.0 * primL[end, j])) .* h0L[:, :, j]
        end

        EL = zeros(3)
        BL = zeros(3)
        BL[1] = 0.75
        BL[2] = 1.0

        # downstream
        primR = zeros(5, 2)
        primR[1, 1] = 0.125 * mi
        primR[2, 1] = 0.0
        primR[3, 1] = 0.0
        primR[4, 1] = 0.0
        primR[5, 1] = mi * 1.25
        primR[1, 2] = 0.125 * me
        primR[2, 2] = 0.0
        primR[3, 2] = 0.0
        primR[4, 2] = 0.0
        primR[5, 2] = me * 1.25

        wR = mixture_prim_conserve(primR, γ)
        h0R = mixture_maxwellian(vs.u, vs.v, primR)

        h1R = similar(h0R)
        h2R = similar(h0R)
        for j in axes(h0R, 3)
            h1R[:, :, j] .= primR[4, j] .* h0R[:, :, j]
            h2R[:, :, j] .= (primR[4, j]^2 + 1.0 / (2.0 * primR[end, j])) .* h0R[:, :, j]
        end

        ER = zeros(3)
        BR = zeros(3)
        BR[1] = 0.75
        BR[2] = -1.0

        lorenzL = zeros(3, 2)
        lorenzR = zeros(3, 2)
        bcL = zeros(5, 2)
        bcR = zeros(5, 2)

        p = (
            x0 = x0,
            x1 = x1,
            wL = wL,
            wR = wR,
            primL = primL,
            primR = primR,
            h0L = h0L,
            h1L = h1L,
            h2L = h2L,
            h0R = h0R,
            h1R = h1R,
            h2R = h2R,
            EL = EL,
            ER = ER,
            BL = BL,
            BR = BR,
            lorenzL = lorenzL,
            lorenzR = lorenzR,
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
                return p.h0L, p.h1L, p.h2L
            else
                return p.h0R, p.h1R, p.h2R
            end
        end
        fE = function (x, p)
            if x <= (p.x0 + p.x1) / 2
                return p.EL
            else
                return p.ER
            end
        end
        fB = function (x, p)
            if x <= (p.x0 + p.x1) / 2
                return p.BL
            else
                return p.BR
            end
        end
        fL = function (x, p)
            if x <= (p.x0 + p.x1) / 2
                return p.lorenzL
            else
                return p.lorenzR
            end
        end
        bc = function (x, p)
            if x <= (p.x0 + p.x1) / 2
                return p.primL
            else
                return p.primR
            end
        end
        ib = IB3F(fw, ff, fE, fB, fL, bc, p)
    end
    ks = SolverSet(set, ps, vs, plasma, ib)
end

ctr, face = init_fvm(ks)

begin
    t = 0.
    dt = timestep(ks, ctr, t)
    nt = Int(floor(ks.set.maxTime / dt))+1
    res = zeros(5, 2)
end

@showprogress for iter in 1:nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt; mode=:kcu, isPlasma=:true)
    update!(ks, ctr, face, dt, res, isMHD=true)
end

soluiton = zeros(ks.ps.nx, 10, 2)
for i in 1:ks.ps.nx
    soluiton[i, 1, 1] = ctr[i].prim[1,1]
    soluiton[i, 1, 2] = ctr[i].prim[1,2] / ks.gas.me
    soluiton[i, 2:4, 1] .= ctr[i].prim[2:4,1]
    soluiton[i, 2:4, 2] .= ctr[i].prim[2:4,2]
    soluiton[i, 5, 1] = 1. / ctr[i].prim[5,1]
    soluiton[i, 5, 2] = ks.gas.me / ctr[i].prim[5,2]

    soluiton[i, 6, 1] = ctr[i].B[2]
    soluiton[i, 6, 2] = ctr[i].E[1]
end

using Plots
plot(ks.ps.x[1:ks.ps.nx], soluiton[:,1,1:2])
plot(ks.ps.x[1:ks.ps.nx], soluiton[:,6,1:2])
