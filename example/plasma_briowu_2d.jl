using Revise, ProgressMeter, OffsetArrays, Kinetic

cd(@__DIR__)
D = read_dict("briowu_2d.txt")
for key in keys(D)
    s = Symbol(key)
    @eval $s = $(D[key])
end

begin
    γ = heat_capacity_ratio(inK, 3)
    set = Setup(case, space, nSpecies, interpOrder, limiter, cfl, maxTime)
    pSpace = PSpace1D(x0, x1, nx, pMeshType, nxg)

    ue0 = umin * sqrt(mi / me)
    ue1 = umax * sqrt(mi / me)
    ve0 = vmin * sqrt(mi / me)
    ve1 = vmax * sqrt(mi / me)
    kne = knudsen * (me / mi)

    vSpace = MVSpace2D(umin, umax, ue0, ue1, nu, vmin, vmax, ve0, ve1, nv, vMeshType, nug, nvg)
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
        h0L = mixture_maxwellian(vSpace.u, vSpace.v, primL)

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
        h0R = mixture_maxwellian(vSpace.u, vSpace.v, primR)

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
    end

    ib = IB3F(
            wL,
            primL,
            h0L,
            h1L,
            h2L,
            bcL,
            EL,
            BL,
            lorenzL,
            wR,
            primR,
            h0R,
            h1R,
            h2R,
            bcR,
            ER,
            BR,
            lorenzR,
        )

    outputFolder = pwd()

    ks = SolverSet(set, pSpace, vSpace, plasma, ib, outputFolder)
    KS = ks

    ctr = OffsetArray{ControlVolume1D3F}(undef, axes(KS.pSpace.x, 1))
    face = Array{Interface1D3F}(undef, KS.pSpace.nx + 1)
end

begin
    idx0 = (eachindex(pSpace.x) |> collect)[1]
    idx1 = (eachindex(pSpace.x) |> collect)[end]

    for i in eachindex(ctr)
        if i <= KS.pSpace.nx ÷ 2                
            ctr[i] = ControlVolume1D3F( KS.pSpace.x[i], KS.pSpace.dx[i], KS.ib.wL, KS.ib.primL, 
            KS.ib.h0L, KS.ib.h1L, KS.ib.h2L, KS.ib.EL, KS.ib.BL, KS.ib.lorenzL )
        else
            ctr[i] = ControlVolume1D3F( KS.pSpace.x[i], KS.pSpace.dx[i], KS.ib.wR, KS.ib.primR, 
            KS.ib.h0R, KS.ib.h1R, KS.ib.h2R, KS.ib.ER, KS.ib.BR, KS.ib.lorenzR )
        end
    end

    face = Array{Interface1D3F}(undef, KS.pSpace.nx+1)
    for i=1:KS.pSpace.nx+1
        face[i] = Interface1D3F(KS.ib.wL, KS.ib.h0L, KS.ib.EL)
    end
end

begin
    simTime = 0.
    dt = Kinetic.timestep(KS, ctr, simTime)
    nt = Int(floor(ks.set.maxTime / dt))+1
    res = zeros(5, 2)
end

@showprogress for iter in 1:nt
    #dt = Kinetic.timestep(KS, ctr, simTime)
    #Kinetic.reconstruct!(ks, ctr)

    Kinetic.evolve!(ks, ctr, face, dt; mode=:kfvs, isPlasma=:true)
    # it's equivalent to the following process
    #=
    @inbounds Threads.@threads for i = 1:KS.pSpace.nx+1
        flux_kfvs!(
            face[i].fw,
            face[i].fh0,
            face[i].fh1,
            face[i].fh2,
            ctr[i-1].h0,
            ctr[i-1].h1,
            ctr[i-1].h2,
            ctr[i].h0,
            ctr[i].h1,
            ctr[i].h2,
            KS.vSpace.u,
            KS.vSpace.v,
            KS.vSpace.weights,
            dt,
            1.0,
        )

        flux_em!(
            face[i].femL,
            face[i].femR,
            ctr[i-2].E,
            ctr[i-2].B,
            ctr[i-1].E,
            ctr[i-1].B,
            ctr[i].E,
            ctr[i].B,
            ctr[i+1].E,
            ctr[i+1].B,
            ctr[i-1].ϕ,
            ctr[i].ϕ,
            ctr[i-1].ψ,
            ctr[i].ψ,
            ctr[i-1].dx,
            ctr[i].dx,
            KS.gas.A1p,
            KS.gas.A1n,
            KS.gas.D1,
            KS.gas.sol,
            KS.gas.χ,
            KS.gas.ν,
            dt,
        )
    end
    =#

    Kinetic.update!(ks, ctr, face, dt, res, isMHD=true)
end

soluiton = zeros(ks.pSpace.nx, 10, 2)
for i in 1:ks.pSpace.nx
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
plot(ks.pSpace.x[1:ks.pSpace.nx], soluiton[:,1,1:2])
plot(ks.pSpace.x[1:ks.pSpace.nx], soluiton[:,6,1:2])

using JLD2
@save "refsol.jld2" ks ctr