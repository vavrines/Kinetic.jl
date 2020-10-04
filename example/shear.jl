using Revise, ProgressMeter, OffsetArrays, Kinetic

cd(@__DIR__)
D = read_dict("shear.txt")
for key in keys(D)
    s = Symbol(key)
    @eval $s = $(D[key])
end

begin
    γ = heat_capacity_ratio(inK, 2)
    set = Setup(case, space, flux, collision, nSpecies, interpOrder, limiter, cfl, maxTime)
    pSpace = PSpace1D(x0, x1, nx, pMeshType, nxg)
    vSpace = VSpace2D(umin, umax, nu, vmin, vmax, nv, vMeshType, nug, nvg)
    
    μᵣ = ref_vhs_vis(knudsen, alphaRef, omegaRef)
    gas = Gas(
        knudsen,
        mach,
        prandtl,
        inK,
        γ,
        omega,
        alphaRef,
        omegaRef,
        μᵣ,
    )

    primL = [1.0, 0.0, 1.0, 1.0]
    primR = [1.0, 0.0, -1.0, 2.0]

    wL = prim_conserve(primL, γ)
    wR = prim_conserve(primR, γ)

    HL = maxwellian(vSpace.u, vSpace.v, primL)
    HR = maxwellian(vSpace.u, vSpace.v, primR)

    BL = HL .* inK ./ (2.0 * primL[end])
    BR = HR .* inK ./ (2.0 * primR[end])

    bc = zeros(4)
    
    ib = IB2F(
            wL,
            primL,
            HL,
            BL,
            bc,
            wR,
            primR,
            HR,
            BR,
            bc,
        )

    outputFolder = pwd()

    ks = SolverSet(set, pSpace, vSpace, gas, ib, outputFolder)
    KS = ks
end

begin
    ctr = OffsetArray{ControlVolume1D2F}(undef, axes(KS.pSpace.x, 1))
    face = Array{Interface1D2F}(undef, KS.pSpace.nx + 1)

    idx0 = (eachindex(pSpace.x) |> collect)[1]
    idx1 = (eachindex(pSpace.x) |> collect)[end]

    for i in eachindex(ctr)
        if i <= KS.pSpace.nx ÷ 2                
            ctr[i] = ControlVolume1D2F( KS.pSpace.x[i], KS.pSpace.dx[i], ks.ib.wL, ks.ib.primL, 
            ks.ib.hL, ks.ib.bL )
        else
            ctr[i] = ControlVolume1D2F( KS.pSpace.x[i], KS.pSpace.dx[i], ks.ib.wR, ks.ib.primR, 
            ks.ib.hR, ks.ib.bR )
        end
    end

    face = Array{Interface1D2F}(undef, KS.pSpace.nx+1)
    for i=1:KS.pSpace.nx+1
        face[i] = Interface1D2F(ks.ib.wL, ks.ib.hL)
    end
end

begin
    iter = 0
    res = zeros(4)
    simTime = 0.0
    dt = Kinetic.timestep(KS, ctr, simTime)
    maxTime = vhs_collision_time(ks.ib.primL, μᵣ, omega)
    nt = Int(floor(maxTime / dt))
end

# There're no default solver for 1D simulation with 2D setting
# Let's do it manually
@showprogress for iter in 1:nt
    #Kinetic.reconstruct!(KS, ctr)

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
            1.0
        )
    end

    @inbounds Threads.@threads for i in 1:ks.pSpace.nx
        #--- store W^n and calculate shakhov term ---#
        w_old = deepcopy(ctr[i].w)

        #--- update W^{n+1} ---#
        @. ctr[i].w += (face[i].fw - face[i+1].fw) / ctr[i].dx
        ctr[i].prim .= conserve_prim(ctr[i].w, ks.gas.γ)

        #--- calculate M^{n+1} and tau^{n+1} ---#
        MH = maxwellian(ks.vSpace.u, ks.vSpace.v, ctr[i].prim)
        MB = MH .* ks.gas.K ./ (2.0 * ctr[i].prim[end])
        τ = vhs_collision_time(ctr[i].prim, ks.gas.μᵣ, ks.gas.ω)

        #--- update distribution function ---#
        for q in axes(MH,2), p in axes(MH,1)
            ctr[i].h[p,q] = (ctr[i].h[p,q] + (face[i].fh[p,q] - face[i+1].fh[p,q]) / ctr[i].dx + dt / τ * MH[p,q]) / (1.0 + dt / τ)
            ctr[i].b[p,q] = (ctr[i].b[p,q] + (face[i].fb[p,q] - face[i+1].fb[p,q]) / ctr[i].dx + dt / τ * MB[p,q]) / (1.0 + dt / τ)
        end
    end
end

sol = zeros(ks.pSpace.nx, 10)
for i in 1:ks.pSpace.nx
    sol[i, 1:3] = ctr[i].prim[1:3]
    sol[i, 4] = 1. / ctr[i].prim[4]
end

using Plots
plot(ks.pSpace.x[1:ks.pSpace.nx], sol[:,1])
plot(ks.pSpace.x[1:ks.pSpace.nx], sol[:,2])
plot(ks.pSpace.x[1:ks.pSpace.nx], sol[:,3])
plot(ks.pSpace.x[1:ks.pSpace.nx], sol[:,4])

u1d = VSpace1D(umin, umax, nu, vMeshType)
f = zeros(ks.vSpace.nv)
for j in axes(ctr[1].h, 2), i in axes(ctr[1].h, 1)
    f[j] = 0.5 * (sum(@. u1d.weights * ctr[ks.pSpace.nx÷2].h[:,j]) + 
        sum(@. u1d.weights * ctr[ks.pSpace.nx÷2+1].h[:,j]))
end
plot(ks.vSpace.v[end÷2,:,1], f)