using ProgressMeter, OffsetArrays, Kinetic

cd(@__DIR__)
D = read_dict("mixture_shock.txt")
for key in keys(D)
    s = Symbol(key)
    @eval $s = $(D[key])
end

begin
    γ = heat_capacity_ratio(inK, 1)
    set = Setup(case, space, flux, collision, nSpecies, interpOrder, limiter, cfl, maxTime)
    pSpace = PSpace1D(x0, x1, nx, pMeshType, nxg)
    ue0 = umin * sqrt(mi / me)
    ue1 = umax * sqrt(mi / me)
    vSpace = MVSpace1D(umin, umax, ue0, ue1, nu, vMeshType, nug)
    kne = knudsen * (me / mi)
    gas = Mixture(
        [knudsen, kne],
        mach,
        prandtl,
        inK,
        γ,
        mi,
        ni,
        me,
        ne,
    )
    wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR = ib_rh(gas.Ma, gas.γ, gas.K, mi, me, ni, ne, vSpace.u)
    ib = IB2F(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR)

    ks = SolverSet(set, pSpace, vSpace, gas, ib, pwd())
end

begin
    ctr = OffsetArray{ControlVolume1D2F}(undef, axes(ks.pSpace.x, 1))
    face = Array{Interface1D2F}(undef, ks.pSpace.nx + 1)

    idx0 = (eachindex(pSpace.x) |> collect)[1]
    idx1 = (eachindex(pSpace.x) |> collect)[end]

    for i in eachindex(ctr)
        if i <= ks.pSpace.nx ÷ 2                
            ctr[i] = ControlVolume1D2F(ks.pSpace.x[i], ks.pSpace.dx[i], ks.ib.wL, ks.ib.primL, 
            ks.ib.hL, ks.ib.bL)
        else
            ctr[i] = ControlVolume1D2F(ks.pSpace.x[i], ks.pSpace.dx[i], ks.ib.wR, ks.ib.primR, 
            ks.ib.hR, ks.ib.bR)
        end
    end

    face = Array{Interface1D2F}(undef, ks.pSpace.nx+1)
    for i=1:ks.pSpace.nx+1
        face[i] = Interface1D2F(ks.ib.wL, ks.ib.hL)
    end
end

begin
    iter = 0
    res = zeros(3)
    simTime = 0.0
    dt = Kinetic.timestep(ks, ctr, simTime)
    nt = Int(floor(ks.set.maxTime / dt))
end

Kinetic.reconstruct!(ks, ctr)

Kinetic.evolve!(ks, ctr, face, dt)

res = zeros(3, 2)
Kinetic.update!(ks, ctr, face, dt, res; bc=:fix)
