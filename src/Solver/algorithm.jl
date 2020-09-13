# ============================================================
# Solution Algorithm
# ============================================================

"""
Structure of solver setup

"""
struct SolverSet <: AbstractSolverSet

    # setup
    set::AbstractSetup

    # physical space
    pSpace::AbstractPhysicalSpace

    # velocity space
    vSpace::AbstractVelocitySpace

    # gas property
    gas::AbstractProperty

    # initial and boundary condition
    ib::AbstractCondition

    # file system
    outputFolder::String

    # constructor
    #SolverSet() = SolverSet("./config/config.txt")

    SolverSet(set, pSpace, vSpace, gas, ib, outputFolder) =
        new(set, pSpace, vSpace, gas, ib, outputFolder)

    function SolverSet(configfilename::String)

        # read settings from configuration file
        D = read_dict(configfilename)

        # generate variables from dictionary
        for key in keys(D)
            s = Symbol(key)
            @eval $s = $(D[key])
        end

        # generate data structure
        dim = parse(Int, space[1])
        gasD = map(parse(Int, space[3]), parse(Int, space[5])) do x, y # (x)f(y)v
            if x == 1
                if y >= 3
                    return 3
                else
                    return dim
                end
            elseif x == 2
                return dim
            elseif x >= 3
                return 3
            else
                return nothing
            end
        end
        γ = heat_capacity_ratio(inK, gasD)
        set = Setup(case, space, nSpecies, interpOrder, limiter, cfl, maxTime)

        if dim == 1

            pSpace = PSpace1D(x0, x1, nx, pMeshType, nxg)

            if case == "shock"

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

                if space == "1d1f1v"
                    vSpace = VSpace1D(umin, umax, nu, vMeshType, nug)

                    wL, primL, fL, bcL, wR, primR, fR, bcR = ib_rh(mach, γ, vSpace.u)

                    ib = IB1F(wL, primL, fL, bcL, wR, primR, fR, bcR)
                elseif space == "1d2f1v"
                    vSpace = VSpace1D(umin, umax, nu, vMeshType, nug)

                    wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR =
                        ib_rh(mach, γ, vSpace.u, inK)

                    ib = IB2F(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR)
                elseif space == "1d1f3v"
                    vSpace = VSpace3D(
                        umin,
                        umax,
                        nu,
                        vmin,
                        vmax,
                        nv,
                        wmin,
                        wmax,
                        nw,
                        vMeshType,
                        nug,
                        nvg,
                        nwg,
                    )

                    wL, primL, fL, bcL, wR, primR, fR, bcR =
                        ib_rh(mach, γ, vSpace.u, vSpace.v, vSpace.w)

                    ib = IB1F(wL, primL, fL, bcL, wR, primR, fR, bcR)
                end

            elseif case == "sod"

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

                if space == "1d1f1v"
                    vSpace = VSpace1D(umin, umax, nu, vMeshType, nug)

                    wL, primL, fL, bcL, wR, primR, fR, bcR = ib_sod(γ, vSpace.u)

                    ib = IB1F(wL, primL, fL, bcL, wR, primR, fR, bcR)
                elseif space == "1d2f1v"
                    vSpace = VSpace1D(umin, umax, nu, vMeshType, nug)

                    wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR =
                        ib_sod(γ, vSpace.u, inK)

                    ib = IB2F(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR)
                elseif space == "1d1f3v"
                    vSpace = VSpace3D(
                        umin,
                        umax,
                        nu,
                        vmin,
                        vmax,
                        nv,
                        wmin,
                        wmax,
                        nw,
                        vMeshType,
                        nug,
                        nvg,
                        nwg,
                    )

                    wL, primL, fL, bcL, wR, primR, fR, bcR =
                        ib_sod(γ, vSpace.u, vSpace.v, vSpace.w)

                    ib = IB1F(wL, primL, fL, bcL, wR, primR, fR, bcR)
                end

            elseif case == "brio-wu"

                v0 = umin * sqrt(mi / me)
                v1 = umax * sqrt(mi / me)
                kne = knudsen * (me / mi)

                vSpace = MVSpace1D(umin, umax, v0, v1, nu, vMeshType, nug)

                gas = Plasma1D(
                    [knudsen, kne],
                    mach,
                    prandtl,
                    inK,
                    γ,
                    mi,
                    ni,
                    me,
                    ne,
                    lD,
                    rL,
                    sol,
                    echi,
                    bnu,
                )

                wL,
                primL,
                h0L,
                h1L,
                h2L,
                h3L,
                bcL,
                EL,
                BL,
                lorenzL,
                wR,
                primR,
                h0R,
                h1R,
                h2R,
                h3R,
                bcR,
                ER,
                BR,
                lorenzR = ib_briowu(γ, vSpace.u, mi, me)

                ib = IB4F(
                    wL,
                    primL,
                    h0L,
                    h1L,
                    h2L,
                    h3L,
                    bcL,
                    EL,
                    BL,
                    lorenzL,
                    wR,
                    primR,
                    h0R,
                    h1R,
                    h2R,
                    h3R,
                    bcR,
                    ER,
                    BR,
                    lorenzR,
                )

            else

                throw("No preset available, please set up manually.")

            end

        elseif dim == 2

            pSpace = PSpace2D(x0, x1, nx, y0, y1, ny, pMeshType, nxg, nyg)

            if case == "cavity"

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

                if space == "2d1f2v"
                    vSpace = VSpace2D(umin, umax, nu, vmin, vmax, nv, vMeshType, nug, nvg)

                    wL, primL, fL, bcL, wR, primR, fR, bcR, bcU, bcD =
                        ib_cavity(γ, uLid, vLid, tLid, vSpace.u, vSpace.v)

                    ib = IB1F(wL, primL, fL, bcL, wR, primR, fR, bcR, bcU, bcD)
                elseif space == "2d2f2v"
                    vSpace = VSpace2D(umin, umax, nu, vmin, vmax, nv, vMeshType, nug, nvg)

                    wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR, bcU, bcD =
                        ib_cavity(γ, uLid, vLid, tLid, vSpace.u, vSpace.v, inK)

                    ib = IB2F(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR, bcU, bcD)
                end

            else

                throw("No preset available, please set up manually.")

            end

        else

            throw("No preset available for 3D simulation, please set up manually.")

        end

        # create working directory
        identifier = string(Dates.now(), "/")
        #outputFolder = string("../out/", replace(identifier, ":"=>"."))
        outputFolder = replace(identifier, ":" => ".")
        mkdir(outputFolder)
        mkdir(string(outputFolder, "data/"))
        cp(configfilename, string(outputFolder, "config.txt"))

        # create new struct
        new(set, pSpace, vSpace, gas, ib, outputFolder)

    end # function

end # struct


"""
Solution algorithm

"""
function solve!(
    KS::SolverSet,
    ctr::AbstractArray{<:AbstractControlVolume1D,1},
    face::Array{<:AbstractInterface1D,1},
    simTime::Float64,
)

    #--- setup ---#
    iter = 0
    dt = timestep(KS, ctr, simTime)
    nt = Int(Floor(KS.set.maxTime / dt)) + 1
    res = zeros(axes(KS.ib.wL))
    write_jld(KS, ctr, simTime)

    #--- main loop ---#
    #while true
    @showprogress for iter = 1:nt

        #dt = timestep(KS, ctr, simTime)
        reconstruct!(KS, ctr)
        evolve!(KS, ctr, face, dt)
        update!(KS, ctr, face, dt, res)

        #iter += 1
        simTime += dt

        if iter % 100 == 0
            println("iter: $(iter), time: $(simTime), dt: $(dt), res: $(res[1:end])")

            #if iter%1000 == 0
            #write_jld(KS, ctr, iter)
            #end
        end

        if simTime > KS.set.maxTime || maximum(res) < 5.e-7
            break
        end

    end # loop

    write_jld(KS, ctr, simTime)

end # function


# ------------------------------------------------------------
# Calculate time step
# ------------------------------------------------------------
function timestep(
    KS::SolverSet,
    ctr::AbstractArray{<:AbstractControlVolume1D,1},
    simTime::Real,
)

    tmax = 0.0

    if KS.set.nSpecies == 1

        @inbounds Threads.@threads for i = 1:KS.pSpace.nx
            prim = ctr[i].prim
            sos = sound_speed(prim, KS.gas.γ)
            vmax = max(KS.vSpace.u1, abs(prim[2])) + sos
            tmax = max(tmax, vmax / ctr[i].dx)
        end

    elseif KS.set.nSpecies == 2

        @inbounds Threads.@threads for i = 1:KS.pSpace.nx
            prim = ctr[i].prim
            sos = sound_speed(prim, KS.gas.γ)
            vmax = max(maximum(KS.vSpace.u1), maximum(abs.(prim[2, :]))) + sos
            tmax = ifelse(
                KS.set.space == "1d4f",
                max(tmax, vmax / ctr[i].dx, KS.gas.sol / ctr[i].dx),
                max(tmax, vmax / ctr[i].dx),
            )
        end

    end

    dt = KS.set.cfl / tmax
    dt = ifelse(dt < (KS.set.maxTime - simTime), dt, KS.set.maxTime - simTime)

    return dt

end


# ------------------------------------------------------------
# Reconstruction
# ------------------------------------------------------------
function reconstruct!(KS::SolverSet, ctr::AbstractArray{<:AbstractControlVolume1D,1})

    if KS.set.interpOrder == 1
        return
    end

    # macroscopic variables
    @inbounds Threads.@threads for i = 2:KS.pSpace.nx-1
        reconstruct3!(
            ctr[i].sw,
            ctr[i-1].w,
            ctr[i].w,
            ctr[i+1].w,
            0.5 * (ctr[i-1].dx + ctr[i].dx),
            0.5 * (ctr[i].dx + ctr[i+1].dx),
            Symbol(KS.set.limiter),
        )
    end

    # particle distribution function
    @inbounds Threads.@threads for i = 2:KS.pSpace.nx-1
        if KS.set.space[3:4] == "1f"
            reconstruct3!(
                ctr[i].sf,
                ctr[i-1].f,
                ctr[i].f,
                ctr[i+1].f,
                0.5 * (ctr[i-1].dx + ctr[i].dx),
                0.5 * (ctr[i].dx + ctr[i+1].dx),
                Symbol(KS.set.limiter),
            )
        elseif KS.set.space[3:4] == "2f"
            reconstruct3!(
                ctr[i].sh,
                ctr[i-1].h,
                ctr[i].h,
                ctr[i+1].h,
                0.5 * (ctr[i-1].dx + ctr[i].dx),
                0.5 * (ctr[i].dx + ctr[i+1].dx),
                Symbol(KS.set.limiter),
            )
            reconstruct3!(
                ctr[i].sb,
                ctr[i-1].b,
                ctr[i].b,
                ctr[i+1].b,
                0.5 * (ctr[i-1].dx + ctr[i].dx),
                0.5 * (ctr[i].dx + ctr[i+1].dx),
                Symbol(KS.set.limiter),
            )
        elseif KS.set.space[3:4] == "4f"
            reconstruct3!(
                ctr[i].sh0,
                ctr[i-1].h0,
                ctr[i].h0,
                ctr[i+1].h0,
                0.5 * (ctr[i-1].dx + ctr[i].dx),
                0.5 * (ctr[i].dx + ctr[i+1].dx),
                Symbol(KS.set.limiter),
            )
            reconstruct3!(
                ctr[i].sh1,
                ctr[i-1].h1,
                ctr[i].h1,
                ctr[i+1].h1,
                0.5 * (ctr[i-1].dx + ctr[i].dx),
                0.5 * (ctr[i].dx + ctr[i+1].dx),
                Symbol(KS.set.limiter),
            )
            reconstruct3!(
                ctr[i].sh2,
                ctr[i-1].h2,
                ctr[i].h2,
                ctr[i+1].h2,
                0.5 * (ctr[i-1].dx + ctr[i].dx),
                0.5 * (ctr[i].dx + ctr[i+1].dx),
                Symbol(KS.set.limiter),
            )
            reconstruct3!(
                ctr[i].sh3,
                ctr[i-1].h3,
                ctr[i].h3,
                ctr[i+1].h3,
                0.5 * (ctr[i-1].dx + ctr[i].dx),
                0.5 * (ctr[i].dx + ctr[i+1].dx),
                Symbol(KS.set.limiter),
            )
        elseif KS.set.space[3:4] == "3f"
            reconstruct3!(
                ctr[i].sh0,
                ctr[i-1].h0,
                ctr[i].h0,
                ctr[i+1].h0,
                0.5 * (ctr[i-1].dx + ctr[i].dx),
                0.5 * (ctr[i].dx + ctr[i+1].dx),
                Symbol(KS.set.limiter),
            )
            reconstruct3!(
                ctr[i].sh1,
                ctr[i-1].h1,
                ctr[i].h1,
                ctr[i+1].h1,
                0.5 * (ctr[i-1].dx + ctr[i].dx),
                0.5 * (ctr[i].dx + ctr[i+1].dx),
                Symbol(KS.set.limiter),
            )
            reconstruct3!(
                ctr[i].sh2,
                ctr[i-1].h2,
                ctr[i].h2,
                ctr[i+1].h2,
                0.5 * (ctr[i-1].dx + ctr[i].dx),
                0.5 * (ctr[i].dx + ctr[i+1].dx),
                Symbol(KS.set.limiter),
            )
        end
    end

end

function reconstruct!(KS::SolverSet, ctr::AbstractArray{ControlVolume2D2F,2})

    if KS.set.interpOrder == 1
        return
    end

    #--- macroscopic variables ---#
    # boundary
    @inbounds Threads.@threads for j in 1:KS.pSpace.ny
        swL = extract_last(ctr[1, j].sw, 1, mode=:view)
        reconstruct2!(
            swL,
            ctr[1, j].w,
            ctr[2, j].w,
            0.5 * (ctr[1, j].dx + ctr[2, j].dx),
        )

        swR = extract_last(ctr[KS.pSpace.nx, j].sw, 1, mode=:view)
        reconstruct2!(
            swR,
            ctr[KS.pSpace.nx-1, j].w,
            ctr[KS.pSpace.nx, j].w,
            0.5 * (ctr[KS.pSpace.nx-1, j].dx + ctr[KS.pSpace.nx, j].dx),
        )
    end

    @inbounds Threads.@threads for i in 1:KS.pSpace.nx
        swD = extract_last(ctr[i, 1].sw, 2, mode=:view)
        reconstruct2!(
            swD,
            ctr[i, 1].w,
            ctr[i, 2].w,
            0.5 * (ctr[i, 1].dy + ctr[i, 2].dy),
        )

        swU = extract_last(ctr[i, KS.pSpace.ny].sw, 2, mode=:view)
        reconstruct2!(
            swU,
            ctr[i, KS.pSpace.ny-1].w,
            ctr[i, KS.pSpace.ny].w,
            0.5 * (ctr[i, KS.pSpace.ny-1].dx + ctr[i, KS.pSpace.ny].dx),
        )
    end

    # inner
    @inbounds Threads.@threads for j in 1:KS.pSpace.ny
        for i in 2:KS.pSpace.nx-1
            sw = extract_last(ctr[i, j].sw, 1, mode=:view)
            reconstruct3!(
                sw,
                ctr[i-1, j].w,
                ctr[i, j].w,
                ctr[i+1, j].w,
                0.5 * (ctr[i-1, j].dx + ctr[i, j].dx),
                0.5 * (ctr[i, j].dx + ctr[i+1, j].dx),
                Symbol(KS.set.limiter),
            )
        end
    end

    @inbounds Threads.@threads for j in 2:KS.pSpace.ny-1
        for i in 1:KS.pSpace.nx
            sw = extract_last(ctr[i, j].sw, 2, mode=:view)
            reconstruct3!(
                sw,
                ctr[i, j-1].w,
                ctr[i, j].w,
                ctr[i, j+1].w,
                0.5 * (ctr[i, j-1].dy + ctr[i, j].dy),
                0.5 * (ctr[i, j].dy + ctr[i, j+1].dy),
                Symbol(KS.set.limiter),
            )
        end
    end

    #--- particle distribution function ---#
    # boundary
    @inbounds Threads.@threads for j in 1:KS.pSpace.ny
        shL = extract_last(ctr[1, j].sh, 1, mode=:view)
        reconstruct2!(
            shL,
            ctr[1, j].h,
            ctr[2, j].h,
            0.5 * (ctr[1, j].dx + ctr[2, j].dx),
        )
        sbL = extract_last(ctr[1, j].sb, 1, mode=:view)
        reconstruct2!(
            sbL,
            ctr[1, j].b,
            ctr[2, j].b,
            0.5 * (ctr[1, j].dx + ctr[2, j].dx),
        )

        shR = extract_last(ctr[KS.pSpace.nx, j].sh, 1, mode=:view)
        reconstruct2!(
            shR,
            ctr[KS.pSpace.nx-1, j].h,
            ctr[KS.pSpace.nx, j].h,
            0.5 * (ctr[KS.pSpace.nx-1, j].dx + ctr[KS.pSpace.nx, j].dx),
        )
        sbR = extract_last(ctr[KS.pSpace.nx, j].sb, 1, mode=:view)
        reconstruct2!(
            sbR,
            ctr[KS.pSpace.nx-1, j].b,
            ctr[KS.pSpace.nx, j].b,
            0.5 * (ctr[KS.pSpace.nx-1, j].dx + ctr[KS.pSpace.nx, j].dx),
        )
    end

    @inbounds Threads.@threads for i in 1:KS.pSpace.nx
        shD = extract_last(ctr[i, 1].sh, 2, mode=:view)
        reconstruct2!(
            shD,
            ctr[i, 1].h,
            ctr[i, 2].h,
            0.5 * (ctr[i, 1].dy + ctr[i, 2].dy),
        )
        sbD = extract_last(ctr[i, 1].sb, 2, mode=:view)
        reconstruct2!(
            sbD,
            ctr[i, 1].b,
            ctr[i, 2].b,
            0.5 * (ctr[i, 1].dy + ctr[i, 2].dy),
        )

        shU = extract_last(ctr[i, KS.pSpace.ny].sh, 2, mode=:view)
        reconstruct2!(
            shU,
            ctr[i, KS.pSpace.ny-1].h,
            ctr[i, KS.pSpace.ny].h,
            0.5 * (ctr[i, KS.pSpace.ny-1].dx + ctr[i, KS.pSpace.ny].dx),
        )
        sbU = extract_last(ctr[i, KS.pSpace.ny].sb, 2, mode=:view)
        reconstruct2!(
            sbU,
            ctr[i, KS.pSpace.ny-1].b,
            ctr[i, KS.pSpace.ny].b,
            0.5 * (ctr[i, KS.pSpace.ny-1].dx + ctr[i, KS.pSpace.ny].dx),
        )
    end

    # inner
    @inbounds Threads.@threads for j in 1:KS.pSpace.ny
        for i in 2:KS.pSpace.nx-1
            sh = extract_last(ctr[i, j].sh, 1, mode=:view)
            reconstruct3!(
                sh,
                ctr[i-1, j].h,
                ctr[i, j].h,
                ctr[i+1, j].h,
                0.5 * (ctr[i-1, j].dx + ctr[i, j].dx),
                0.5 * (ctr[i, j].dx + ctr[i+1, j].dx),
                Symbol(KS.set.limiter),
            )

            sb = extract_last(ctr[i, j].sb, 1, mode=:view)
            reconstruct3!(
                sb,
                ctr[i-1, j].b,
                ctr[i, j].b,
                ctr[i+1, j].b,
                0.5 * (ctr[i-1, j].dx + ctr[i, j].dx),
                0.5 * (ctr[i, j].dx + ctr[i+1, j].dx),
                Symbol(KS.set.limiter),
            )
        end
    end

    @inbounds Threads.@threads for j in 2:KS.pSpace.ny-1
        for i in 1:KS.pSpace.nx
            sh = extract_last(ctr[i, j].sh, 2, mode=:view)
            reconstruct3!(
                sh,
                ctr[i, j-1].h,
                ctr[i, j].h,
                ctr[i, j+1].h,
                0.5 * (ctr[i, j-1].dy + ctr[i, j].dy),
                0.5 * (ctr[i, j].dy + ctr[i, j+1].dy),
                Symbol(KS.set.limiter),
            )

            sb = extract_last(ctr[i, j].sb, 2, mode=:view)
            reconstruct3!(
                sb,
                ctr[i, j-1].b,
                ctr[i, j].b,
                ctr[i, j+1].b,
                0.5 * (ctr[i, j-1].dy + ctr[i, j].dy),
                0.5 * (ctr[i, j].dy + ctr[i, j+1].dy),
                Symbol(KS.set.limiter),
            )
        end
    end

end

function reconstruct!(KS::SolverSet, ctr::AbstractArray{ControlVolume2D3F,2})

    if KS.set.interpOrder == 1
        return
    end

    #--- macroscopic variables ---#
    # boundary
    @inbounds Threads.@threads for j in 1:KS.pSpace.ny
        swL = extract_last(ctr[1, j].sw, 1, mode=:view)
        reconstruct2!(
            swL,
            ctr[1, j].w,
            ctr[2, j].w,
            0.5 * (ctr[1, j].dx + ctr[2, j].dx),
        )

        swR = extract_last(ctr[KS.pSpace.nx, j].sw, 1, mode=:view)
        reconstruct2!(
            swR,
            ctr[KS.pSpace.nx-1, j].w,
            ctr[KS.pSpace.nx, j].w,
            0.5 * (ctr[KS.pSpace.nx-1, j].dx + ctr[KS.pSpace.nx, j].dx),
        )
    end

    @inbounds Threads.@threads for i in 1:KS.pSpace.nx
        swD = extract_last(ctr[i, 1].sw, 2, mode=:view)
        reconstruct2!(
            swD,
            ctr[i, 1].w,
            ctr[i, 2].w,
            0.5 * (ctr[i, 1].dy + ctr[i, 2].dy),
        )

        swU = extract_last(ctr[i, KS.pSpace.ny].sw, 2, mode=:view)
        reconstruct2!(
            swU,
            ctr[i, KS.pSpace.ny-1].w,
            ctr[i, KS.pSpace.ny].w,
            0.5 * (ctr[i, KS.pSpace.ny-1].dx + ctr[i, KS.pSpace.ny].dx),
        )
    end

    # inner
    @inbounds Threads.@threads for j in 1:KS.pSpace.ny
        for i in 2:KS.pSpace.nx-1
            sw = extract_last(ctr[i, j].sw, 1, mode=:view)
            reconstruct3!(
                sw,
                ctr[i-1, j].w,
                ctr[i, j].w,
                ctr[i+1, j].w,
                0.5 * (ctr[i-1, j].dx + ctr[i, j].dx),
                0.5 * (ctr[i, j].dx + ctr[i+1, j].dx),
                Symbol(KS.set.limiter),
            )
        end
    end

    @inbounds Threads.@threads for j in 2:KS.pSpace.ny-1
        for i in 1:KS.pSpace.nx
            sw = extract_last(ctr[i, j].sw, 2, mode=:view)
            reconstruct3!(
                sw,
                ctr[i, j-1].w,
                ctr[i, j].w,
                ctr[i, j+1].w,
                0.5 * (ctr[i, j-1].dy + ctr[i, j].dy),
                0.5 * (ctr[i, j].dy + ctr[i, j+1].dy),
                Symbol(KS.set.limiter),
            )
        end
    end

    #--- particle distribution function ---#
    # boundary
    @inbounds Threads.@threads for j in 1:KS.pSpace.ny
        s0L = extract_last(ctr[1, j].sh0, 1, mode=:view)
        reconstruct2!(
            s0L,
            ctr[1, j].h0,
            ctr[2, j].h0,
            0.5 * (ctr[1, j].dx + ctr[2, j].dx),
        )
        s1L = extract_last(ctr[1, j].sh1, 1, mode=:view)
        reconstruct2!(
            s1L,
            ctr[1, j].h1,
            ctr[2, j].h1,
            0.5 * (ctr[1, j].dx + ctr[2, j].dx),
        )
        s2L = extract_last(ctr[1, j].sh2, 1, mode=:view)
        reconstruct2!(
            s2L,
            ctr[1, j].h2,
            ctr[2, j].h2,
            0.5 * (ctr[1, j].dx + ctr[2, j].dx),
        )

        s0R = extract_last(ctr[KS.pSpace.nx, j].sh0, 1, mode=:view)
        reconstruct2!(
            s0R,
            ctr[KS.pSpace.nx-1, j].h0,
            ctr[KS.pSpace.nx, j].h0,
            0.5 * (ctr[KS.pSpace.nx-1, j].dx + ctr[KS.pSpace.nx, j].dx),
        )
        s1R = extract_last(ctr[KS.pSpace.nx, j].sh1, 1, mode=:view)
        reconstruct2!(
            s1R,
            ctr[KS.pSpace.nx-1, j].h1,
            ctr[KS.pSpace.nx, j].h1,
            0.5 * (ctr[KS.pSpace.nx-1, j].dx + ctr[KS.pSpace.nx, j].dx),
        )
        s2R = extract_last(ctr[KS.pSpace.nx, j].sh2, 1, mode=:view)
        reconstruct2!(
            s2R,
            ctr[KS.pSpace.nx-1, j].h2,
            ctr[KS.pSpace.nx, j].h2,
            0.5 * (ctr[KS.pSpace.nx-1, j].dx + ctr[KS.pSpace.nx, j].dx),
        )
    end

    @inbounds Threads.@threads for i in 1:KS.pSpace.nx
        s0D = extract_last(ctr[i, 1].sh0, 2, mode=:view)
        reconstruct2!(
            s0D,
            ctr[i, 1].h0,
            ctr[i, 2].h0,
            0.5 * (ctr[i, 1].dy + ctr[i, 2].dy),
        )
        s1D = extract_last(ctr[i, 1].sh1, 2, mode=:view)
        reconstruct2!(
            s1D,
            ctr[i, 1].h1,
            ctr[i, 2].h1,
            0.5 * (ctr[i, 1].dy + ctr[i, 2].dy),
        )
        s2D = extract_last(ctr[i, 1].sh2, 2, mode=:view)
        reconstruct2!(
            s2D,
            ctr[i, 1].h2,
            ctr[i, 2].h2,
            0.5 * (ctr[i, 1].dy + ctr[i, 2].dy),
        )

        s0U = extract_last(ctr[i, KS.pSpace.ny].sh0, 2, mode=:view)
        reconstruct2!(
            s0U,
            ctr[i, KS.pSpace.ny-1].h0,
            ctr[i, KS.pSpace.ny].h0,
            0.5 * (ctr[i, KS.pSpace.ny-1].dx + ctr[i, KS.pSpace.ny].dx),
        )
        s1U = extract_last(ctr[i, KS.pSpace.ny].sh1, 2, mode=:view)
        reconstruct2!(
            s1U,
            ctr[i, KS.pSpace.ny-1].h1,
            ctr[i, KS.pSpace.ny].h1,
            0.5 * (ctr[i, KS.pSpace.ny-1].dx + ctr[i, KS.pSpace.ny].dx),
        )
        s2U = extract_last(ctr[i, KS.pSpace.ny].sh2, 2, mode=:view)
        reconstruct2!(
            s2U,
            ctr[i, KS.pSpace.ny-1].h2,
            ctr[i, KS.pSpace.ny].h2,
            0.5 * (ctr[i, KS.pSpace.ny-1].dx + ctr[i, KS.pSpace.ny].dx),
        )
    end

    # inner
    @inbounds Threads.@threads for j in 1:KS.pSpace.ny
        for i in 2:KS.pSpace.nx-1
            sh0 = extract_last(ctr[i, j].sh0, 1, mode=:view)
            reconstruct3!(
                sh0,
                ctr[i-1, j].h0,
                ctr[i, j].h0,
                ctr[i+1, j].h0,
                0.5 * (ctr[i-1, j].dx + ctr[i, j].dx),
                0.5 * (ctr[i, j].dx + ctr[i+1, j].dx),
                Symbol(KS.set.limiter),
            )

            sh1 = extract_last(ctr[i, j].sh1, 1, mode=:view)
            reconstruct3!(
                sh1,
                ctr[i-1, j].h1,
                ctr[i, j].h1,
                ctr[i+1, j].h1,
                0.5 * (ctr[i-1, j].dx + ctr[i, j].dx),
                0.5 * (ctr[i, j].dx + ctr[i+1, j].dx),
                Symbol(KS.set.limiter),
            )

            sh2 = extract_last(ctr[i, j].sh2, 1, mode=:view)
            reconstruct3!(
                sh2,
                ctr[i-1, j].h2,
                ctr[i, j].h2,
                ctr[i+1, j].h2,
                0.5 * (ctr[i-1, j].dx + ctr[i, j].dx),
                0.5 * (ctr[i, j].dx + ctr[i+1, j].dx),
                Symbol(KS.set.limiter),
            )
        end
    end

    @inbounds Threads.@threads for j in 2:KS.pSpace.ny-1
        for i in 1:KS.pSpace.nx
            sh0 = extract_last(ctr[i, j].sh0, 2, mode=:view)
            reconstruct3!(
                sh0,
                ctr[i, j-1].h0,
                ctr[i, j].h0,
                ctr[i, j+1].h0,
                0.5 * (ctr[i, j-1].dy + ctr[i, j].dy),
                0.5 * (ctr[i, j].dy + ctr[i, j+1].dy),
                Symbol(KS.set.limiter),
            )

            sh1 = extract_last(ctr[i, j].sh1, 2, mode=:view)
            reconstruct3!(
                sh1,
                ctr[i, j-1].h1,
                ctr[i, j].h1,
                ctr[i, j+1].h1,
                0.5 * (ctr[i, j-1].dy + ctr[i, j].dy),
                0.5 * (ctr[i, j].dy + ctr[i, j+1].dy),
                Symbol(KS.set.limiter),
            )

            sh2 = extract_last(ctr[i, j].sh2, 2, mode=:view)
            reconstruct3!(
                sh2,
                ctr[i, j-1].h2,
                ctr[i, j].h2,
                ctr[i, j+1].h2,
                0.5 * (ctr[i, j-1].dy + ctr[i, j].dy),
                0.5 * (ctr[i, j].dy + ctr[i, j+1].dy),
                Symbol(KS.set.limiter),
            )
        end
    end

end


"""
Evolution

"""
function evolve!(
    KS::SolverSet,
    ctr::AbstractArray{<:AbstractControlVolume1D,1},
    face::Array{<:AbstractInterface1D,1},
    dt::Real;
    mode = :kfvs::Symbol,
)

    #if KS.set.case == "heat"
    #		flux_maxwell!(KS.ib.bcL, face[1], ctr[1], 1, dt)
    #    end

    if KS.set.space[3:end] == "1f1v"

        @inbounds Threads.@threads for i = 1:KS.pSpace.nx+1
            flux_kfvs!(
                face[i].fw,
                face[i].ff,
                ctr[i-1].f .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sf,
                ctr[i].f .- 0.5 .* ctr[i].dx .* ctr[i].sf,
                KS.vSpace.u,
                KS.vSpace.weights,
                dt,
                ctr[i-1].sf,
                ctr[i].sf,
            )
        end

    elseif KS.set.space[3:end] == "1f3v"

        @inbounds Threads.@threads for i = 1:KS.pSpace.nx+1
            flux_kfvs!(
                face[i].fw,
                face[i].ff,
                ctr[i-1].f .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sf,
                ctr[i].f .- 0.5 .* ctr[i].dx .* ctr[i].sf,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.w,
                KS.vSpace.weights,
                dt,
                ctr[i-1].sf,
                ctr[i].sf,
            )
        end

    elseif KS.set.space[3:end] == "2f1v"

        @inbounds Threads.@threads for i = 1:KS.pSpace.nx+1
            flux_kcu!(
                face[i].fw,
                face[i].fh,
                face[i].fb,
                ctr[i-1].w .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sw,
                ctr[i-1].h .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh,
                ctr[i-1].b .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sb,
                ctr[i].w .- 0.5 .* ctr[i].dx .* ctr[i].sw,
                ctr[i].h .- 0.5 .* ctr[i].dx .* ctr[i].sh,
                ctr[i].b .- 0.5 .* ctr[i].dx .* ctr[i].sb,
                KS.vSpace.u,
                KS.vSpace.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                dt,
            )
        end

    elseif KS.set.space[3:end] == "4f1v"

        if mode == :kcu
            @inbounds Threads.@threads for i = 1:KS.pSpace.nx+1
                flux_kcu!(
                    face[i].fw,
                    face[i].fh0,
                    face[i].fh1,
                    face[i].fh2,
                    face[i].fh3,
                    ctr[i-1].w .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sw,
                    ctr[i-1].h0 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh0,
                    ctr[i-1].h1 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh1,
                    ctr[i-1].h2 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh2,
                    ctr[i-1].h3 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh3,
                    ctr[i].w .- 0.5 .* ctr[i].dx .* ctr[i].sw,
                    ctr[i].h0 .- 0.5 .* ctr[i].dx .* ctr[i].sh0,
                    ctr[i].h1 .- 0.5 .* ctr[i].dx .* ctr[i].sh1,
                    ctr[i].h2 .- 0.5 .* ctr[i].dx .* ctr[i].sh2,
                    ctr[i].h3 .- 0.5 .* ctr[i].dx .* ctr[i].sh3,
                    KS.vSpace.u,
                    KS.vSpace.weights,
                    KS.gas.K,
                    KS.gas.γ,
                    KS.gas.mi,
                    KS.gas.ni,
                    KS.gas.me,
                    KS.gas.ne,
                    KS.gas.Kn[1],
                    dt,
                )
            end
        elseif mode == :kfvs
            @inbounds Threads.@threads for i = 1:KS.pSpace.nx+1
                flux_kfvs!(
                    face[i].fw,
                    face[i].fh0,
                    face[i].fh1,
                    face[i].fh2,
                    face[i].fh3,
                    ctr[i-1].h0 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh0,
                    ctr[i-1].h1 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh1,
                    ctr[i-1].h2 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh2,
                    ctr[i-1].h3 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh3,
                    ctr[i].h0 .- 0.5 .* ctr[i].dx .* ctr[i].sh0,
                    ctr[i].h1 .- 0.5 .* ctr[i].dx .* ctr[i].sh1,
                    ctr[i].h2 .- 0.5 .* ctr[i].dx .* ctr[i].sh2,
                    ctr[i].h3 .- 0.5 .* ctr[i].dx .* ctr[i].sh3,
                    KS.vSpace.u,
                    KS.vSpace.weights,
                    dt,
                    ctr[i-1].sh0,
                    ctr[i-1].sh1,
                    ctr[i-1].sh2,
                    ctr[i-1].sh3,
                    ctr[i].sh0,
                    ctr[i].sh1,
                    ctr[i].sh2,
                    ctr[i].sh3,
                )
            end
        end

        @inbounds Threads.@threads for i = 1:KS.pSpace.nx+1
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
                    KS.gas.Ap,
                    KS.gas.An,
                    KS.gas.D,
                    KS.gas.sol,
                    KS.gas.χ,
                    KS.gas.ν,
                    dt,
                )
            
        end

    elseif KS.set.space[3:end] == "3f2v"


    end

end


"""
Update flow variables

"""
function update!(
    KS::SolverSet,
    ctr::AbstractArray{ControlVolume1D1F,1},
    face::Array{Interface1D1F,1},
    dt::Real,
    residual::Array{<:AbstractFloat}; # 1D / 2D
    mode = :bgk::Symbol,
    bc = :extra::Symbol,
)

    sumRes = zeros(axes(KS.ib.wL))
    sumAvg = zeros(axes(KS.ib.wL))

    phase = KS.set.space[3:end]
    if phase == "1f1v"
        @inbounds Threads.@threads for i = 2:KS.pSpace.nx-1
            step!(
                face[i].fw,
                face[i].ff,
                ctr[i].w,
                ctr[i].prim,
                ctr[i].f,
                face[i+1].fw,
                face[i+1].ff,
                KS.vSpace.u,
                KS.vSpace.weights,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                ctr[i].dx,
                dt,
                sumRes,
                sumAvg,
                mode,
            )
        end
    elseif phase == "1f3v"
        @inbounds Threads.@threads for i = 2:KS.pSpace.nx-1
            step!(
                face[i].fw,
                face[i].ff,
                ctr[i].w,
                ctr[i].prim,
                ctr[i].f,
                face[i+1].fw,
                face[i+1].ff,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.w,
                KS.vSpace.weights,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                ctr[i].dx,
                dt,
                sumRes,
                sumAvg,
                mode,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.pSpace.nx) / (sumAvg[i] + 1.e-7)
    end

    if bc == :extra
        ng = 1 - first(eachindex(KS.pSpace.x))

        for i in 1:ng
            @. ctr[1-i].w = ctr[1].w
            @. ctr[1-i].prim = ctr[1].prim
            @. ctr[1-i].f = ctr[1].f

            @. ctr[KS.pSpace.nx+i].w = ctr[KS.pSpace.nx].w
            @. ctr[KS.pSpace.nx+i].prim = ctr[KS.pSpace.nx].prim
            @. ctr[KS.pSpace.nx+i].f = ctr[KS.pSpace.nx].f
        end
    end

end

function update!(
    KS::SolverSet,
    ctr::AbstractArray{ControlVolume1D2F,1},
    face::Array{Interface1D2F,1},
    dt::Real,
    residual::Array{<:AbstractFloat}; # 1D / 2D
    mode = :bgk::Symbol,
    bc = :extra::Symbol,
)

    sumRes = zeros(axes(KS.ib.wL))
    sumAvg = zeros(axes(KS.ib.wL))

    @inbounds Threads.@threads for i = 2:KS.pSpace.nx-1
        step!(
            face[i].fw,
            face[i].fh,
            face[i].fb,
            ctr[i].w,
            ctr[i].prim,
            ctr[i].h,
            ctr[i].b,
            face[i+1].fw,
            face[i+1].fh,
            face[i+1].fb,
            KS.vSpace.u,
            KS.vSpace.weights,
            KS.gas.K,
            KS.gas.γ,
            KS.gas.μᵣ,
            KS.gas.ω,
            KS.gas.Pr,
            ctr[i].dx,
            dt,
            sumRes,
            sumAvg,
            mode,
        )
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.pSpace.nx) / (sumAvg[i] + 1.e-7)
    end

    if bc == :extra
        ng = 1 - first(eachindex(KS.pSpace.x))

        for i in 1:ng
            @. ctr[1-i].w .= ctr[1].w
            @. ctr[1-i].prim .= ctr[1].prim
            @. ctr[1-i].h .= ctr[1].h
            @. ctr[1-i].b .= ctr[1].b

            @. ctr[KS.pSpace.nx+i].w .= ctr[KS.pSpace.nx].w
            @. ctr[KS.pSpace.nx+i].prim .= ctr[KS.pSpace.nx].prim
            @. ctr[KS.pSpace.nx+i].h .= ctr[KS.pSpace.nx].h
            @. ctr[KS.pSpace.nx+i].b .= ctr[KS.pSpace.nx].b
        end
    end

end

function update!(
    KS::SolverSet,
    ctr::AbstractArray{ControlVolume1D3F,1},
    face::Array{Interface1D3F,1},
    dt::Real,
    residual::Array{<:AbstractFloat}; # 1D / 2D
    mode = :bgk::Symbol,
    bc = :extra::Symbol,
)

    sumRes = zeros(axes(KS.ib.wL))
    sumAvg = zeros(axes(KS.ib.wL))

    @inbounds Threads.@threads for i = 1:KS.pSpace.nx
        step!(
            KS,
            face[i],
            ctr[i],
            face[i+1],
            dt,
            sumRes,
            sumAvg,
        )
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.pSpace.nx) / (sumAvg[i] + 1.e-7)
    end

    if bc == :extra
        ng = 1 - first(eachindex(KS.pSpace.x))

        for i in 1:ng
            @. ctr[1-i].w .= ctr[1].w
            @. ctr[1-i].prim .= ctr[1].prim

            @. ctr[KS.pSpace.nx+i].w .= ctr[KS.pSpace.nx].w
            @. ctr[KS.pSpace.nx+i].prim .= ctr[KS.pSpace.nx].prim
        end

        for i in 1:ng
            ctr[1-i].h0 = deepcopy(ctr[1].h0)
            ctr[1-i].h1 = deepcopy(ctr[1].h1)
            ctr[1-i].h2 = deepcopy(ctr[1].h2)
            ctr[1-i].E = deepcopy(ctr[1].E)
            ctr[1-i].B = deepcopy(ctr[1].B)
            ctr[1-i].ϕ = deepcopy(ctr[1].ϕ)
            ctr[1-i].ψ = deepcopy(ctr[1].ψ)
            ctr[1-i].lorenz = deepcopy(ctr[1].lorenz)

            ctr[KS.pSpace.nx+i].h0 = deepcopy(ctr[KS.pSpace.nx].h0)
            ctr[KS.pSpace.nx+i].h1 = deepcopy(ctr[KS.pSpace.nx].h1)
            ctr[KS.pSpace.nx+i].h2 = deepcopy(ctr[KS.pSpace.nx].h2)
            ctr[KS.pSpace.nx+i].E = deepcopy(ctr[KS.pSpace.nx].E)
            ctr[KS.pSpace.nx+i].B = deepcopy(ctr[KS.pSpace.nx].B)
            ctr[KS.pSpace.nx+i].ϕ = deepcopy(ctr[KS.pSpace.nx].ϕ)
            ctr[KS.pSpace.nx+i].ψ = deepcopy(ctr[KS.pSpace.nx].ψ)
            ctr[KS.pSpace.nx+i].lorenz = deepcopy(ctr[KS.pSpace.nx].lorenz)
        end

    end

end

function update!(
    KS::SolverSet,
    ctr::AbstractArray{ControlVolume1D4F,1},
    face::Array{Interface1D4F,1},
    dt::Real,
    residual::Array{<:AbstractFloat}; # 1D / 2D
    mode = :bgk::Symbol,
    bc = :extra::Symbol,
)

    sumRes = zeros(axes(KS.ib.wL))
    sumAvg = zeros(axes(KS.ib.wL))

    @inbounds Threads.@threads for i = 1:KS.pSpace.nx
        step!(
            KS,
            face[i],
            ctr[i],
            face[i+1],
            dt,
            sumRes,
            sumAvg,
        )
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.pSpace.nx) / (sumAvg[i] + 1.e-7)
    end

    if bc == :extra
        ng = 1 - first(eachindex(KS.pSpace.x))

        for i in 1:ng
            @. ctr[1-i].w .= ctr[1].w
            @. ctr[1-i].prim .= ctr[1].prim

            @. ctr[KS.pSpace.nx+i].w .= ctr[KS.pSpace.nx].w
            @. ctr[KS.pSpace.nx+i].prim .= ctr[KS.pSpace.nx].prim
        end

        for i in 1:ng
            ctr[1-i].h0 = deepcopy(ctr[1].h0)
            ctr[1-i].h1 = deepcopy(ctr[1].h1)
            ctr[1-i].h2 = deepcopy(ctr[1].h2)
            ctr[1-i].h3 = deepcopy(ctr[1].h3)
            ctr[1-i].E = deepcopy(ctr[1].E)
            ctr[1-i].B = deepcopy(ctr[1].B)
            ctr[1-i].ϕ = deepcopy(ctr[1].ϕ)
            ctr[1-i].ψ = deepcopy(ctr[1].ψ)
            ctr[1-i].lorenz = deepcopy(ctr[1].lorenz)

            ctr[KS.pSpace.nx+i].h0 = deepcopy(ctr[KS.pSpace.nx].h0)
            ctr[KS.pSpace.nx+i].h1 = deepcopy(ctr[KS.pSpace.nx].h1)
            ctr[KS.pSpace.nx+i].h2 = deepcopy(ctr[KS.pSpace.nx].h2)
            ctr[KS.pSpace.nx+i].h3 = deepcopy(ctr[KS.pSpace.nx].h3)
            ctr[KS.pSpace.nx+i].E = deepcopy(ctr[KS.pSpace.nx].E)
            ctr[KS.pSpace.nx+i].B = deepcopy(ctr[KS.pSpace.nx].B)
            ctr[KS.pSpace.nx+i].ϕ = deepcopy(ctr[KS.pSpace.nx].ϕ)
            ctr[KS.pSpace.nx+i].ψ = deepcopy(ctr[KS.pSpace.nx].ψ)
            ctr[KS.pSpace.nx+i].lorenz = deepcopy(ctr[KS.pSpace.nx].lorenz)
        end

    end

end
