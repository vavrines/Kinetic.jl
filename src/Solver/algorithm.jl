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
            @inbounds Threads.@threads for i = 2:KS.pSpace.nx
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
            @inbounds Threads.@threads for i = 2:KS.pSpace.nx
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

        @inbounds Threads.@threads for i = 2:KS.pSpace.nx
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

    end

end


"""
Update flow variables

"""
function update!(
    KS::SolverSet,
    ctr::AbstractArray{<:AbstractControlVolume1D,1},
    face::Array{<:AbstractInterface1D,1},
    dt::Real,
    residual::Array{<:AbstractFloat}; # 1D / 2D
    collision = :bgk::Symbol,
)

    sumRes = zeros(axes(KS.ib.wL))
    sumAvg = zeros(axes(KS.ib.wL))

    @inbounds Threads.@threads for i = 2:KS.pSpace.nx-1
        if KS.set.space == "1d1f1v"
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
                collision,
            )
        elseif KS.set.space == "1d1f3v"
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
                collision,
            )
        elseif KS.set.space == "1d2f1v"
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
                collision,
            )
        elseif KS.set.space[3:end] == "4f1v"
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

    end

    #if KS.set.case == "heat"
    #    ctr[KS.pSpace.nx].w = deepcopy(ctr[KS.pSpace.nx-1].w)
    #    ctr[KS.pSpace.nx].prim = deepcopy(ctr[KS.pSpace.nx-1].prim)
    #    ctr[KS.pSpace.nx].h = deepcopy(ctr[KS.pSpace.nx-1].h)
    #end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.pSpace.nx) / (sumAvg[i] + 1.e-7)
    end

end


# ------------------------------------------------------------
# Time stepping
# ------------------------------------------------------------

"""
Update flow variables with finite volume formulation

`step!(fwL, w, prim, fwR, γ, dx, RES, AVG)`

"""
function step!(
    fwL::Array{<:AbstractFloat,1},
    w::Array{<:AbstractFloat,1},
    prim::Array{<:AbstractFloat,1},
    fwR::Array{<:AbstractFloat,1},
    γ::Real,
    dx::Real,
    RES::Array{<:AbstractFloat,1},
    AVG::Array{<:AbstractFloat,1},
)

    #--- store W^n and calculate H^n,\tau^n ---#
    w_old = deepcopy(w)

    #--- update W^{n+1} ---#
    @. w += (fwL - fwR) / dx
    prim .= conserve_prim(w, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

end


#--- 1D1F1V ---#
function step!(
    fwL::Array{<:AbstractFloat,1},
    ffL::AbstractArray{<:AbstractFloat,1},
    w::Array{<:AbstractFloat,1},
    prim::Array{<:AbstractFloat,1},
    f::AbstractArray{<:AbstractFloat,1},
    fwR::Array{<:AbstractFloat,1},
    ffR::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    weights::AbstractArray{<:AbstractFloat,1},
    γ::Real,
    μᵣ::Real,
    ω::Real,
    Pr::Real,
    dx::Real,
    dt::Real,
    RES::Array{<:AbstractFloat,1},
    AVG::Array{<:AbstractFloat,1},
    collision = :bgk::Symbol,
)

    #--- store W^n and calculate H^n,\tau^n ---#
    w_old = deepcopy(w)

    if collision == :shakhov
        q = heat_flux(f, prim, u, weights)
        M_old = maxwellian(u, prim)
        S = shakhov(u, M_old, q, prim, Pr)
    else
        S = zeros(axes(f))
    end

    #--- update W^{n+1} ---#
    @. w += (fwL - fwR) / dx
    prim .= conserve_prim(w, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    #--- calculate M^{n+1} and tau^{n+1} ---#
    M = maxwellian(u, prim)
    M .+= S
    τ = vhs_collision_time(prim, μᵣ, ω)

    #--- update distribution function ---#
    for i in eachindex(u)
        f[i] = (f[i] + (ffL[i] - ffR[i]) / dx + dt / τ * M[i]) / (1.0 + dt / τ)
    end

end


#--- 1D1F3V ---#
function step!(
    fwL::Array{<:AbstractFloat,1},
    ffL::AbstractArray{<:AbstractFloat,3},
    w::Array{<:AbstractFloat,1},
    prim::Array{<:AbstractFloat,1},
    f::AbstractArray{<:AbstractFloat,3},
    fwR::Array{<:AbstractFloat,1},
    ffR::AbstractArray{<:AbstractFloat,3},
    uVelo::AbstractArray{<:AbstractFloat,3},
    vVelo::AbstractArray{<:AbstractFloat,3},
    wVelo::AbstractArray{<:AbstractFloat,3}, # avoid conflict with w
    weights::AbstractArray{<:AbstractFloat,3},
    γ::Real,
    μᵣ::Real,
    ω::Real,
    Pr::Real,
    dx::Real,
    dt::Real,
    RES::Array{<:AbstractFloat,1},
    AVG::Array{<:AbstractFloat,1},
    collision = :bgk::Symbol,
)

    #--- store W^n and calculate shakhov term ---#
    w_old = deepcopy(w)

    if collision == :shakhov
        q = heat_flux(f, prim, uVelo, vVelo, wVelo, weights)
        M_old = maxwellian(uVelo, vVelo, wVelo, prim)
        S = shakhov(uVelo, vVelo, wVelo, M_old, q, prim, Pr, K)
    else
        S = zeros(axes(f))
    end

    #--- update W^{n+1} ---#
    @. w += (fwL - fwR) / dx
    prim .= conserve_prim(w, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    #--- calculate M^{n+1} and tau^{n+1} ---#
    M = maxwellian(uVelo, vVelo, wVelo, prim)
    M .+= S
    τ = vhs_collision_time(prim, μᵣ, ω)

    #--- update distribution function ---#
    for k in axes(wVelo, 3), j in axes(vVelo, 2), i in axes(uVelo, 1)
        f[i, j, k] =
            (f[i, j, k] + (ffL[i, j, k] - ffR[i, j, k]) / dx + dt / τ * M[i, j, k]) /
            (1.0 + dt / τ)
    end

end


#--- 1D2F1V ---#
function step!(
    fwL::Array{<:AbstractFloat,1},
    fhL::AbstractArray{<:AbstractFloat,1},
    fbL::AbstractArray{<:AbstractFloat,1},
    w::Array{<:AbstractFloat,1},
    prim::Array{<:AbstractFloat,1},
    h::AbstractArray{<:AbstractFloat,1},
    b::AbstractArray{<:AbstractFloat,1},
    fwR::Array{<:AbstractFloat,1},
    fhR::AbstractArray{<:AbstractFloat,1},
    fbR::AbstractArray{<:AbstractFloat,1},
    u::AbstractArray{<:AbstractFloat,1},
    weights::AbstractArray{<:AbstractFloat,1},
    K::Real,
    γ::Real,
    μᵣ::Real,
    ω::Real,
    Pr::Real,
    dx::Real,
    dt::Real,
    RES::Array{<:AbstractFloat,1},
    AVG::Array{<:AbstractFloat,1},
    collision = :bgk::Symbol,
)

    #--- store W^n and calculate shakhov term ---#
    w_old = deepcopy(w)

    if collision == :shakhov
        q = heat_flux(h, b, prim, u, weights)
        MH_old = maxwellian(u, prim)
        MB_old = MH_old .* K ./ (2.0 * prim[end])
        SH, SB = shakhov(u, MH_old, MB_old, q, prim, Pr, K)
    else
        SH = zeros(axes(h))
        SB = zeros(axes(b))
    end

    #--- update W^{n+1} ---#
    @. w += (fwL - fwR) / dx
    prim .= conserve_prim(w, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    #--- calculate M^{n+1} and tau^{n+1} ---#
    MH = maxwellian(u, prim)
    MB = MH .* K ./ (2.0 * prim[end])
    MH .+= SH
    MB .+= SB
    τ = vhs_collision_time(prim, μᵣ, ω)

    #--- update distribution function ---#
    for i in eachindex(u)
        h[i] = (h[i] + (fhL[i] - fhR[i]) / dx + dt / τ * MH[i]) / (1.0 + dt / τ)
        b[i] = (b[i] + (fbL[i] - fbR[i]) / dx + dt / τ * MB[i]) / (1.0 + dt / τ)
    end

end


#--- 1D4F1V ---#
function step!(
    KS::SolverSet,
    faceL::Interface1D4F,
    cell::ControlVolume1D4F,
    faceR::Interface1D4F,
    dt::AbstractFloat,
    RES::Array{<:AbstractFloat,2},
    AVG::Array{<:AbstractFloat,2},
)

    #--- update conservative flow variables: step 1 ---#
    # w^n
    w_old = deepcopy(cell.w)
    prim_old = deepcopy(cell.prim)

    # flux -> w^{n+1}
    @. cell.w += (faceL.fw - faceR.fw) / cell.dx
    cell.prim .= mixture_conserve_prim(cell.w, KS.gas.γ)

    # temperature protection
    if min(minimum(cell.prim[5, 1]), minimum(cell.prim[5, 2])) < 0
        println("warning: temperature update is negative")
        cell.w .= w_old
        cell.prim .= prim_old
    end

    #=
    # source -> w^{n+1}
    # DifferentialEquations.jl
    tau = get_tau(cell.prim, KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.Kn[1])
    for j in axes(wRan, 2)
    prob = ODEProblem( mixture_source, 
        vcat(cell.w[1:5,j,1], cell.w[1:5,j,2]),
        dt,
        (tau[1], tau[2], KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.Kn[1], KS.gas.γ) )
    sol = solve(prob, Rosenbrock23())

    cell.w[1:5,j,1] .= sol[end][1:5]
    cell.w[1:5,j,2] .= sol[end][6:10]
    for k=1:2
    cell.prim[:,j,k] .= Kinetic.conserve_prim(cell.w[:,j,k], KS.gas.γ)
    end
    end
    =#
    #=
    # explicit
    tau = get_tau(cell.prim, KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.Kn[1])
    mprim = get_mixprim(cell.prim, tau, KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.Kn[1])
    mw = get_conserved(mprim, KS.gas.γ)
    for k=1:2
    cell.w[:,:,k] .+= (mw[:,:,k] .- w_old[:,:,k]) .* dt ./ tau[k]
    end
    cell.prim .= get_primitive(cell.w, KS.gas.γ)
    =#

    #--- update electromagnetic variables ---#
    # flux -> E^{n+1} & B^{n+1}
    cell.E[1] -= dt * (faceL.femR[1] + faceR.femL[1]) / cell.dx
    cell.E[2] -= dt * (faceL.femR[2] + faceR.femL[2]) / cell.dx
    cell.E[3] -= dt * (faceL.femR[3] + faceR.femL[3]) / cell.dx
    cell.B[1] -= dt * (faceL.femR[4] + faceR.femL[4]) / cell.dx
    cell.B[2] -= dt * (faceL.femR[5] + faceR.femL[5]) / cell.dx
    cell.B[3] -= dt * (faceL.femR[6] + faceR.femL[6]) / cell.dx
    cell.ϕ -= dt * (faceL.femR[7] + faceR.femL[7]) / cell.dx
    cell.ψ -= dt * (faceL.femR[8] + faceR.femL[8]) / cell.dx

    # source -> ϕ
    #@. cell.ϕ += dt * (cell.w[1,:,1] / KS.gas.mi - cell.w[1,:,2] / KS.gas.me) / (KS.gas.lD^2 * KS.gas.rL)

    # source -> U^{n+1}, E^{n+1} and B^{n+1}
    mr = KS.gas.mi / KS.gas.me
    A, b = em_coefficients(cell.prim, cell.E, cell.B, mr, KS.gas.lD, KS.gas.rL, dt)
    x = A \ b

    #--- calculate lorenz force ---#
    cell.lorenz[1,1] = 0.5 * (x[1] + cell.E[1] + (cell.prim[3,1] + x[5]) * cell.B[3] - (cell.prim[4,1] + x[6]) * cell.B[2]) / KS.gas.rL
    cell.lorenz[2,1] = 0.5 * (x[2] + cell.E[2] + (cell.prim[4,1] + x[6]) * cell.B[1] - (cell.prim[2,1] + x[4]) * cell.B[3]) / KS.gas.rL
    cell.lorenz[3,1] = 0.5 * (x[3] + cell.E[3] + (cell.prim[2,1] + x[4]) * cell.B[2] - (cell.prim[3,1] + x[5]) * cell.B[1]) / KS.gas.rL
    cell.lorenz[1,2] = -0.5 * (x[1] + cell.E[1] + (cell.prim[3,2] + x[8]) * cell.B[3] - (cell.prim[4,2] + x[9]) * cell.B[2]) * mr / KS.gas.rL
    cell.lorenz[2,2] = -0.5 * (x[2] + cell.E[2] + (cell.prim[4,2] + x[9]) * cell.B[1] - (cell.prim[2,2] + x[7]) * cell.B[3]) * mr / KS.gas.rL
    cell.lorenz[3,2] = -0.5 * (x[3] + cell.E[3] + (cell.prim[2,2] + x[7]) * cell.B[2] - (cell.prim[3,2] + x[8]) * cell.B[1]) * mr / KS.gas.rL

    cell.E[1] = x[1]
    cell.E[2] = x[2]
    cell.E[3] = x[3]

    #--- update conservative flow variables: step 2 ---#
    cell.prim[2, 1] = x[4]
    cell.prim[3, 1] = x[5]
    cell.prim[4, 1] = x[6]
    cell.prim[2, 2] = x[7]
    cell.prim[3, 2] = x[8]
    cell.prim[4, 2] = x[9]

    cell.w .= Kinetic.mixture_prim_conserve(cell.prim, KS.gas.γ)

    #--- update particle distribution function ---#
    # flux -> f^{n+1}
    @. cell.h0 += (faceL.fh0 - faceR.fh0) / cell.dx
    @. cell.h1 += (faceL.fh1 - faceR.fh1) / cell.dx
    @. cell.h2 += (faceL.fh2 - faceR.fh2) / cell.dx
    @. cell.h3 += (faceL.fh3 - faceR.fh3) / cell.dx

    # force -> f^{n+1} : step 1
    for j in axes(cell.h0, 2)
        _h0 = @view cell.h0[:, j]
        _h1 = @view cell.h1[:, j]
        _h2 = @view cell.h2[:, j]
        _h3 = @view cell.h3[:, j]

        shift_pdf!(_h0, cell.lorenz[1, j], KS.vSpace.du[1, j], dt)
        shift_pdf!(_h1, cell.lorenz[1, j], KS.vSpace.du[1, j], dt)
        shift_pdf!(_h2, cell.lorenz[1, j], KS.vSpace.du[1, j], dt)
        shift_pdf!(_h3, cell.lorenz[1, j], KS.vSpace.du[1, j], dt)
    end

    # force -> f^{n+1} : step 2
    for k in axes(cell.h1, 3)
        @. cell.h3[:,k] += 2. * dt * cell.lorenz[2,k] * cell.h1[:,k] + (dt * cell.lorenz[2,k])^2 * cell.h0[:,k] +
                                2. * dt * cell.lorenz[3,k] * cell.h2[:,k] + (dt * cell.lorenz[3,k])^2 * cell.h0[:,k]
        @. cell.h2[:,k] += dt * cell.lorenz[3,k] * cell.h0[:,k]
        @. cell.h1[:,k] += dt * cell.lorenz[2,k] * cell.h0[:,k]
    end

    # source -> f^{n+1}
    tau = aap_hs_collision_time(
        cell.prim,
        KS.gas.mi,
        KS.gas.ni,
        KS.gas.me,
        KS.gas.ne,
        KS.gas.Kn[1],
    )

    # interspecies interaction
    prim = deepcopy(cell.prim)
    #for j in axes(prim, 2)
    #    prim[:,j,:] .= aap_hs_prim(cell.prim[:,j,:], tau, KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.Kn[1])
    #end

    g = mixture_maxwellian(KS.vSpace.u, prim)

    # BGK term
    Mu, Mv, Mw, MuL, MuR = mixture_gauss_moments(prim, KS.gas.K)
    for k in axes(cell.h0, 3)
        @. cell.h0[:, k] =
            (cell.h0[:, k] + dt / tau[k] * g[:, k]) / (1.0 + dt / tau[k])
        @. cell.h1[:, k] =
            (cell.h1[:, k] + dt / tau[k] * Mv[1, k] * g[:, k]) /
            (1.0 + dt / tau[k])
        @. cell.h2[:, k] =
            (cell.h2[:, k] + dt / tau[k] * Mw[1, k] * g[:, k]) /
            (1.0 + dt / tau[k])
        @. cell.h3[:, k] =
            (cell.h3[:, k] + dt / tau[k] * (Mv[2, k] + Mw[2, k]) * g[:, k]) /
            (1.0 + dt / tau[k])
    end

    #--- record residuals ---#
    @. RES += (w_old - cell.w)^2
    @. AVG += abs(cell.w)

end
