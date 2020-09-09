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

                gas = GasProperty(
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

                gas = GasProperty(
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

                gas = PlasmaProperty(
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

                gas = GasProperty(
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
            vmax = max(maximum(KS.quad.u1), maximum(abs.(prim[2, :]))) + sos
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
        ctr[i].sw .= reconstruct3(
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
        if KS.set.space[1:4] == "1d1f"
            ctr[i].sf .= reconstruct3(
                ctr[i-1].f,
                ctr[i].f,
                ctr[i+1].f,
                0.5 * (ctr[i-1].dx + ctr[i].dx),
                0.5 * (ctr[i].dx + ctr[i+1].dx),
                Symbol(KS.set.limiter),
            )
        elseif KS.set.space[1:4] == "1d2f"
            ctr[i].sh .= reconstruct3(
                ctr[i-1].h,
                ctr[i].h,
                ctr[i+1].h,
                0.5 * (ctr[i-1].dx + ctr[i].dx),
                0.5 * (ctr[i].dx + ctr[i+1].dx),
                Symbol(KS.set.limiter),
            )
            ctr[i].sb .= reconstruct3(
                ctr[i-1].b,
                ctr[i].b,
                ctr[i+1].b,
                0.5 * (ctr[i-1].dx + ctr[i].dx),
                0.5 * (ctr[i].dx + ctr[i+1].dx),
                Symbol(KS.set.limiter),
            )
        elseif KS.set.space[1:4] == "1d4f"
            ctr[i].sh0 .= reconstruct3(
                ctr[i-1].h0,
                ctr[i].h0,
                ctr[i+1].h0,
                0.5 * (ctr[i-1].dx + ctr[i].dx),
                0.5 * (ctr[i].dx + ctr[i+1].dx),
                Symbol(KS.set.limiter),
            )
            ctr[i].sh1 .= reconstruct3(
                ctr[i-1].h1,
                ctr[i].h1,
                ctr[i+1].h1,
                0.5 * (ctr[i-1].dx + ctr[i].dx),
                0.5 * (ctr[i].dx + ctr[i+1].dx),
                Symbol(KS.set.limiter),
            )
            ctr[i].sh2 .= reconstruct3(
                ctr[i-1].h2,
                ctr[i].h2,
                ctr[i+1].h2,
                0.5 * (ctr[i-1].dx + ctr[i].dx),
                0.5 * (ctr[i].dx + ctr[i+1].dx),
                Symbol(KS.set.limiter),
            )
            ctr[i].sh3 .= reconstruct3(
                ctr[i-1].h3,
                ctr[i].h3,
                ctr[i+1].h3,
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
    dt::Real,
)

    #if KS.set.case == "heat"
    #		flux_maxwell!(KS.ib.bcL, face[1], ctr[1], 1, dt)
    #    end

    if KS.set.space == "1d1f1v"

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

    elseif KS.set.space == "1d1f3v"

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

    elseif KS.set.space[1:4] == "1d2f"

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

    elseif KS.set.space[1:4] == "1d4f"

        if KS.set.nSpecies == 2
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
                    KS.gas.knudsen[1],
                    dt,
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
    residual::Array{<:AbstractFloat,1};
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
