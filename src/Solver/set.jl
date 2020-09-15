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
