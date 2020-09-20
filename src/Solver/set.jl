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

        # generate variables from configuration file
        dict = read_dict(configfilename)
        for key in keys(dict)
            s = Symbol(key)
            @eval $s = $(dict[key])
        end

        # set
        set = Setup(case, space, nSpecies, interpOrder, limiter, cfl, maxTime)

        # physical space
        Dx = parse(Int, space[1])
        if Dx == 1
            pSpace = PSpace1D(x0, x1, nx, pMeshType, nxg)
        elseif Dx == 2
            pSpace = PSpace2D(x0, x1, nx, y0, y1, ny, pMeshType, nxg, nyg)
        else
            throw("No preset available for 3D simulation, please set it up manually.")
        end

        # velocity space
        Dv = parse(Int, space[5])
        if Dv == 1
            if nSpecies == 1
                vSpace = VSpace1D(umin, umax, nu, vMeshType, nug)
            elseif nSpecies == 2
                ue0 = umin * sqrt(mi / me)
                ue1 = umax * sqrt(mi / me)
                vSpace = MVSpace1D(umin, umax, ue0, ue1, nu, vMeshType, nug)
            else
                throw("The velocity space only supports up to two species.")
            end
        elseif Dv == 2
            if nSpecies == 1
                vSpace = VSpace2D(umin, umax, nu, vmin, vmax, nv, vMeshType, nug, nvg)
            elseif nSpecies == 2
                ue0 = umin * sqrt(mi / me)
                ue1 = umax * sqrt(mi / me)
                ve0 = vmin * sqrt(mi / me)
                ve1 = vmax * sqrt(mi / me)
                vSpace = MVSpace2D(
                            umin, 
                            umax, 
                            ue0, 
                            ue1, 
                            nu, 
                            vmin, 
                            vmax, 
                            ve0, 
                            ve1, 
                            nv, 
                            vMeshType, 
                            nug, 
                            nvg,
                        )
            else
                throw("The velocity space only supports up to two species.")
            end
        elseif Dv == 3
            if nSpecies == 1
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
            elseif nSpecies == 2
                ue0 = umin * sqrt(mi / me)
                ue1 = umax * sqrt(mi / me)
                ve0 = vmin * sqrt(mi / me)
                ve1 = vmax * sqrt(mi / me)
                we0 = wmin * sqrt(mi / me)
                we1 = wmax * sqrt(mi / me)
                vSpace = MVSpace3D(
                            umin, 
                            umax, 
                            ue0, 
                            ue1, 
                            nu, 
                            vmin, 
                            vmax, 
                            ve0, 
                            ve1, 
                            nv, 
                            wmin,
                            wmax,
                            we0,
                            we1,
                            nw,
                            vMeshType, 
                            nug, 
                            nvg,
                            nwg,
                        )
            else
                throw("The velocity space only supports up to two species.")
            end
        end

        # gas property
        γD = map(parse(Int, space[3]), parse(Int, space[5])) do x, y # (x)f(y)v
            if x == 1 # 1f
                if y >= 3 # 3v
                    return 3
                else # 1v
                    return Dx
                end
            elseif x == 2 # 2f
                return Dx
            elseif x >= 3 # 3f / 4f
                return 3
            else
                return nothing
            end
        end
        γ = heat_capacity_ratio(inK, γD)

        if nSpecies == 1
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
        elseif nSpecies == 2
            if !(@isdefined rL) && !(@isdefined echi) # undefined plasma args
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
            else
                if Dx == 1
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
                elseif Dx == 2
                    gas = Plasma2D(
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
                else
                    throw("The plasma property only supports up to 2D case.")
                end
            end
        else
            throw("The gas property only supports up to two species.")
        end

        if @isdefined uLid
            ib = set_ib(set, vSpace, gas, uLid, vLid, tLid)
        else
            ib = set_ib(set, vSpace, gas)
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
Generate AbstractIB

"""
function set_ib(set, vSpace, gas, uLid = 0.15, vLid = 0.0, tLid = 1.0)

    if set.case == "shock"

        if set.space[3:end] == "1f1v"
            wL, primL, fL, bcL, wR, primR, fR, bcR = ib_rh(gas.Ma, gas.γ, vSpace.u)
            ib = IB1F(wL, primL, fL, bcL, wR, primR, fR, bcR)
        elseif set.space[3:end] == "1f3v"
            wL, primL, fL, bcL, wR, primR, fR, bcR =
                ib_rh(gas.Ma, gas.γ, vSpace.u, vSpace.v, vSpace.w)
            ib = IB1F(wL, primL, fL, bcL, wR, primR, fR, bcR)
        elseif set.space[3:end] == "2f1v"
            wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR =
                ib_rh(gas.mach, gas.γ, vSpace.u, gas.K)
            ib = IB2F(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR)
        end

    elseif set.case == "sod"

        if set.space[3:end] == "1f1v"
            wL, primL, fL, bcL, wR, primR, fR, bcR = ib_sod(gas.γ, vSpace.u)
            ib = IB1F(wL, primL, fL, bcL, wR, primR, fR, bcR)
        elseif set.space[3:end] == "1f3v"
            wL, primL, fL, bcL, wR, primR, fR, bcR =
                ib_sod(gas.γ, vSpace.u, vSpace.v, vSpace.w)
            ib = IB1F(wL, primL, fL, bcL, wR, primR, fR, bcR)
        elseif set.space[3:end] == "2f1v"
            wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR =
                ib_sod(gas.γ, vSpace.u, gas.K)
            ib = IB2F(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR)
        end

    elseif set.case == "brio-wu"

        @assert set.space[3:end] == "4f1v"

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
        lorenzR = ib_briowu(gas.γ, vSpace.u, gas.mi, gas.me)

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

    elseif set.case == "cavity"

        if set.space[3:end] == "1f2v"
            wL, primL, fL, bcL, wR, primR, fR, bcR, bcU, bcD =
                ib_cavity(gas.γ, uLid, vLid, tLid, vSpace.u, vSpace.v)
            ib = IB1F(wL, primL, fL, bcL, wR, primR, fR, bcR, bcU, bcD)
        elseif set.space[3:end] == "2f2v"
            wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR, bcU, bcD =
                ib_cavity(gas.γ, uLid, vLid, tLid, vSpace.u, vSpace.v, gas.K)
            ib = IB2F(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR, bcU, bcD)
        end

    else

        throw("No default ib available, please set it up manually.")

    end

    return ib

end