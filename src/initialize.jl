# ============================================================
# Initializer of Simulation
# ============================================================


export initialize,
       init_fvm


# ------------------------------------------------------------
# Initialize solver from input file
# ------------------------------------------------------------
function initialize(configfilename::String)

    println("==============================================================")
	println("Kinetic.jl")
	println("A Software Package for Kinetic Theory and Simulation")
    println("==============================================================")
    println("")
    println("reading configurations from $(configfilename)")
    println("")
    print("initializeing solver: ")

    if configfilename[end-2:end] == "txt"
        
        allowed = ["space"]
        D = read_dict(configfilename, allowed)
        space = D["space"]
        dim = parse(Int, space[1])

        if dim == 1
            println("1D solver")
            ks = SolverSet1D(configfilename)
        elseif dim == 2
            println("2D solver")
            #ks = SolverSet2D(configfilename)
        elseif dim == 3
            # coming soon
        end

        ctr, face = init_fvm(ks)

        return ks, ctr, face, 0.

    elseif configfilename[end-3:end] == "jld2"
        
        _1, _2, _3 = @load configfilename KS ctr t
        ks, ctr, simTime = eval(_1), eval(_2), eval(_3)

        face = init_fvm(ks)[2]
        
        return ks, ctr, face, simTime

    end

end


# ------------------------------------------------------------
# Initialize finite volume method
# ------------------------------------------------------------
function init_fvm(KS::SolverSet1D)

    dim = parse(Int, KS.set.space[1])

    if KS.set.space == "1d1f"
        #ctr = Array{ControlVolume1D1F}(undef, KS.pSpace.nx)
        ctr = OffsetArray{ControlVolume1D1F}(undef, eachindex(KS.pSpace.x)) # with ghost cells
        face = Array{Interface1D1F}(undef, KS.pSpace.nx+1)
        
        for i in eachindex(ctr)
            # shock problems
            if KS.set.case == "shock"
                if i <= KS.pSpace.nx÷2
                    ctr[i] = ControlVolume1D1F(KS.pSpace.x[i], KS.pSpace.dx[i], KS.ib.wL, KS.ib.primL, KS.ib.hL)
                else
                    ctr[i] = ControlVolume1D1F(KS.pSpace.x[i], KS.pSpace.dx[i], KS.ib.wR, KS.ib.primR, KS.ib.hR)
                end
            end
        end
    
        for i=1:KS.pSpace.nx+1
            face[i] = Interface1D1F(KS.ib.hL)
        end

    else
    end

    return ctr, face

end


