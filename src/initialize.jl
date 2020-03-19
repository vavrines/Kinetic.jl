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
        
        ks = SolverSet(configfilename)
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
function init_fvm(KS::SolverSet)

    if KS.set.nSpecies == 1
    
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

    elseif KS.set.nSpecies == 2

        if KS.set.space == "1d4f"
            #ctr = Array{ControlVolume1D1F}(undef, KS.pSpace.nx)
            ctr = OffsetArray{MControlVolume1D4F}(undef, eachindex(KS.pSpace.x)) # with ghost cells
            face = Array{MInterface1D4F}(undef, KS.pSpace.nx+1)
            
            for i in eachindex(ctr)
                # shock problems
                if KS.set.case == "brio-wu"
                    if i <= KS.pSpace.nx÷2
                        ctr[i] = MControlVolume1D4F( KS.pSpace.x[i], KS.pSpace.dx[i], 
                                                     KS.ib.wL, KS.ib.primL, KS.ib.h0L, KS.ib.h1L, 
                                                     KS.ib.h2L, KS.ib.h3L, KS.ib.EL, KS.ib.BL, KS.ib.lorenzL )
                    else
                        ctr[i] = MControlVolume1D4F( KS.pSpace.x[i], KS.pSpace.dx[i], 
                                                     KS.ib.wR, KS.ib.primR, KS.ib.h0R, KS.ib.h1R, 
                                                     KS.ib.h2R, KS.ib.h3R, KS.ib.ER, KS.ib.BR, KS.ib.lorenzR )
                    end
                end
            end
        
            for i=1:KS.pSpace.nx+1
                face[i] = MInterface1D4F(KS.ib.h0L)
            end

        else
        end
    
    end

    return ctr, face

end
