# ============================================================
# Initializer of Simulation
# ============================================================


"""
Initialize solver from input file

"""
function initialize(configfilename::T) where {T<:AbstractString}

    println("==============================================================")
    println("Kinetic.jl")
    println("A Lightweight Toolbox for Kinetic Modeling and Simulation")
    println("==============================================================")
    println("")
    println("reading configurations from $(configfilename)")
    println("")
    println("initializing solver: ")

    if configfilename[end-2:end] == "txt"

        ks = SolverSet(configfilename)
        ctr, face = init_fvm(ks)

        return ks, ctr, face, 0.0

    elseif configfilename[end-3:end] == "jld2"

        _1, _2, _3 = @load configfilename KS ctr t
        ks, ctr, simTime = eval(_1), eval(_2), eval(_3)

        face = init_fvm(ks)[2]

        return ks, ctr, face, simTime

    end

end


"""
Initialize finite volume method

"""
function init_fvm(KS::T) where {T<:AbstractSolverSet}

    if KS.set.space[1:2] == "1d"

        if KS.set.space[3:4] == "1f"

            ctr = OffsetArray{ControlVolume1D1F}(undef, eachindex(KS.pSpace.x))
            face = Array{Interface1D1F}(undef, KS.pSpace.nx + 1)

            for i in eachindex(ctr)
                if i <= KS.pSpace.nx ÷ 2
                    ctr[i] = ControlVolume1D1F(
                        KS.pSpace.x[i],
                        KS.pSpace.dx[i],
                        KS.ib.wL,
                        KS.ib.primL,
                        KS.ib.fL,
                    )
                else
                    ctr[i] = ControlVolume1D1F(
                        KS.pSpace.x[i],
                        KS.pSpace.dx[i],
                        KS.ib.wR,
                        KS.ib.primR,
                        KS.ib.fR,
                    )
                end
            end

            for i = 1:KS.pSpace.nx+1
                face[i] = Interface1D1F(KS.ib.wL, KS.ib.fL)
            end

        elseif KS.set.space[3:4] == "2f"

            ctr = OffsetArray{ControlVolume1D2F}(undef, eachindex(KS.pSpace.x))
            face = Array{Interface1D2F}(undef, KS.pSpace.nx + 1)

            for i in eachindex(ctr)
                if i <= KS.pSpace.nx ÷ 2
                    ctr[i] = ControlVolume1D2F(
                        KS.pSpace.x[i],
                        KS.pSpace.dx[i],
                        KS.ib.wL,
                        KS.ib.primL,
                        KS.ib.hL,
                        KS.ib.bL,
                    )
                else
                    ctr[i] = ControlVolume1D2F(
                        KS.pSpace.x[i],
                        KS.pSpace.dx[i],
                        KS.ib.wR,
                        KS.ib.primR,
                        KS.ib.hR,
                        KS.ib.bR,
                    )
                end
            end

            for i = 1:KS.pSpace.nx+1
                face[i] = Interface1D2F(KS.ib.wL, KS.ib.hL)
            end

        elseif KS.set.space[3:4] == "3f"

            ctr = OffsetArray{ControlVolume1D3F}(undef, eachindex(KS.pSpace.x))
            face = Array{Interface1D3F}(undef, KS.pSpace.nx + 1)

            for i in eachindex(ctr)
                if i <= KS.pSpace.nx ÷ 2
                    ctr[i] = ControlVolume1D3F(
                        KS.pSpace.x[i],
                        KS.pSpace.dx[i],
                        KS.ib.wL,
                        KS.ib.primL,
                        KS.ib.h0L,
                        KS.ib.h1L,
                        KS.ib.h2L,
                        KS.ib.EL,
                        KS.ib.BL,
                        KS.ib.lorenzL,
                    )
                else
                    ctr[i] = ControlVolume1D3F(
                        KS.pSpace.x[i],
                        KS.pSpace.dx[i],
                        KS.ib.wR,
                        KS.ib.primR,
                        KS.ib.h0R,
                        KS.ib.h1R,
                        KS.ib.h2R,
                        KS.ib.ER,
                        KS.ib.BR,
                        KS.ib.lorenzR,
                    )
                end
            end

            for i = 1:KS.pSpace.nx+1
                face[i] = Interface1D3F(KS.ib.wL, KS.ib.h0L, KS.ib.EL)
            end

        elseif KS.set.space[3:4] == "4f"

            ctr = OffsetArray{ControlVolume1D4F}(undef, eachindex(KS.pSpace.x))
            face = Array{Interface1D4F}(undef, KS.pSpace.nx + 1)

            for i in eachindex(ctr)
                if i <= KS.pSpace.nx ÷ 2
                    ctr[i] = ControlVolume1D4F(
                        KS.pSpace.x[i],
                        KS.pSpace.dx[i],
                        KS.ib.wL,
                        KS.ib.primL,
                        KS.ib.h0L,
                        KS.ib.h1L,
                        KS.ib.h2L,
                        KS.ib.h3L,
                        KS.ib.EL,
                        KS.ib.BL,
                        KS.ib.lorenzL,
                    )
                else
                    ctr[i] = ControlVolume1D4F(
                        KS.pSpace.x[i],
                        KS.pSpace.dx[i],
                        KS.ib.wR,
                        KS.ib.primR,
                        KS.ib.h0R,
                        KS.ib.h1R,
                        KS.ib.h2R,
                        KS.ib.h3R,
                        KS.ib.ER,
                        KS.ib.BR,
                        KS.ib.lorenzR,
                    )
                end
            end

            for i = 1:KS.pSpace.nx+1
                face[i] = Interface1D4F(KS.ib.wL, KS.ib.h0L, KS.ib.EL)
            end

        end

        return ctr, face

    elseif KS.set.space[1:2] == "2d"

        if KS.set.space[3:4] == "1f"

            ctr = OffsetArray{ControlVolume2D1F}(
                undef,
                axes(KS.pSpace.x, 1),
                axes(KS.pSpace.y, 2),
            )
            a1face = Array{Interface2D1F}(undef, KS.pSpace.nx + 1, KS.pSpace.ny)
            a2face = Array{Interface2D1F}(undef, KS.pSpace.nx, KS.pSpace.ny + 1)

            for j in axes(ctr, 2), i in axes(ctr, 1)
                if i <= KS.pSpace.nx ÷ 2
                    ctr[i, j] = ControlVolume2D1F(
                        KS.pSpace.x[i, j],
                        KS.pSpace.y[i, j],
                        KS.pSpace.dx[i, j],
                        KS.pSpace.dy[i, j],
                        KS.ib.wL,
                        KS.ib.primL,
                        KS.ib.fL,
                    )
                else
                    ctr[i] = ControlVolume2D1F(
                        KS.pSpace.x[i, j],
                        KS.pSpace.y[i, j],
                        KS.pSpace.dx[i, j],
                        KS.pSpace.dy[i, j],
                        KS.ib.wR,
                        KS.ib.primR,
                        KS.ib.fR,
                    )
                end
            end

            for j = 1:KS.pSpace.ny, i = 1:KS.pSpace.nx+1
                a1face[i] = Interface2D1F(KS.pSpace.dy[i, j], 0.0, 1.0, KS.ib.wL, KS.ib.hL)
            end
            for j = 1:KS.pSpace.ny+1, i = 1:KS.pSpace.nx
                a2face[i] = Interface2D1F(KS.pSpace.dx[i, j], 1.0, 0.0, KS.ib.wL, KS.ib.hL)
            end

        elseif KS.set.space[3:4] == "2f"

            ctr = OffsetArray{ControlVolume2D2F}(
                undef,
                axes(KS.pSpace.x, 1),
                axes(KS.pSpace.y, 2),
            )
            a1face = Array{Interface2D2F}(undef, KS.pSpace.nx + 1, KS.pSpace.ny)
            a2face = Array{Interface2D2F}(undef, KS.pSpace.nx, KS.pSpace.ny + 1)

            for j in axes(ctr, 2), i in axes(ctr, 1)
                if i <= KS.pSpace.nx ÷ 2
                    ctr[i, j] = ControlVolume2D2F(
                        KS.pSpace.x[i, j],
                        KS.pSpace.y[i, j],
                        KS.pSpace.dx[i, j],
                        KS.pSpace.dy[i, j],
                        KS.ib.wL,
                        KS.ib.primL,
                        KS.ib.hL,
                        KS.ib.bL,
                    )
                else
                    ctr[i, j] = ControlVolume2D2F(
                        KS.pSpace.x[i, j],
                        KS.pSpace.y[i, j],
                        KS.pSpace.dx[i, j],
                        KS.pSpace.dy[i, j],
                        KS.ib.wR,
                        KS.ib.primR,
                        KS.ib.hR,
                        KS.ib.bR,
                    )
                end
            end

            for j = 1:KS.pSpace.ny, i = 1:KS.pSpace.nx+1
                a1face[i] = Interface2D2F(KS.pSpace.dy[i, j], 0.0, 1.0, KS.ib.wL, KS.ib.hL)
            end
            for j = 1:KS.pSpace.ny+1, i = 1:KS.pSpace.nx
                a2face[i] = Interface2D2F(KS.pSpace.dx[i, j], 1.0, 0.0, KS.ib.wL, KS.ib.hL)
            end

        end

        return ctr, a1face, a2face

    end

end
