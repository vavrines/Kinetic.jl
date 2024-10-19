using KitBase, LinearAlgebra
using KitBase.JLD2
using KitBase.ProgressMeter: @showprogress

cd(@__DIR__)
D = KitBase.read_dict("cylinder_unstr.txt")
set = KitBase.set_setup(D)

ps = KitBase.set_geometry(D)
for i in eachindex(ps.faceType)
    i1 = ps.facePoints[i, 1]
    i2 = ps.facePoints[i, 2]

    if i1 in [ps.cells.index[1]; ps.cells.index[2]] &&
       i2 in [ps.cells.index[1]; ps.cells.index[2]]
        ps.faceType[i] = 2
        c1 = ps.faceCells[i, 1]
        c2 = ps.faceCells[i, 2]

        if c1 != -1
            ps.cellType[c1] = 2
        elseif c2 != -1
            ps.cellType[c2] = 2
        else
            throw("index error")
        end
    end
end
for i in eachindex(ps.cellType)
    if ps.cellType[i] == 1 && ps.cellCenter[i, 1] > 0.0
        ps.cellType[i] = 3 # outlet boundary
    end
end

vs = KitBase.set_velocity(D)
gas = KitBase.set_property(D)

begin
    primL = [1.0, KitBase.sound_speed(1.0, gas.γ) * gas.Ma, 0.0, 1.0]
    wL = KitBase.prim_conserve(primL, gas.γ)
    hL = KitBase.maxwellian(vs.u, vs.v, primL)
    bL = @. hL * gas.K / 2 / primL[end]
    primR = [1.0, 0.0, 0.0, 1.0]
    wR = KitBase.prim_conserve(primR, gas.γ)
    hR = KitBase.maxwellian(vs.u, vs.v, primR)
    bR = @. hR * gas.K / 2 / primR[end]

    p = (wL=wL, wR=wR, primL=primL, primR=primR, HL=hL, HR=hR, BL=bL, BR=bR)

    fw = function (x, y, p)
        return p.wL
    end
    ff = function (x, y, p)
        return p.HL, p.BL
    end
    bc = function (x, y, p)
        return p.primR
    end

    ib = IB2F(fw, ff, bc, p)
end

ks = KitBase.SolverSet(set, ps, vs, gas, ib, @__DIR__)
ctr, face = KitBase.init_fvm(ks, ks.ps)
dt = KitBase.timestep(ks, ctr, 0.0)
nt = ks.set.maxTime ÷ dt |> Int

@showprogress for iter in 1:nt
    @inbounds Threads.@threads for i in eachindex(face)
        vn = ks.vs.u .* face[i].n[1] .+ ks.vs.v .* face[i].n[2]
        vt = ks.vs.v .* face[i].n[1] .- ks.vs.u .* face[i].n[2]

        if !(-1 in ps.faceCells[i, :])
            KitBase.flux_kfvs!(
                face[i].fw,
                face[i].fh,
                face[i].fb,
                ctr[ps.faceCells[i, 1]].h,
                ctr[ps.faceCells[i, 1]].b,
                ctr[ps.faceCells[i, 2]].h,
                ctr[ps.faceCells[i, 2]].b,
                vn,
                vt,
                ks.vs.weights,
                dt,
                face[i].len,
            )
            face[i].fw .= KitBase.global_frame(face[i].fw, face[i].n[1], face[i].n[2])
        else
            idx = ifelse(ps.faceCells[i, 1] != -1, 1, 2)

            if ps.cellType[ps.faceCells[i, idx]] == 2
                _prim = ks.ib.bc(ps.faceCenter[i, 1], ps.faceCenter[i, 2], ks.ib.p)
                bc = KB.local_frame(_prim, face[i].n[1], face[i].n[2])

                KitBase.flux_boundary_maxwell!(
                    face[i].fw,
                    face[i].fh,
                    face[i].fb,
                    bc,
                    ctr[ps.faceCells[i, idx]].h,
                    ctr[ps.faceCells[i, idx]].b,
                    vn,
                    vt,
                    ks.vs.weights,
                    ks.gas.K,
                    dt,
                    face[i].len,
                    1,
                )

                face[i].fw .= KitBase.global_frame(face[i].fw, face[i].n[1], face[i].n[2])
            end
        end
    end

    sumres = zeros(4)
    sumavg = zeros(4)
    @inbounds Threads.@threads for i in eachindex(ctr)
        if ps.cellType[i] in (0, 2)
            dirc = [sign(dot(ctr[i].n[j], face[ps.cellFaces[i, j]].n)) for j in 1:3]

            KitBase.step!(
                ctr[i].w,
                ctr[i].prim,
                ctr[i].h,
                ctr[i].b,
                face[ps.cellFaces[i, 1]].fw,
                face[ps.cellFaces[i, 1]].fh,
                face[ps.cellFaces[i, 1]].fb,
                face[ps.cellFaces[i, 2]].fw,
                face[ps.cellFaces[i, 2]].fh,
                face[ps.cellFaces[i, 2]].fb,
                face[ps.cellFaces[i, 3]].fw,
                face[ps.cellFaces[i, 3]].fh,
                face[ps.cellFaces[i, 3]].fb,
                ks.vs.u,
                ks.vs.v,
                ks.vs.weights,
                ks.gas.K,
                ks.gas.γ,
                ks.gas.μᵣ,
                ks.gas.ω,
                ks.gas.Pr,
                ks.ps.cellArea[i],
                dirc,
                dt,
                sumres,
                sumavg,
                :bgk,
            )
        end
    end

    for i in eachindex(ps.cellType)
        if ps.cellType[i] == 3
            ids = ps.cellNeighbors[i, :]
            deleteat!(ids, findall(x -> x == -1, ids))
            id1, id2 = ids
            ctr[i].w .= 0.5 .* (ctr[id1].w .+ ctr[id2].w)
            ctr[i].h .= 0.5 .* (ctr[id1].h .+ ctr[id2].h)
            ctr[i].b .= 0.5 .* (ctr[id1].b .+ ctr[id2].b)
            ctr[i].prim .= conserve_prim(ctr[i].w, ks.gas.γ)
        end
    end
end

KB.write_vtk(ks, ctr)
