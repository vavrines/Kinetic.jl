@everywhere using Distributed, SharedArrays
addprocs(3)
@everywhere using KitBase

begin
    vars = Dict{String,Any}()
    vars["case"] = "sod"
    vars["space"] = "1d0f0v"
    vars["flux"] = "kfvs"
    vars["collision"] = "bgk"
    vars["nSpecies"] = 1
    vars["interpOrder"] = 1
    vars["limiter"] = "vanleer"
    vars["boundary"] = "fix"
    vars["cfl"] = 0.5
    vars["maxTime"] = 0.2
    vars["x0"] = 0.0
    vars["x1"] = 1.0
    vars["nx"] = 2000
    vars["pMeshType"] = "uniform"
    vars["nxg"] = 0
    vars["knudsen"] = 0.001
    vars["mach"] = 0.0
    vars["prandtl"] = 1.0
    vars["inK"] = 0.0
    vars["omega"] = 0.81
    vars["alphaRef"] = 1.0
    vars["omegaRef"] = 0.5
end

set = KitBase.set_setup(vars)
pSpace = KitBase.set_geometry(vars)
vSpace = KitBase.set_velocity(vars)
gas = KitBase.set_property(vars)
ib = KitBase.set_ib(vars, set, vSpace, gas)
folder = @__DIR__
ks = KitBase.SolverSet(set, pSpace, vSpace, gas, ib, folder)

dt = ks.pSpace.dx[1] / (5.0 + KitBase.sound_speed(ks.ib.primL, ks.gas.γ))
nt = floor(ks.set.maxTime / dt) |> Int

#--- parallel ---#
wp = SharedArray{Float64}((ks.pSpace.nx, 3), init = A -> (A = zeros(ks.pSpace.nx, 3)))
for i = 1:ks.pSpace.nx
    if i <= ks.pSpace.nx ÷ 2
        wp[i, :] .= ks.ib.wL
    else
        wp[i, :] .= ks.ib.wR
    end
end
fwp = SharedArray{Float64}(
    (ks.pSpace.nx + 1, 3),
    init = A -> (A = zeros(ks.pSpace.nx + 1, 3)),
)

@time for iter = 1:nt÷3
    @sync @distributed for i = 2:ks.pSpace.nx
        flux = @view fwp[i, :]
        KitBase.flux_gks!(
            flux,
            wp[i-1, :],
            wp[i, :],
            ks.gas.γ,
            ks.gas.K,
            ks.gas.μᵣ,
            ks.gas.ω,
            dt,
            0.5 * ks.pSpace.dx[i-1],
            0.5 * ks.pSpace.dx[i],
        )
    end

    @sync @distributed for i = 2:ks.pSpace.nx-1
        for j = 1:3
            wp[i, j] += (fwp[i, j] - fwp[i+1, j]) / ks.pSpace.dx[i]
        end
    end
end
"""~13.620491 seconds (2.26 M allocations: 101.219 MiB, 0.22% gc time)"""

#--- serial ---#
w = zeros(ks.pSpace.nx, 3)
for i = 1:ks.pSpace.nx
    if i <= ks.pSpace.nx ÷ 2
        w[i, :] .= ks.ib.wL
    else
        w[i, :] .= ks.ib.wR
    end
end
fw = zeros(ks.pSpace.nx + 1, 3)

@time for iter = 1:nt÷3
    for i = 2:ks.pSpace.nx
        flux = @view fw[i, :]
        KitBase.flux_gks!(
            flux,
            w[i-1, :],
            w[i, :],
            ks.gas.γ,
            ks.gas.K,
            ks.gas.μᵣ,
            ks.gas.ω,
            dt,
            0.5 * ks.pSpace.dx[i-1],
            0.5 * ks.pSpace.dx[i],
        )
    end

    for i = 2:ks.pSpace.nx-1
        for j = 1:3
            w[i, j] += (fw[i, j] - fw[i+1, j]) / ks.pSpace.dx[i]
        end
    end
end
"""~20.830331 seconds (323.96 M allocations: 24.472 GiB, 16.89% gc time)"""
