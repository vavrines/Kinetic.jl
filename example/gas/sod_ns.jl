using KitBase, Plots
using KitBase.ProgressMeter: @showprogress

set =
    Setup(; case="sod", space="1d0f0v", flux="gks", limiter="minmod", cfl=0.3, maxTime=0.1)
ps = PSpace1D(0.0, 1.0, 100, 1)
vs = nothing
gas = Gas(; Kn=1e-4, K=2.0, γ=5 / 3)
fw, bc, p = KB.config_ib(set, ps, vs, gas)
ib = IB(fw, bc, p)

ks = SolverSet(set, ps, vs, gas, ib)
ctr, face = init_fvm(ks)

t = 0.0
dt = KB.timestep(ks, ctr, t)
nt = ks.set.maxTime ÷ dt |> Int
res = zeros(3)

@showprogress for iter in 1:nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, res)
end

plot(ks, ctr)

# low-level gks navier-stokes flux
@showprogress for iter in 1:nt
    reconstruct!(ks, ctr)
    for i in 1:ps.nx+1
        flux_gks!(
            face[i].fw,
            ctr[i-1].w .+ ps.dx[i-1] / 2 .* ctr[i-1].sw,
            ctr[i].w .- ps.dx[i] / 2 .* ctr[i].sw,
            gas.K,
            gas.γ,
            gas.μᵣ,
            gas.ω,
            dt,
            ps.dx[i-1] / 2,
            ps.dx[i] / 2,
            ctr[i-1].sw,
            ctr[i].sw,
        )
    end
    update!(ks, ctr, face, dt, res)
end

plot(ks, ctr)
