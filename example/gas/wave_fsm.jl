using KitBase, Plots
using KitBase.ProgressMeter: @showprogress

set =
    Setup(; space="1d1f3v", collision="fsm", interpOrder=1, boundary="period", maxTime=1.0)
ps = PSpace1D(0.0, 1.0, 100, 1)
vs = VSpace3D(-8, 8, 48, -8, 8, 28, -8, 8, 28)

knudsen = 1e-2
μᵣ = ref_vhs_vis(knudsen, 1.0, 0.5)
fsm = fsm_kernel(vs, μᵣ, 5, 1.0)
gas = Gas(; Kn=knudsen, K=0.0, fsm=fsm)

begin
    fw = function (x, p)
        ρ = 1 + 0.1 * sin(2 * π * x)
        u = 1.0
        λ = ρ
        return prim_conserve([ρ, u, 0, 0, λ], gas.γ)
    end

    ff = function (x, p)
        w = fw(x, p)
        prim = conserve_prim(w, gas.γ)
        return maxwellian(vs.u, vs.v, vs.w, prim)
    end

    bc = function (x, p)
        w = fw(x, p)
        return conserve_prim(w, gas.γ)
    end

    p = NamedTuple()
    ib = IB1F(fw, ff, bc, p)
end

ks = SolverSet(set, ps, vs, gas, ib)
ctr, face = init_fvm(ks)
plot(ks, ctr)

t = 0.0
dt = KitBase.timestep(ks, ctr, t)
nt = ks.set.maxTime ÷ dt |> Int
res = zeros(3)
@showprogress for iter in 1:100#nt
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, res)
end

plot!(ks, ctr; line=:dash)
