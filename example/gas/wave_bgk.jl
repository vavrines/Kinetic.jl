using KitBase, Plots
using KitBase.ProgressMeter: @showprogress

set =
    Setup(; space="1d2f1v", collision="bgk", interpOrder=1, boundary="period", maxTime=1.0)
ps = PSpace1D(0.0, 1.0, 100, 1)
vs = VSpace1D(-5, 5, 60)

knudsen = 1e-2
μᵣ = ref_vhs_vis(knudsen, 1.0, 0.5)
gas = Gas(; Kn=knudsen, K=2.0)

begin
    fw = function (x, p)
        ρ = 1 + 0.1 * sin(2.0 * π * x)
        u = 1.0
        λ = ρ
        return prim_conserve([ρ, u, λ], gas.γ)
    end

    ff = function (x, p)
        w = fw(x, p)
        prim = conserve_prim(w, gas.γ)
        h = maxwellian(p.u, prim)
        b = h * p.K / 2 / prim[end]
        return h, b
    end

    bc = function (x, p)
        w = fw(x, p)
        return conserve_prim(w, gas.γ)
    end

    p = (u=vs.u, K=gas.K)
    ib = IB2F(fw, ff, bc, p)
end

ks = SolverSet(set, ps, vs, gas, ib)
ctr, face = init_fvm(ks)
plot(ks, ctr)

t = 0.0
dt = KitBase.timestep(ks, ctr, t)
nt = ks.set.maxTime ÷ dt |> Int
res = zeros(3)
@showprogress for iter in 1:nt
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, res)
end

plot!(ks, ctr; line=:dash)
