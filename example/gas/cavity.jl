using KitBase, Plots
using KitBase: ib_cavity
using KitBase.ProgressMeter: @showprogress

set = Setup(
    case = "cacity",
    space = "2d2f2v",
    boundary = ["maxwell", "maxwell", "maxwell", "maxwell"],
    limiter = "minmod",
    cfl = 0.5,
    maxTime = 10, # time
)
ps = PSpace2D(0, 1, 45, 0, 1, 45)
vs = VSpace2D(-5, 5, 28, -5, 5, 28)
gas = Gas(Kn = 0.075, K = 1.0)
ib = IB2F(ib_cavity(set, ps, vs, gas)...)

ks = SolverSet(set, ps, vs, gas, ib)
ctr, a1face, a2face = init_fvm(ks)

t = 0.0
dt = timestep(ks, ctr, 0.0)
nt = ks.set.maxTime ÷ dt |> Int
res = zeros(4)

@showprogress for iter = 1:nt
    evolve!(ks, ctr, a1face, a2face, dt)
    update!(ks, ctr, a1face, a2face, dt, res)

    if maximum(res) < 1e-6
        break
    end
end

plot(ks, ctr)
