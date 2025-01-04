using KitBase, Plots
using KitBase.ProgressMeter: @showprogress

set = Setup(; case="shock", space="1d1f3v", collision="fsm", maxTime=10)
ps = PSpace1D(-25, 25, 80, 1)
vs = VSpace3D(-10, 10, 64, -10, 10, 28, -10, 10, 28)
gas = begin
    _g = Gas(; Kn=1.0, Ma=3.0, K=0.0)
    Gas(; Kn=_g.Kn, Ma=_g.Ma, K=_g.K, fsm=fsm_kernel(vs, _g.μᵣ))
end
ib = IB1F(KB.config_ib(set, ps, vs, gas)...)

ks = SolverSet(set, ps, vs, gas, ib)
ctr, face = init_fvm(ks, :dynamic_array)

t = 0.0
dt = timestep(ks, ctr, t)
nt = ks.set.maxTime ÷ dt |> Int
res = zeros(5)

@showprogress for iter in 1:nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, res)
end

plot(ks, ctr)

