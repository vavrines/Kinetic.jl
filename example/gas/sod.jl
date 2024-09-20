using KitBase, Plots

cd(@__DIR__)
ks, ctr, face, t = initialize("sod.txt")
t = solve!(ks, ctr, face, t)

plot(ks, ctr)

#--- low-level procedures ---#
using KitBase.ProgressMeter: @showprogress

set = Setup(case = "sod", space = "1d2f1v", maxTime = 0.2)
ps = PSpace1D(0.0, 1.0, 100, 1)
vs = VSpace1D(-5.0, 5.0, 72)
gas = Gas(Kn = 1e-4, K = 2.0, γ = 5 / 3)
fw, ff, bc, p = KB.config_ib(set, ps, vs, gas)
ib = IB2F(fw, ff, bc, p)

ks = SolverSet(set, ps, vs, gas, ib)
ctr, face = init_fvm(ks, :static_array; structarray = true)

t = 0.0
dt = KB.timestep(ks, ctr, t)
nt = ks.set.maxTime ÷ dt |> Int
res = zeros(3)

@showprogress for iter = 1:nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, res)
end

plot(ks, ctr)

#--- pure 1D setting ---#
set = Setup(case = "sod", space = "1d1f1v", maxTime = 0.12)
gas = Gas(Kn = 1e-4, K = 0.0, γ = 3.0)
fw, ff, bc, p = config_ib(set, ps, vs, gas)
ib = IB1F(fw, ff, bc, p)

ks = SolverSet(set, ps, vs, gas, ib)
ctr, face = init_fvm(ks, :static_array; structarray = true)

t = 0.0
dt = KitBase.timestep(ks, ctr, t)
nt = ks.set.maxTime ÷ dt |> Int
res = zeros(3)
@showprogress for iter = 1:nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, res)
end

plot(ks, ctr)
