using KitBase, Plots

set = Setup(;
    matter="scalar", # material
    case="burgers", # test case
    space="1d0f0v", # phase space
    flux="gks", # flux
    collision="", # collision: for scalar conservation laws there are none
    interpOrder=1, # interpolation order
    boundary="period", # boundary
    cfl=0.5, # cfl
    maxTime=1.0, # simulation time
)
pSpace = PSpace1D(0.0, 1.0, 100, 1)
vSpace = nothing
property = Scalar(0.0, 1e-6)
ib = IB((x, p) -> sin(2π * x), property)

ks = SolverSet(set, pSpace, vSpace, property, ib)
ctr, face = init_fvm(ks)

t = 0.0
dt = timestep(ks, ctr, t)
nt = ks.set.maxTime ÷ dt |> Int
anim = @animate for iter in 1:nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, 0.0)

    plot(ks, ctr; xlabel="x", label="u", ylims=[-1, 1])
end

gif(anim, "burgers.gif"; fps=45)
