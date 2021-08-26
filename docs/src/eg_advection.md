# Advection diffusion

The first example is the scalar advection-diffusion equation.
It's a one dimensional problem in spatial domain ``x``.
Let's first configure the solver setup.
```julia
using Kinetic, Plots

set = Setup(
    matter = "scalar", # material
    case = "advection", # test case
    space = "1d0f0v", # phase space
    flux = "gks", # flux
    collision = "", # collision: for scalar conservation laws there are none
    interpOrder = 1, # interpolation order
    boundary = "period", # boundary condition
    cfl = 0.5, # cfl
    maxTime = 1.0, # simulation time
)
```

Then we generate the computational mesh.
Since we solve the macroscopic transport equation, the phase space is set to be nothing.
```julia
ps = PSpace1D(0.0, 1.0, 100, 1)
vs = nothing
```

And we define the physical property of material.
For the advection-diffusion equation, the two fields are the advection speed and viscosity respectively.
```julia
property = Scalar(1.0, 1e-6)
```

A sine wave is used as the initial condition.
```julia
ib = IB(x -> sin(2π * x), property)
```

For brevity, the above setups can be integrated into a single structure.
We also allocate the structures for cell-centered solutions and interface fluxes.
```julia
ks = SolverSet(set, ps, vs, property, ib)
ctr, face = init_fvm(ks)
```

The solution algorithm can be processed together with visualization.
```julia
t = 0.0
dt = KitBase.timestep(ks, ctr, t)
nt = ks.set.maxTime ÷ dt |> Int

anim = @animate for iter = 1:nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, 0.0)

    plot(ks, ctr, xlabel="x", label="u", ylims=[-1, 1])
end

gif(anim, "advection.gif", fps = 45)
```

![](./assets/advection.gif)
