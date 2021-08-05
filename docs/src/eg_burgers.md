# Burgers

Now we could turn to the Burgers equation.
It's a typical hyperbolic conservation law, where discontinuous solution can emerge in a self-evolving system.
Let's consider the same initial configuration as advection-diffusion example.
```julia
using KitBase, Plots

set = Setup(
    "scalar", # matter
    "burgers", # case
    "1d0f0v", # space
    "gks", # flux
    "", # collision: for scalar conservation laws there are none
    1, # species
    1, # interpolation order
    "vanleer", # limiter
    "period", # boundary
    0.5, # cfl
    1.0, # simulation time
)
pSpace = PSpace1D(0.0, 1.0, 100, 1)
vSpace = nothing
property = Scalar(0, 1e-4)
w0 = 1.0
prim0 = conserve_prim(w0, property.a)
ib = IB(w0, prim0, prim0, w0, prim0, prim0)
ks = SolverSet(set, pSpace, vSpace, property, ib)
```

The we allocate the data structure needed.
```julia
ctr, face = init_fvm(ks, ks.ps)
for i in eachindex(ctr)
    ctr[i].w = sin(2π * ks.pSpace.x[i]) # initial condition
    ctr[i].prim .= conserve_prim(ctr[i].w)
end
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

    sol = zeros(ks.pSpace.nx)
    for i in 1:ks.pSpace.nx
        sol[i] = ctr[i].w
    end
    plot(ks.pSpace.x[1:ks.pSpace.nx], sol, xlabel="x", label="u", ylims=[-1, 1])
end

gif(anim, "burgers.gif", fps = 45)
```

![](./assets/burgers.gif)
