# Advection diffusion

The first example is the scalar advection-diffusion equation.
It's a one dimensional problem in physical domain ``x``.
Let's first configure the solver setup.
```julia
set = Setup(
    "scalar", # matter
    "advection", # case
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
property = Scalar(1.0, 1e-6)
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

gif(anim, "advection.gif", fps = 45)
```

![](./assets/advection.gif)
