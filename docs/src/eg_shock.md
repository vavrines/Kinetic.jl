# Tutorial 

Thanks to the brilliant expressiveness and low-overhead abstraction, we provide different levels of solution algorithm for modeling and simulating advection-diffusion dynamics.
The high-level solver is able to solve complex physics in a few lines, while the low-level APIs keep all the detailed implementations and are easy to be called from Python and C.
In the following, we present two quick tutorials to illustrate the usage of Kinetic.

## Shock tube problem

We first solve the shock tube problem in gas dynamics.
It's a two dimensional problem, with one in physical domain ``x`` and another in particle velocity domain ``u``.
First let us prepare the configuration file as
```
# case
case = sod
space = 1d2f1v
nSpecies = 1
flux = kfvs
collision = bgk
interpOrder = 2
limiter = vanleer
boundary = fix
cfl = 0.5
maxTime = 0.2

# physical space
x0 = 0
x1 = 1
nx = 200
pMeshType = uniform
nxg = 1

# velocity space
vMeshType = rectangle
umin = -5
umax = 5
nu = 28
nug = 0

# gas
knudsen = 0.0001
mach = 0.0
prandtl = 1
inK = 2
omega = 0.81
alphaRef = 1.0
omegaRef = 0.5
```

The configuration file can be understood as follows:
- The simulation case is the standard Sod shock tube
- A phase space in 1D physical and 1D velocity space is created with two particle distribution functions inside
- The numerical flux function is the kinetic flux vector splitting method and the collision term uses the BGK relaxation
- The reconstruction step employs van Leer limiter to create 2nd-order interpolation
- The two boundaries are fixed with Dirichlet boundary condition
- The timestep is determined with a CFL number of 0.5
- The maximum simulation time is 0.2
- The physical space spans in [0, 1] with 200 uniform cells
- The velocity space spans in [-5, 5] with 28 uniform cells
- The reference Knudsen number is set as 1e-4
- The reference Mach number is absent
- The reference Prandtl number is 1
- The gas molecule contains two internal degrees of freedom
- The viscosity is evaluated with the following formulas
```math
\mu = \mu_{ref} \left(\frac{T}{T_{ref}}\right)^{\omega}
```
```math
\mu_{ref}=\frac{5(\alpha+1)(\alpha+2) \sqrt{\pi}}{4 \alpha(5-2 \omega)(7-2 \omega)} Kn_{ref}
```

The configuration file directly generate variables during runtime via meta-programming in Julia,
and it can be stored in any text format (txt, toml, cfg, etc.). 
For example, if `config.txt` is created, 
we then execute the following codes to conduct a simulation
```julia
using Kinetic
set, ctr, face, t = initialize("config.txt")
t = solve!(set, ctr, face, t)
```

The computational setup is stored in `set` and the control volume solutions are stored in `ctr` and `face`. 
The high-level solver `solve!` is equivalent as the following low-level procedures
```julia
dt = timestep(ks, ctr, t)
nt = Int(floor(ks.set.maxTime / dt))
res = zeros(3)
for iter = 1:nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, res)
end
```

The result can be visualized with built-in function `plot_line`, which presents the profiles of gas density, velocity and temperature inside the tube.
```julia
plot_line(set, ctr)
```
![](./assets/sod.png)


## Lid-driven cavity

We then show the lid-driven cavity.
It's a four dimensional problem, with two in physical domain ``(x,y)`` and another in particle velocity domain ``(u,v)``.
Similarly, we prepare the configuration file as
```
# setup
case = cavity
space = 2d2f2v
flux = kfvs
collision = bgk
nSpecies = 1
interpOrder = 2
limiter = vanleer
boundary = maxwell
cfl = 0.8
maxTime = 2.0

# phase space
x0 = 0.0
x1 = 1.0
nx = 45
y0 = 0.0
y1 = 1.0
ny = 45
pMeshType = uniform
nxg = 0
nyg = 0

# velocity space
umin = -5.0
umax = 5.0
nu = 28
vmin = -5.0
vmax = 5.0
nv = 28
vMeshType = rectangle
nug = 0
nvg = 0

# gas
knudsen = 0.075
mach = 2.0
prandtl = 1.0
inK = 1.0
omega = 0.72
alphaRef = 1.0
omegaRef = 0.5

# boundary
uLid = 0.15
vLid = 0.0
tLid = 1.0
```

We then execute the following codes to conduct a simulation
```julia
using Kinetic
set, ctr, a1face, a2face, t = initialize("config.txt")
t = solve!(set, ctr, a1face, a2face, t)
```

The high-level solver `solve!` is equivalent as the following low-level procedures
```julia
using ProgressMeter
res = zeros(4)
dt = timestep(ks, ctr, simTime)
nt = floor(ks.set.maxTime / dt) |> Int
@showprogress for iter = 1:nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, a1face, a2face, dt; mode = Symbol(ks.set.flux), bc = Symbol(ks.set.boundary))
    update!(ks, ctr, a1face, a2face, dt, res; coll = Symbol(ks.set.collision), bc = Symbol(ks.set.boundary))
end
```

It can be further expanded into the lower-level backend.
```julia
# lower-level backend 
@showprogress for iter = 1:nt
    # horizontal flux
    @inbounds Threads.@threads for j = 1:ks.pSpace.ny
        for i = 2:ks.pSpace.nx
            KitBase.flux_kfvs!(
                a1face[i, j].fw,
                a1face[i, j].fh,
                a1face[i, j].fb,
                ctr[i-1, j].h,
                ctr[i-1, j].b,
                ctr[i, j].h,
                ctr[i, j].b,
                ks.vSpace.u,
                ks.vSpace.v,
                ks.vSpace.weights,
                dt,
                a1face[i, j].len,
            )
        end
    end
    
    # vertical flux
    vn = ks.vSpace.v
    vt = -ks.vSpace.u
    @inbounds Threads.@threads for j = 2:ks.pSpace.ny
        for i = 1:ks.pSpace.nx
            KitBase.flux_kfvs!(
                a2face[i, j].fw,
                a2face[i, j].fh,
                a2face[i, j].fb,
                ctr[i, j-1].h,
                ctr[i, j-1].b,
                ctr[i, j].h,
                ctr[i, j].b,
                vn,
                vt,
                ks.vSpace.weights,
                dt,
                a2face[i, j].len,
            )
            a2face[i, j].fw .= KitBase.global_frame(a2face[i, j].fw, 0., 1.)
        end
    end
    
    # boundary flux
    @inbounds Threads.@threads for j = 1:ks.pSpace.ny
        KitBase.flux_boundary_maxwell!(
            a1face[1, j].fw,
            a1face[1, j].fh,
            a1face[1, j].fb,
            ks.ib.bcL,
            ctr[1, j].h,
            ctr[1, j].b,
            ks.vSpace.u,
            ks.vSpace.v,
            ks.vSpace.weights,
            ks.gas.K,
            dt,
            ctr[1, j].dy,
            1.,
        )

        KitBase.flux_boundary_maxwell!(
            a1face[ks.pSpace.nx+1, j].fw,
            a1face[ks.pSpace.nx+1, j].fh,
            a1face[ks.pSpace.nx+1, j].fb,
            ks.ib.bcR,
            ctr[ks.pSpace.nx, j].h,
            ctr[ks.pSpace.nx, j].b,
            ks.vSpace.u,
            ks.vSpace.v,
            ks.vSpace.weights,
            ks.gas.K,
            dt,
            ctr[ks.pSpace.nx, j].dy,
            -1.,
        )
    end
    
    @inbounds Threads.@threads for i = 1:ks.pSpace.nx
        KitBase.flux_boundary_maxwell!(
            a2face[i, 1].fw,
            a2face[i, 1].fh,
            a2face[i, 1].fb,
            ks.ib.bcD,
            ctr[i, 1].h,
            ctr[i, 1].b,
            vn,
            vt,
            ks.vSpace.weights,
            ks.gas.K,
            dt,
            ctr[i, 1].dx,
            1,
        )
        a2face[i, 1].fw .= KitBase.global_frame(a2face[i, 1].fw, 0., 1.)
        
        KitBase.flux_boundary_maxwell!(
            a2face[i, ks.pSpace.ny+1].fw,
            a2face[i, ks.pSpace.ny+1].fh,
            a2face[i, ks.pSpace.ny+1].fb,
            [1., 0.0, -0.15, 1.0],
            ctr[i, ks.pSpace.ny].h,
            ctr[i, ks.pSpace.ny].b,
            vn,
            vt,
            ks.vSpace.weights,
            ks.gas.K,
            dt,
            ctr[i, ks.pSpace.ny].dy,
            -1,
        )
        a2face[i, ks.pSpace.ny+1].fw .= KitBase.global_frame(
            a2face[i, ks.pSpace.ny+1].fw,
            0.,
            1.,
        )
    end

    # update
    @inbounds for j = 1:ks.pSpace.ny
        for i = 1:ks.pSpace.nx
            KitBase.step!(
                ctr[i, j].w,
                ctr[i, j].prim,
                ctr[i, j].h,
                ctr[i, j].b,
                a1face[i, j].fw,
                a1face[i, j].fh,
                a1face[i, j].fb,
                a1face[i+1, j].fw,
                a1face[i+1, j].fh,
                a1face[i+1, j].fb,
                a2face[i, j].fw,
                a2face[i, j].fh,
                a2face[i, j].fb,
                a2face[i, j+1].fw,
                a2face[i, j+1].fh,
                a2face[i, j+1].fb,
                ks.vSpace.u,
                ks.vSpace.v,
                ks.vSpace.weights,
                ks.gas.K,
                ks.gas.γ,
                ks.gas.μᵣ,
                ks.gas.ω,
                ks.gas.Pr,
                ctr[i, j].dx * ctr[i, j].dy,
                dt,
                zeros(4),
                zeros(4),
                :bgk,
            )
        end
    end
end
```

The result can be visualized with built-in function `plot_contour`, which presents the contours of gas density, U-velocity, V-velocity and temperature inside the cavity.
```julia
# visulization
KitBase.plot_contour(ks, ctr)
```

![](./assets/cavity.png)

It is equivalent as the following backend.
```julia
# low-level backend
begin
    using Plots
    sol = zeros(4, ks.pSpace.nx, ks.pSpace.ny)
    for i in axes(ρ, 1)
        for j in axes(ρ, 2)
            sol[1:3, i, j] .= ctr[i, j].prim[1:3]
            sol[4, i, j] = 1.0 / ctr[i, j].prim[4]
        end
    end
    contourf(ks.pSpace.x[1:ks.pSpace.nx, 1], ks.pSpace.y[1, 1:ks.pSpace.ny], sol[3, :, :]')
end
```