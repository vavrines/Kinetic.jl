# Particle velocity space

A struct `vSpace <: AbstractSetup` defines the particle velocity setup of a simulation.
Structs for 1-3 dimensional particle velocity space are built.
```@docs
VSpace1D
VSpace2D
VSpace3D
```
It contains
- u0 (v0, w0): location of starting point
- u1 (v1, w1): location of ending point
- nu (nv, nw): number of cells in one direction
- u (v, w): locations of middle points of all cells
- du (dv, dw): intervals of all cell points
- weights: quadrature weights for numerical integral

Note that the one-dimensional velocity space can be used to handle 1-3 dimensional unstructured topology as well.
In addition, velocity space structs for multi-component substance are implemented.
```@docs
MVSpace1D
MVSpace2D
MVSpace3D
```

For the simulation cases where no phase-space evolution is involved, `vSpace` can be set as `nothing` directly.
```julia
vSpace = nothing
```