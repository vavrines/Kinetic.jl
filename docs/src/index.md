# Home

Kinetic.jl is a lightweight toolbox for modeling and simulation on the basis of kinetic theory of gases, photons, plasmas and neutrons. 
It can be used to solve either Boltzmann and related model equations

```math
\frac{\partial f}{\partial t}+ \mathbf u \cdot \nabla_\mathbf x f + \mathbf a \cdot \nabla_\mathbf u f = \int_{\mathcal R^3} \int_{\mathcal S^2} B (f(\mathbf u_*')f(\mathbf u')-f(\mathbf u_*) f(\mathbf u)) d\Omega d\mathbf u_*
```

or their upscaling moment system

```math
\frac{\partial \mathbf W}{\partial t} + \nabla_\mathbf x \cdot \mathbf F = \mathbf S
```

The package is registered in the Julialang official entry and compatible with Julia 1.3 or newer version. 
To make use of it, execute `Julia` and type
```julia
julia> ]
(v1.3) pkg> add Kinetic
```
This will install Kinetic and all its dependencies.
After that, load the package,
```julia
julia> using Kinetic
```