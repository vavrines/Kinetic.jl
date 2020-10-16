# Installation Instructions

Kinetic.jl is a registered Julialang package in the official entry.
You can install the latest version using the built-in package manager (accessed by pressing `]` in the `Julia` command prompt) to add the package and instantiate/build all dependencies.

```julia
julia> ]
(v1.5) pkg> add Kinetic
(v1.5) pkg> build Kinetic
```
This will install Kinetic and all its dependencies.
After that, load the package,
```julia
julia> using Kinetic
```

Kinetic.jl can be updated to the latest tagged release from the package manager by typing
```julia
(v1.5) pkg> update Kinetic
```
