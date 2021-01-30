# Installation Instructions

Kinetic is a registered Julia package in the official entry.
We recommend installing it with the built-in Julia package manager.
It automatically installs a currently stable and tagged release. 
From the Julia REPL, you can add the package, build and instantiate the environment.
```julia
julia> ]
(v1.6) pkg> add Kinetic
(v1.6) pkg> build Kinetic
(v1.6) pkg> instantiate
```

This will automatically install Kinetic and all its dependencies.
After that, we can `using` or `import` the package.
`using` will load the module and make its exported names available for direct use.
```julia
julia> using Kinetic
julia> linspace(0, 1, 5)
5-element Vector{Float64}:
 0.0
 0.25
 0.5
 0.75
 1.0
```
Correspondingly, `import` only loads the module while the names needs to be accessed with dot syntax.
```julia
julia> import Kinetic
julia> Kinetic.linspace(0, 1, 5)
5-element Vector{Float64}:
 0.0
 0.25
 0.5
 0.75
 1.0
```

Kinetic can be updated to the latest tagged release from the package manager.
```julia
(v1.5) pkg> update Kinetic
```

!!! tip "Use Julia 1.3 or newer"
    Kinetic matches perfectly with Julia 1.3 and newer versions.
    Installing it with an older version of Julia will locate incomplete functionality.