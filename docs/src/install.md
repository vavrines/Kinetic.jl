# Installation Instructions

Kinetic is a registered Julia package in the official entry.
We recommend installing it with the built-in Julia package manager, because this automatically installs a stable, tagged release. 
From the Julia REPL, you could add the package and instantiate/build all dependencies via
```julia
julia> ]
(v1.5) pkg> add Kinetic
(v1.5) pkg> build Kinetic
```
This will install Kinetic and all its dependencies.
After that, using/import the package,
```julia
julia> using Kinetic
```

Kinetic can be updated to the latest tagged release from the package manager by executing
```julia
(v1.5) pkg> update Kinetic
```

!!! tip "Use Julia 1.3 or newer"
    Kinetic matches perfectly with Julia 1.3 and newer versions.
    Installing it with an older version of Julia will locate incomplete functionality.
