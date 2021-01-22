# General

Kinetic is organized with the data structures and methods of both generality and convenience. 
While the low-level methods can be applied to multi-dimensional arrays directly, we provide a set of domain-specific structs that handles multiple dispatch in an elegant way.

The struct for setting up the solver is
```@docs
SolverSet
```

The inner structs of `AbstractSolverSet` are
```@autodocs
Modules = [KitBase]
Order = [:type]
```