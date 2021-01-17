# Data Structure

Kinetic is organized with the data structures and methods of both generality and convenience. 
While most of the methods can be applied to multi-dimensional arrays directly, we provide a set of domain-specific structs that handles multiple dispatch in an elegant way.
In the finite volume method, the data is stored separately throughout the cells.
Therefore, we provide AbstractControlVolume structs for solving different equations that are used as arrays of structs (AoS) in the simulations.
The 1D control volume structs are
```@docs
ControlVolume1D
ControlVolume1D1F
ControlVolume1D2F
ControlVolume1D3F
ControlVolume1D4F
```

The 2D control volume structs are
```@docs
ControlVolume2D
ControlVolume2D1F
ControlVolume2D2F
ControlVolume2D3F
```

The numerical fluxes are evaluated through AbstractInterface structs
```@docs
Interface1D
Interface1D1F
Interface1D2F
Interface1D3F
Interface1D4F
Interface2D
Interface2D1F
Interface2D2F
```

The rest structs for saving general computational setups are
```@autodocs
Modules = [KitBase]
Order = [:type]
```