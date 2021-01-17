# Data Structure

Kinetic is organized with the data structures and methods of both generality and convenience. 
While most of the methods can be applied to multi-dimensional arrays directly, we provide a set of domain-specific structs that handles multiple dispatch in an elegant way.
In the finite volume method, the data is stored separately throughout the cells.
Therefore, we provide AbstractControlVolume structs for solving different equations that are used as arrays of structs (AoS) in the simulations.
The 1D control volume structs are
```@docs
KitBase.ControlVolume1D
KitBase.ControlVolume1D1F
KitBase.ControlVolume1D2F
KitBase.ControlVolume1D3F
KitBase.ControlVolume1D4F
```

The 2D control volume structs are
```@docs
KitBase.ControlVolume2D
KitBase.ControlVolume2D1F
KitBase.ControlVolume2D2F
KitBase.ControlVolume2D3F
```

The numerical fluxes are evaluated through AbstractInterface structs
```@docs
KitBase.Interface1D
KitBase.Interface1D1F
KitBase.Interface1D2F
KitBase.Interface1D3F
KitBase.Interface1D4F
KitBase.Interface2D
KitBase.Interface2D1F
KitBase.Interface2D2F
```

The rest structs for saving general computational setups are
```@autodocs
Modules = [KitBase]
Order = [:type]
```