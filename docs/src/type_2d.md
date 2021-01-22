# Two-dimensional data

The 2D control volume structs are
```@docs
ControlVolume2D
ControlVolume2D1F
ControlVolume2D2F
ControlVolume2D3F
```

The numerical fluxes are evaluated through `AbstractInterface` structs
```@docs
Interface2D
Interface2D1F
Interface2D2F
```

The rest structs for saving general computational setups are
```@docs
SolverSet
```

```@autodocs
Modules = [KitBase]
Order = [:type]
```