# Finite volume data

In the finite volume method, the data is stored separately throughout the cells.
Therefore, we provide `AbstractControlVolume` and `AbstractInterface` structs for processing in-cell and edge information,
which are used as arrays of structs (AoS) in numerical simulations.
Considering one-dimensional physical space ``x``, we provide the following control volume structs.
The structs differs from the number of particle distribution functions.
```@docs
ControlVolume1D
ControlVolume1D1F
ControlVolume1D2F
ControlVolume1D3F
ControlVolume1D4F
```

Within each cell, different numbers of particle distribution function can be tracked.
The interface data is stored correspondingly.
```@docs
Interface1D
Interface1D1F
Interface1D2F
Interface1D3F
Interface1D4F
```

The 2D control volume structs are implemented as well.
```@docs
ControlVolume2D
ControlVolume2D1F
ControlVolume2D2F
ControlVolume2D3F
```

The numerical fluxes are evaluated through `AbstractInterface` structs.
```@docs
Interface2D
Interface2D1F
Interface2D2F
```
