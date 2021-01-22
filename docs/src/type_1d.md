# One-dimensional data

In the finite volume method, the data is stored separately throughout the cells.
Therefore, we provide `AbstractControlVolume` and `AbstractInterface` structs for processing in-cell and edge information,
which are used as arrays of structs (AoS) in numerical simulations.
Considering one-dimensional physical space ``x``, we provide the following control volume structs.

```@docs
ControlVolume1D
ControlVolume1D1F
ControlVolume1D2F
ControlVolume1D3F
ControlVolume1D4F
```

Within each cell, different numbers of particle distribution function can be tracked.
Correspondingly, the interface data is stored with

```@docs
Interface1D
Interface1D1F
Interface1D2F
Interface1D3F
Interface1D4F
```