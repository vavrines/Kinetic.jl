# Parameter Settings

A struct `set <: AbstractSetup` defines the general setup of a simulation.
```@docs
Setup
```
It contains
- matter: fluid substance
- case: simulation case name
- space: ``n_1 d n_2 f n_3 v``, which denotes the physical dimensionality, numbers of particle distribution functions and velocity dimensionality
- flux: numerical flux function name
- collision: collision operator of kinetic equation
- nSpecies: number of species
- interpOrder: order of accuracy for reconstruction
- limiter: limiter function name
- boundary: boundary condition
- cfl: Courant-Friedrichs-Lewy number for determining time step
- maxTime: maximum simulation time