# Initial and boundary conditions

A struct `ib <: AbstractCondition` defines the initial and boundary conditions of a simulation.
It contains the values of conservative and primitive variables, and particle distribution functions at left and right (up and down) domain for both initial and boundary conditions.
It is set this way to easily deal with discontinuous initial conditions.

```@docs
IB
IB1F
IB2F
IB3F
IB4F
```