# Configuration of solver

Kinetic is organized with the data structures and methods of both generality and convenience. 
While the low-level methods can be applied to multi-dimensional arrays directly, we provide a set of domain-specific structs that handles multiple dispatch in an elegant way.

For a solver pending for execution, its configurations can be handled in a `SolverSet <: AbstractSolverSet` struct.
```@docs
SolverSet
```
It contains six fields:
- set: general setup of a simulation
- pSpace: physical space settings
- vSpace: particle velocity space settings
- gas: properties of the simulated substance
- ib: initial and boundary conditions
- outputFolder: file directory for the output results

This struct plays an key role in the solution algorithm.