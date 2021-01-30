# Update

The update solver calculate the variables at n+1 step based on numerical fluxes and in-cell collisions.
```@docs
update!
```

The current solver supports different collision models, for example:
- `:bgk`: BGK relaxation model
- `:shakhov`: Shakhov relaxation model
- `:boltzmann`: Boltzmann: original Boltzmann collision integral

The boundary conditions vary.
- `:fix`: fixed Dirichlet boundary
- `:period`: periodic boundary
- `:extra`: extrapolation
- `:maxwell`: Maxwell's diffusive boundary

The current solver adopts implicit-explicit (IMEX) uniformly.
Further Multi-step time integrators can be used in conjunction with method of lines in [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).
