# General framework

Kinetic employs the finite volume method (FVM) for modeling and simulation. 
The general solution algorithm can be conclude as follows, where both explicit and implicit methods are implemented.

![](./assets/solver_process.png)

The high-level solver function is 
```@docs
solve!
```

The detailed solution procedures can be concluded as follows
- pre-process
- timestep calculation
- reconstruction
- evolution
- update
- post-process
