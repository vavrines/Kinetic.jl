# General

The finite volume method (FVM) is employed in Kinetic. 
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
