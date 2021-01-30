# Evolution

```@docs
evolve!
```

The evolution solver calculate the interface numerical fluxes based on two neighbor cells.
Different flux functions can be used with the option `model`.
- macroscopic: Godunov, Lax, Roe, HLL, wave-propagation
- mesoscopic: upwind, central-upwind, gas-kinetic scheme

The available flux solvers are listed as follows.
```@docs
flux_lax!
flux_hll!
flux_roe!
flux_gks
flux_gks!
flux_kfvs!
flux_kcu!
flux_ugks!
flux_boundary_maxwell!
flux_boundary_specular!
flux_em!
flux_emx!
flux_emy!
```
