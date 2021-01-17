# Evolution

```@docs
KitBase.evolve!
```

The evolution solver calculate the interface numerical fluxes based on two neighbor cells.
- macroscopic: Godunov, Lax, Roe, HLL, wave-propagation
- mesoscopic: upwind, central-upwind, gas-kinetic scheme

The available flux solvers are
```@docs
KitBase.flux_lax!
KitBase.flux_hll!
KitBase.flux_roe!
KitBase.flux_gks
KitBase.flux_gks!
KitBase.flux_kfvs!
KitBase.flux_kcu!
KitBase.flux_ugks!
KitBase.flux_boundary_maxwell!
KitBase.flux_boundary_specular!
KitBase.flux_em!
KitBase.flux_emx!
KitBase.flux_emy!
```
