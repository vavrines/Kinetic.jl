# Reconstruction

```@docs
KitBase.reconstruct!
```

The reconstruction solver interpolates piecewise solutions with the desirable order of accuracy.
The reconstruction stencils can be based on 2 or 3 cells
```@docs
KitBase.reconstruct2
KitBase.reconstruct2!
KitBase.reconstruct3
KitBase.reconstruct3!
```

The available schemes are
```@docs
KitBase.vanleer
KitBase.minmod
KitBase.superbee
KitBase.vanalbaba
KitBase.weno5
```
