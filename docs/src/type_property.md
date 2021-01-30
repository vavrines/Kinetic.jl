# Particle properties

A struct `gas <: AbstractProperty` defines the properties of particle model.
It currently supports the following models:
- scalar
- gas-type molecule
- plasma

```@docs
Scalar
Gas
Mixture
Plasma1D
Plasma2D
```

The fields denote, for example:
- Kn: reference Knudsen number
- Ma: reference Mach number
- Pr: reference Prandtl number
- K: internal degree of freedom of molecule
- γ: adiabatic index
- ω: viscosity index
- αᵣ: reference ``\alpha`` in viscosity evaluation
- ωᵣ: reference ``\omega`` in viscosity evaluation
- μᵣ: reference viscosity
- m: mass of each particle
- np: number of particles

The viscosity is evaluated the following hard-sphere model.
```math
\mu = \mu_{ref} \left(\frac{T}{T_{ref}}\right)^{\omega}
```
```math
\mu_{ref}=\frac{5(\alpha+1)(\alpha+2) \sqrt{\pi}}{4 \alpha(5-2 \omega)(7-2 \omega)} Kn_{ref}
```
