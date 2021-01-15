# Basics

The real world is built upon physics, which show a diverse set of behaviors on different characteristic scales.
Consider the molecular motion of gases as an example.
Down to finest scale of characteristic motion, the Newton's second law depicts particle motions
```math
\mathbf{F} = m \mathbf{a}
```
As a first order system this reads
```math
\frac{d \mathbf x}{dt} = \mathbf v, \frac{d \mathbf v}{dt} = \frac{\mathbf F}{m}
```

There are more than `2e25` molecules per cubic meter in normal atmosphere.
As a result, it is despairing to track all the individual motions every moment.
Instead, ensemble averaging is needed to provide a bottom-up view, from which the coarse-grained modeling happens.
At the mean free path and collision time scale of molecules, particles travel freely during most of time with mild intermolecular collisions.
Such dynamics can be described with an operator splitting approach, i.e. the Boltzmann equation

```math
\frac{\partial f}{\partial t}+ \mathbf u \cdot \nabla_\mathbf x f + \mathbf a \cdot \nabla_\mathbf u f = \int_{\mathcal R^3} \int_{\mathcal S^2} B (f(\mathbf u_*')f(\mathbf u')-f(\mathbf u_*) f(\mathbf u)) d\Omega d\mathbf u_*
```

Meanwhile, with the enlargement of modeling scale to a macroscopic hydrodynamic level, the accumulating effect of particle collisions results in an equalization of local temperature and velocity,
where the moderate non-equilibrium effects can be well described by viscous transport, heat conduction and mass diffusion,
i.e., the so called transport phenomena. 
From microscopic particle transport to macroscopic fluid motion, there is a
continuous variation of flow dynamics. 
In the comtinuum limit, the Euler and Navier-Stokes equations can be routinely used to describe macroscopic fluid evolutions.

```math
\frac{\partial \mathbf W}{\partial t} + \nabla_\mathbf x \cdot \mathbf F = \mathbf S
```