# Kinetic.jl

Kinetic.jl is a lightweight toolbox for modeling and simulation on the basis of kinetic theory of gases, photons, plasmas and neutrons.
It can be used to solve either Boltzmann and its related model equations

```math
\frac{\partial f}{\partial t}+ \mathbf u \cdot \nabla_\mathbf x f + \mathbf a \cdot \nabla_\mathbf u f = \int_{\mathcal R^3} \int_{\mathcal S^2} B (f(\mathbf u_*')f(\mathbf u')-f(\mathbf u_*) f(\mathbf u)) d\Omega d\mathbf u_*
```

or their upscaling moment system

```math
\frac{\partial \mathbf W}{\partial t} + \nabla_\mathbf x \cdot \mathbf F = \mathbf S
```

The solution algorithms are organized using the finite volume method (FVM), and can be run in 1-3 dimensions on CPUs and GPUs.
The user interface is designed as intuitive and neat as possible, allowing users to focus on physics and to cooperate with existing packages in the Julialang ecosystem.




We strive for a user interface that makes _Kinetic.jl_ as friendly and intuitive to use as possible, allowing users to focus on the science. Internally, we have attempted to write the underlying algorithm so that the code runs as fast as possible for the configuration chosen by the user –- from simple two-dimensional setups to complex three-dimensional simulations –- and so that as much code as possible is shared between the CPU and GPU algorithms.

## Getting help

If you are interested in using _Kinetic.jl_ or are trying to figure out how to use it, please feel free to get in touch and raise questions.
Open an issue or pull request if you have questions, suggestions or solutions.
