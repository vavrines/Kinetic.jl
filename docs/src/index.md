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

# Getting help

If you are interested in using Kinetic.jl or are trying to figure out how to use it, please feel free to get in touch and raise questions.
Open an issue or pull request if you have questions, suggestions or solutions.
