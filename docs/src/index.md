# Kinetic.jl

Kinetic is a portable toolbox for the study of computational fluid dynamics and scientific machine learning.
The default module consists of [KitBase.jl](https://github.com/vavrines/KitBase.jl) with basic physics and [KitML.jl](https://github.com/vavrines/KitML.jl) with neural dynamics. 
The high-performance Fortran library [KitFort.jl](https://github.com/vavrines/KitFort.jl) can be manually imported when the executing efficiency becomes the first priority. 
A wrapper [kineticpy](https://github.com/vavrines/kineticpy) has been built to locate the data hierarchies and methods in Python.

## Scope of application

Kinetic is interested in the evolution of many-particle systems, e.g. gases, photons, plasmas, neutrons, electrons, etc.
Based on the finite volume method (FVM), it provides a direct modeling tool in the discrete space,
where 1-3 dimensional theoretical modeling and numerical simulation can be conducted.
Any advection-diffusion type equation can be hooked and solved within the framework.
Special attention has been paid to the [kinetic theory](https://en.wikipedia.org/wiki/Kinetic_theory_of_gases) and the Boltzmann-type equations,
which depicts the time-space evolution of particles via ensemble averaging at the mesoscopic level.
A partial list of current supported models and equations include:
- linear Boltzmann equation;
- nonlinear Boltzmann equation;
- multi-component Boltzmann equations;
- Fokker-Planck-Landau equation;
- direct simulation Monte Carlo;
- advection-diffusion equation;
- Burgers' equation;
- Euler equations;
- Navier-Stokes equations;
- Extended hydrodynamical equations from gas kinetic expansion;
- Magnetohydrodynamical equations;
- Maxwell's equations.

## Design philosophy

The code hierarchy is designed as intuitive and neat as possible.
It's dedicated to providing a friendly interface for educational usage in kinetic theory and rich functionality for scientific research.
Benefiting from the brilliant expressiveness and low-overhead abstraction provided by the [Julia programming language](https://julialang.org/), 
we provide a combination of high-level and low-level APIs to allow the users to focus on physics and to cooperate with the existing packages in the Julia ecosystem.

## What is new?

Finite volume method is a proven approach for simulating conservation laws.
Compared with the existing open-source softwares, e.g. [OpenFOAM](https://openfoam.org/), [SU2](https://su2code.github.io/) and [Clawpack](https://www.clawpack.org/), 
Kinetic holds the novelty through the following points:
- 100% Julia stack that encounters no two-language problem
- Comprehensive support for kinetic theory and phase-space equations
- Lightweight design to ensure the flexibility for secondary development
- Closely coupling with scientific machine learning

## How to get help?

If you are interested in using Kinetic.jl or are trying to figure out how to use it, please feel free to get in touch and raise questions.
Open an issue or pull request if you have questions, suggestions or solutions.
