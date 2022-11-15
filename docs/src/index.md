# Kinetic.jl

Kinetic is a computational fluid dynamics toolbox written in Julia. Based on differentiable programming, mechanical and neural network models are fused and solved in a unified framework. Simultaneous 1-3 dimensional numerical simulations can be performed on CPUs and GPUs.

The ecosystem follows the modular design philosophy. Depending on the specific use case, the main module is split into portable components to reduce the lantency caused by the LLVM [just-in-time](https://llvm.org/docs/tutorial/index.html#building-a-jit-in-llvm) compiler:
- [KitBase.jl](https://github.com/vavrines/KitBase.jl): physical models and numerical schemes
- [KitML.jl](https://github.com/vavrines/KitML.jl): neural models and machine learning methods
- [KitFort.jl](https://github.com/vavrines/KitFort.jl): optional high-performance Fortran backend
- [kineticpy](https://github.com/vavrines/kineticpy): Python interface built on top of [pyjulia](https://github.com/JuliaPy/pyjulia)

## Scope of application

Kinetic models and simulates fluid dynamics problems from the perspective of particle transport.
Any advection-diffusion-type equation of different particles, including molecules, photons, plasmas, neutrons, etc., can be solved within the framework. 
Special attentions have been paid on Hilbert's sixth problem, i.e. to build the numerical passage between [kinetic theory of gases](https://en.wikipedia.org/wiki/Kinetic_theory_of_gases), e.g. the Boltzmann equation, and continuum mechanics, e.g. the Euler and Navier-Stokes equations. A partial list of current supported models and equations include:
- linear Boltzmann equation
- nonlinear Boltzmann equation
- multi-component Boltzmann equation
- Fokker-Planck-Landau equation
- direct simulation Monte Carlo
- advection-diffusion equation
- Burgers' equation
- Euler equations
- Navier-Stokes equations
- Extended hydrodynamical equations from asymptotic expansion
- Magnetohydrodynamical equations
- Maxwell's equations

## Design philosophy

The code hierarchy is designed as intuitive and neat as possible.
It's dedicated to providing a friendly interface for educational usage in kinetic theory and rich functionality for scientific research.
Benefiting from the brilliant expressiveness and low-overhead abstraction provided by the [Julia programming language](https://julialang.org/), 
we provide different levels of APIs to allow the users to focus on physics and to cooperate with the existing packages in the Julia ecosystem.

## What is new?

Finite volume method is a proven approach for simulating conservation laws.
Compared with the existing open-source softwares, e.g. [OpenFOAM](https://openfoam.org/), [SU2](https://su2code.github.io/) and [Clawpack](https://www.clawpack.org/), 
Kinetic holds the novelty through the following points:
- 100% Julia stack that encounters no two-language problem
- Comprehensive support for kinetic theory and phase-space equations
- Lightweight design to ensure the flexibility for secondary development
- Closely coupling with scientific machine learning

## How to get help?

If you are interested in using Kinetic or are trying to figure out how to use it, please feel free to get in touch and raise questions.
Do open an issue or pull request if you have questions, suggestions, or solutions.