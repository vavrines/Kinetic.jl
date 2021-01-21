# Kinetic.jl

Kinetic is a portable Julia toolbox for the study of computational fluid dynamics and scientific machine learning.
The main module consists of [KitBase.jl](https://github.com/vavrines/KitBase.jl) with basic physics and [KitML.jl](https://github.com/vavrines/KitML.jl) with neural dynamics. 
The high-performance Fortran library [KitFort.jl](https://github.com/vavrines/KitFort.jl) is optional and can be manually imported when the executing efficiency becomes the first priority. 
A wrapper [kineticpy](https://github.com/vavrines/kineticpy) has been built to locate the data hierarchies and methods in Python.

## Topics

Kinetic is interested in theoretical modeling and numerical simulation of many-particle systems, e.g. gases, photons, plasmas, neutrons, electrons, etc.
It employs the finite volume method (FVM) to conduct 1-3 dimensional numerical simulations on CPUs and GPUs.
Any advection-diffusion type equation can be hooked and solved within the framework.
The user interface is designed as intuitive and neat as possible.
The combined development of high-level and low-level APIs allows users to focus on physics and to cooperate with existing packages in the Julia ecosystem.

## Novelty

Finite volume method is a proven approach for simulating conservation laws.
Compared with the existing open-source softwares, e.g. OpenFOAM, SU2 and Clawpack, Kinetic holds the following innovative points
- 100% Julia stack with no two-language problem
- Well-established support for the Boltzmann and other integro-differential equations in phase space
- Closely coupled with scientific machine learning

## Getting help

If you are interested in using Kinetic.jl or are trying to figure out how to use it, please feel free to get in touch and raise questions.
Open an issue or pull request if you have questions, suggestions or solutions.
