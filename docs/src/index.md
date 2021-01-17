# Kinetic.jl

Kinetic is a lightweight Julia toolbox for the study of computational fluid dynamics and scientific machine learning.
The main module consists of [KitBase.jl](https://github.com/vavrines/KitBase.jl) with basic physics and [KitML.jl](https://github.com/vavrines/KitML.jl) with neural dynamics. 
The high-performance Fortran library [KitFort.jl](https://github.com/vavrines/KitFort.jl) is optional and can be manually imported when the executing efficiency becomes the first priority. 
A wrapper [kineticpy](https://github.com/vavrines/kineticpy) has been built to locate the data hierarchies and methods in Python.

Kinetic is interested in theoretical modeling and numerical simulation of many-particle systems, e.g. gases, photons, plasmas, neutrons, electrons, etc.
It employs the finite volume method (FVM) to conduct 1-3 dimensional numerical simulations on CPUs and GPUs.
Any advection-diffusion type equation can be hooked and solved within the framework.
The user interface is designed as intuitive and neat as possible, allowing users to focus on physics and to cooperate with existing packages in the Julia ecosystem.

## Getting help

If you are interested in using Kinetic.jl or are trying to figure out how to use it, please feel free to get in touch and raise questions.
Open an issue or pull request if you have questions, suggestions or solutions.
