# KitFort and high performance computing

Numerical simulations of nonlinear models and differential equations are essentially connected with supercomputers and high-performance computing (HPC).
The performance of a supercomputer or a software program is commonly measured in floating-point operations per second (FLOPS).
Through the milestone astronomy research of [Celeste](https://juliacomputing.com/case-studies/celeste/), Julia has entered the PetaFLOPS club (together with C/C++ and Fortran) since 2017.
Julia is experiencing a dramatic Rise in HPC and elsewhere, and that is why we use Julia to organize the Kinetic.
However, compared with the mature C/C++ ecosystem, the equivalent execution efficiency isn't going to happen in all time and situations.
Some existing hardware architecture, e.g. [Sunway TaihuLight](https://en.wikipedia.org/wiki/Sunway_TaihuLight), the previou fastest supercomputer in [TOP500](https://www.top500.org/) list, is built upon 40,960 Chinese-designed SW26010 manycore 64-bit RISC processors, which is not specifically optimized for Julia.
Therefore, we've develop an accompanying package [KitFort.jl](https://github.com/vavrines/KitFort.jl).
The Fortran codes have been linked to the Julia syntax with the built-in `ccall` function.
It's not a default submodule of Kinetic since we believe the Julia codes are sufficient for general users and developers and encounter no two-language problem.
However, it can be manually imported when the executing efficiency becomes the first priority by executing
```julia
julia> ]
(v1.8) pkg> add KitFort
```
After that, using/import the package.
```julia
julia> using KitFort
```

It can be updated to the latest tagged release from the package manager by executing
```julia
(v1.8) pkg> update KitFort
```
