---
title: 'Kinetic.jl: A lightweight finite volume toolbox in Julia'
tags:
  - computational fluid dynamics
  - kinetic theory
  - scientific machine learning
  - julia
authors:
  - name: Tianbai Xiao
    orcid: 0000-0001-9127-9497
    affiliation: 1
affiliations:
 - name: Karlsruhe Institute of Technology
   index: 1
date: 06 January 2021
bibliography: paper.bib
---

# Summary

``Kinetic.jl`` is a lightweight Julia toolbox for computational fluid dynamics and scientific machine learning. It focus on the theoretical and numerical studies of many-particle systems of gases, photons, plasmas, neutrons, etc. The finite volume method (FVM) is implemented to conduct 1-3 dimensional numerical simulations for any advection-diffusion type equations.

The main module consists of ``KitBase.jl`` with basic physics and ``KitML.jl`` with neural dynamics. The high-performance Fortran library ``KitFort.jl`` is not included by default, but can be manually imported when the executing efficiency becomes the priority. A Python wrapper ``kineticpy`` has been built as well to locate the structs and methods through ``pyjulia``. The ecosystem is designed to balance the efficiency for scientific research and the simplicity for educational usage.





An illustrative example of the shock tube problem

configuration file `config.txt`


```
# case
case = sod
space = 1d2f1v
nSpecies = 1
flux = kfvs
collision = bgk
interpOrder = 2
limiter = vanleer
boundary = fix
cfl = 0.5
maxTime = 0.2

# physical space
x0 = 0
x1 = 1
nx = 200
pMeshType = uniform
nxg = 1

# velocity space
vMeshType = rectangle
umin = -5
umax = 5
nu = 28
nug = 0

# gas
knudsen = 0.0001
mach = 0.0
prandtl = 1
inK = 2
omega = 0.81
alphaRef = 1.0
omegaRef = 0.5
```

and execute the following codes

```julia
using Kinetic
ks, ctr, face, t = initialize("config.txt")
t = solve!(ks, ctr, face, t)
plot_line(ks, ctr)
```

`solve!` is equivalent as the low-level solution algorithm

```julia
dt = timestep(ks, ctr, t)
nt = Int(floor(ks.set.maxTime / dt))
res = zeros(3)
for iter = 1:nt
    Kinetic.reconstruct!(ks, ctr)
    Kinetic.evolve!(ks, ctr, face, dt)
    Kinetic.update!(ks, ctr, face, dt, res)
end
```





``Kinetic.jl`` is interested in theoretical and numerical studies of many-particle systems of gases, photons, plasmas, neutrons, etc. It employs the finite volume method (FVM) to conduct 1-3 dimensional numerical simulations on CPUs and GPUs. Any advection-diffusion-type equation can be solved within the framework. Special attentions have been paid on Hilbert's sixth problem, i.e. to build the numerical passage between kinetic theory of gases, e.g. the Boltzmann equation



e.g. compressible gas dynamics.



Hilbert's sixth problem


cubersome
phase field



``Oceananigans.jl`` is designed for high-resolution simulations in idealized
geometries and supports direct numerical simulation, large eddy simulation,
arbitrary numbers of active and passive tracers, and linear and nonlinear
equations of state for seawater. Under the hood, ``Oceananigans.jl`` employs a
finite volume algorithm similar to that used by the Massachusetts Institute of
Technology general circulation model [@Marshall1997].

![Fig. 1](free_convection_and_baroclinic_instability.png)
Fig. 1: (Left) Large eddy simulation of small-scale oceanic boundary layer
turbulence forced by a surface cooling in a horizontally periodic domain using
$256^3$ grid points. The upper layer is well-mixed by turbulent convection and
bounded below by a strong buoyancy interface. (Right) Simulation of
instability of a horizontal density gradient in a rotating channel using
$256\times512\times128$ grid points. A similar process called baroclinic
instability acting on basin-scale temperature gradients fills the oceans with
eddies that stir carbon and heat. Plots made with `matplotlib` [@Hunter2007]
and `cmocean` [@Thyng2016].

``Oceananigans.jl`` leverages the Julia programming language [@Bezanson2017] to
implement high-level, low-cost abstractions, a friendly user interface, and a
high-performance model in one language and a common code base for execution on
the CPU or GPU with Julia’s native GPU compiler [@Besard2019]. Because Julia is
a high-level language, development is streamlined and users can flexibly specify
model configurations, set up arbitrary diagnostics and output, extend the code
base, and implement new features. Configuring a model with `architecture=CPU()`
or `architecture=GPU()` will execute the model on the CPU or GPU. By pinning a
simulation script against a specific version of Oceananigans, simulation results
are reproducible up to hardware differences.

Performance benchmarks show significant speedups when running on a GPU. Large
simulations on an Nvidia Tesla V100 GPU require ~1 nanosecond per grid point per
iteration. GPU simulations are therefore roughly 3x more cost-effective
than CPU simulations on cloud computing platforms such as Google Cloud. A GPU
with 32 GB of memory can time-step models with ~150 million grid points assuming
five fields are being evolved; for example, three velocity components and
tracers for temperature and salinity. These performance gains permit the
long-time integration of realistic simulations, such as large eddy simulation of
oceanic boundary layer turbulence over a seasonal cycle or the generation of
training data for turbulence parameterizations in Earth system models.

``Oceananigans.jl`` is continuously tested on CPUs and GPUs with unit tests,
integration tests, analytic solutions to the incompressible Navier-Stokes
equations, convergence tests, and verification experiments against published
scientific results. Future development plans include support for distributed
parallelism with CUDA-aware MPI as well as topography.

Ocean models that are similar to ``Oceananigans.jl`` include MITgcm
[@Marshall1997] and MOM6 [@Adcroft2019], both written in Fortran. However,
``Oceananigans.jl`` features a more efficient non-hydrostatic pressure solver
than MITgcm (and MOM6 is strictly hydrostatic). PALM [@Maronga2020] is Fortran
software for large eddy simulation of atmospheric and oceanic boundary layers
with complex boundaries on parallel CPU and GPU architectures. ``Oceananigans.jl``
is distinguished by its use of Julia which allows for a script-based interface as
opposed to a configuration-file-based interface used by MITgcm, MOM6, and PALM.
Dedalus [@Burns2020] is Python software with an intuitive script-based interface
that solves general partial differential equations, including the incompressible
Navier-Stokes equations, with spectral methods.

``Kinetic.jl`` is an open-source project hosted on GitHub and distributed under MIT license.

# Acknowledgements

The current work is funded by the Alexander von Humboldt Foundation

# References