---
title: 'Kinetic.jl: A portable finite volume toolbox for scientific and neural computing'
tags:
  - kinetic theory
  - computational fluid dynamics
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

Kinetic.jl is a lightweight finite volume toolbox written in the Julia programming language for the study of computational physics and scientific machine learning.
It is an open-source project hosted on GitHub and distributed under MIT license.
The main module consists of KitBase.jl for basic physics and KitML.jl for neural dynamics.
The library provides a rich set of numerical flux functions and source terms.
Any advection-diffusion type mechanical or neural equation can be set up and solved in the framework.
The techniques of scientific machine learning can be integrated seamlessly to build data-driven closure models and accelerate calculation of nonlinear terms.
The package is designed to balance the programming flexibility for scientific research, the algorithmic efficiency for applications, and the simplicity for educational usage.

# Statement of need

A physical system can perform a wonderfully diverse set of acts on different characteristic scales.
For example, particle transports can be described by kinetic theory of gases at particle mean free path scale [@chapman1990] and fluid mechanics at macroscopic level [@batchelor2000].
With the rapidly advancing computing power, the finite volume method (FVM) is a prevalent method for quantitatively describing physical evolutions.

Most existing FVM libraries, e.g. the OpenFOAM [@jasak2007], are written in compiled languages (C/C++ and Fortran), which enjoy the perfect execution efficiency but sacrifice the flexibility for development.
For example, it is cumbersome to implement the phase-field evolution from the Boltzmann equation [@xiao2017;@xiao2020a] in OpenFOAM or integrate it with scientific machine learning (SciML) packages.
One compromise can be made with a combination of static and dynamic languages [@clawpack2020], where the high-level front-ends and the low-level computational back-ends are split.
Basically it benefits general users, while researchers still need to work on the back-end if a new feature is required. 
Besides, the two-language problem introduces additional trade-off in both development and execution.
Different from these packages, Kinetic.jl is built upon the Julia programming language [@bezanson2017], which is dynamically typed and designed for high performance computing for a broad range of devices. 
Based on type inference and multiple dispatch, it is a promising choice to solve the two-language problem.

Kinetic.jl focuses on the theoretical and numerical studies of many-particle systems of gases, photons, plasmas, neutrons, etc.
A hierarchy of abstractions is implemented.
At the highest level, it is feasible to model and simulate a fluid dynamic problem within 10 lines of code. 
At the lowest level, we design the methods for general numbers and arrays, so that it is possible to cooperate with existing packages in Julia ecosystem.
As an example, It uses Flux.jl [@Flux2018] to create and train scientific machine learning models.
The package holds the following innovations:

- 100% Julia stack that encounters no two-language problem

- Comprehensive support for kinetic theory and phase-space equations

- Lightweight design to ensure the flexibility for secondary development

- Close coupling with scientific machine learning

# KitBase.jl

The main module of Kinetic.jl is splitted into two pieces to reduce the just-in-time (JIT) compilation time for domain specific applications.
The basic physical laws and finite volume method are implemented in KitBase.jl.
It provides a variety of solvers for the Boltzmann equation, Maxwell's equations, advection-diffusion equation, Burgers' equation, Euler and Navier-Stokes equations, etc.
Different parallel computing techniques are provided, e.g. multi-threading, distributed computing and CUDA programming.

In the following, we present an illustrative example of solving lid-driven cavity problem with the Boltzmann equation. 
Two initialization methods, i.e. configuration text and Julia script, are available for setting up the solver.
With the configuration file `config.toml` set as below,
```toml
# setup
matter = gas # material
case = cavity # case
space = 2d2f2v # phase
flux = kfvs # flux function
collision = bgk # intermolecular collision
nSpecies = 1 # number of species
interpOrder = 2 # interpolation order of accuracy
limiter = vanleer # limiter function
boundary = maxwell # boundary condition
cfl = 0.8 # CFL number
maxTime = 5.0 # maximal simulation time

# physical space
x0 = 0.0 # starting point in x
x1 = 1.0 # ending point in x
nx = 45 # number of cells in x
y0 = 0.0 # starting point in y
y1 = 1.0 # ending point in y
ny = 45 # number of cells in y
pMeshType = uniform # mesh type
nxg = 0 # number of ghost cell in x
nyg = 0 # number of ghost cell in y

# velocity space
umin = -5.0 # starting point in u
umax = 5.0 # ending point in u
nu = 28 # number of cells in u
vmin = -5.0 # starting point in v
vmax = 5.0 # ending point in v
nv = 28 # number of cells in v
vMeshType = rectangle # mesh type
nug = 0 # number of ghost cell in u
nvg = 0 # number of ghost cell in v

# gas property
knudsen = 0.075 # Knudsen number
mach = 0.0 # Mach number
prandtl = 1.0 # Prandtl number
inK = 1.0 # molecular inner degree of freedom
omega = 0.72 # viscosity index of hard-sphere gas
alphaRef = 1.0 # viscosity index of hard-sphere gas in reference state
omegaRef = 0.5 # viscosity index of hard-sphere gas ub reference state

# boundary condition
uLid = 0.15 # U-velocity of moving wall
vLid = 0.0 # V-velocity of moving wall
tLid = 1.0 # temperature of wall
```

let us execute the following codes
```julia
using Kinetic
set, ctr, xface, yface, t = initialize("config.toml")
t = solve!(set, ctr, xface, yface, t)
plot_contour(set, ctr)
```

The computational setup is stored in `set`. 
The solutions over control volumes are represented in an array `ctr`, while `xface` and `yface` record the interface fluxes along x and y directions.
The result is visualized with built-in function `plot_contour`, which presents the distributions of gas density, velocity and temperature inside the cavity.

![Fig. 1](cavity.png)
Fig. 1: macroscopic variables in the lid-driven cavity (topleft: density, top right: U-velocity, bottomleft: V-velocity, bottomright: temperature).

# KitML.jl

Machine learning is building its momentum in scientific computing.
Given the nonlinear structure of differential and integral equations, it is promising to incorporate the universal function approximator from machine learning models into the governing equations and achieve the balance between efficiency and accuracy.
In KitML.jl, we provide strategies to construct hybrid mechanical-neural models.
The detailed theory and implementation can be found in the paper [@xiao2020b].

# Extension

Numerical simulations of nonlinear models and differential equations are essentially connected with supercomputers and high-performance computing (HPC). 
Considering that some existing hardware architecture, e.g. Sunway TaihuLight with Chinese-designed SW26010 processors, only provides optimization for specific languages, we have developed an accompanying package KitFort.jl.
It is not a default component of Kinetic.jl, but can be manually imported.
Besides, a wrapper kineticpy has been built as well to locate the structures and methods from the Python ecosystem. 

# Acknowledgements

The current work is funded by the Alexander von Humboldt Foundation (Ref3.5-CHN-1210132-HFST-P).

# References
