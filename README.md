<div align="center"> <img
  src="https://i.postimg.cc/ncXfgjXd/dancing-circles.gif"
  alt="Kinetic Logo" width="300"></img>
</div>
<!--
# Kinetic.jl
<img src="https://i.postimg.cc/ncXfgjXd/dancing-circles.gif" width="300"/>
-->

[![version](https://juliahub.com/docs/Kinetic/version.svg)](https://juliahub.com/ui/Packages/Kinetic/wrVmu)
![](https://travis-ci.com/vavrines/Kinetic.jl.svg?branch=master)
[![](https://img.shields.io/badge/docs-stable-green.svg)](https://xiaotianbai.com/Kinetic.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-green.svg)](https://xiaotianbai.com/Kinetic.jl/dev/)
![](https://zenodo.org/badge/243490351.svg)

**Kinetic.jl** is a lightweight [Julia](https://julialang.org) toolbox for the study of computational fluid dynamics and scientific machine learning. The main module consists of [KitBase.jl](https://github.com/vavrines/KitBase.jl) with basic physics and [KitML.jl](https://github.com/vavrines/KitML.jl) with neural dynamics. The high-performance Fortran library [KitFort.jl](https://github.com/vavrines/KitFort.jl) can be manually imported when the executing efficiency becomes the first priority. A Python wrapper [kineticpy](https://github.com/vavrines/kineticpy) has been built to locate the structs and methods here through [pyjulia](https://github.com/JuliaPy/pyjulia). The status of continuous integration and coverage for the ecosystem is listed below.

| Kinetic | KitBase | KitML | KitFort |
| ----------   | --------- | ---------------- | ------ |
| ![CI](https://github.com/vavrines/Kinetic.jl/workflows/CI/badge.svg) | ![CI](https://github.com/vavrines/KitBase.jl/workflows/CI/badge.svg) | ![CI](https://github.com/vavrines/KitML.jl/workflows/CI/badge.svg) | ![CI](https://github.com/vavrines/KitFort.jl/workflows/CI/badge.svg) |
| [![codecov](https://img.shields.io/codecov/c/github/vavrines/Kinetic.jl.svg)](https://codecov.io/gh/vavrines/Kinetic.jl) | [![codecov](https://img.shields.io/codecov/c/github/vavrines/KitBase.jl.svg)](https://codecov.io/gh/vavrines/KitBase.jl) | [![codecov](https://img.shields.io/codecov/c/github/vavrines/KitML.jl.svg)](https://codecov.io/gh/vavrines/KitML.jl) | [![codecov](https://img.shields.io/codecov/c/github/vavrines/KitFort.jl.svg)](https://codecov.io/gh/vavrines/KitFort.jl) |

## Motivation

Kinetic.jl is interested in theoretical and numerical studies of many-particle systems of gases, photons, plasmas, neutrons, etc.
It employs the finite volume method (FVM) to conduct 1-3 dimensional numerical simulations on CPUs and GPUs.
Any advection-diffusion-type equation can be solved within the framework.
Special attentions have been paid on Hilbert's sixth problem, i.e. to build the numerical passage between kinetic theory of gases, e.g. the Boltzmann equation

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;f}{\partial&space;t}&plus;\mathbf{v}&space;\cdot&space;\nabla_{\mathbf{x}}&space;f&space;=&space;\int_{\mathbb&space;R^3}&space;\int_{\mathcal&space;S^2}&space;\mathcal&space;B(\cos&space;\beta,&space;|\mathbf{v}-\mathbf{v_*}|)&space;\left[&space;f(\mathbf&space;v')f(\mathbf&space;v_*')-f(\mathbf&space;v)f(\mathbf&space;v_*)\right]&space;d\mathbf&space;\Omega&space;d\mathbf&space;v_*" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\frac{\partial&space;f}{\partial&space;t}&plus;\mathbf{v}&space;\cdot&space;\nabla_{\mathbf{x}}&space;f&space;=&space;\int_{\mathbb&space;R^3}&space;\int_{\mathcal&space;S^2}&space;\mathcal&space;B(\cos&space;\beta,&space;|\mathbf{v}-\mathbf{v_*}|)&space;\left[&space;f(\mathbf&space;v')f(\mathbf&space;v_*')-f(\mathbf&space;v)f(\mathbf&space;v_*)\right]&space;d\mathbf&space;\Omega&space;d\mathbf&space;v_*" title="\frac{\partial f}{\partial t}+\mathbf{v} \cdot \nabla_{\mathbf{x}} f = \int_{\mathbb R^3} \int_{\mathcal S^2} \mathcal B(\cos \beta, |\mathbf{v}-\mathbf{v_*}|) \left[ f(\mathbf v')f(\mathbf v_*')-f(\mathbf v)f(\mathbf v_*)\right] d\mathbf \Omega d\mathbf v_*" /></a>

and continuum mechanics, e.g. the Euler and Navier-Stokes equations

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;\mathbf&space;W}{\partial&space;t}&space;&plus;&space;\nabla_\mathbf&space;x&space;\cdot&space;\mathbf&space;F&space;=&space;\mathbf&space;S" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\frac{\partial&space;\mathbf&space;W}{\partial&space;t}&space;&plus;&space;\nabla_\mathbf&space;x&space;\cdot&space;\mathbf&space;F&space;=&space;\mathbf&space;S" title="\frac{\partial \mathbf W}{\partial t} + \nabla_\mathbf x \cdot \mathbf F = \mathbf S" /></a>

A list of current supported equations include:
- linear Boltzmann equation
- nonlinear Boltzmann equation
- nonlinear kinetic model equation
- multi-component Boltzmann equations
- advection-diffusion equation
- Burgers equation
- Euler equations
- Navier-Stokes equations
- Extended hydrodynamical equations from gas kinetic expansion
- Magnetohydrodynamical equations
- Poisson equation
- Fokker-Planck-Landau equation
- Maxwell's equations

## Documentation

For the detailed information on the implementation and usage of the package, please
[check the documentation](https://xiaotianbai.com/Kinetic.jl/dev/).

## Contributing

If you have further questions regarding Kinetic.jl or have got an idea on improving it, please feel free to get in touch. Open an issue or pull request if you'd like to work on a new feature or even if you're new to open-source and want to find a cool little project or issue to work on that fits your interests. We're more than happy to help along the way.
