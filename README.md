<div align="center"> <img
  src="https://i.postimg.cc/ncXfgjXd/dancing-circles.gif"
  alt="Kinetic Logo" width="300"></img>
</div>
<!--
# Kinetic.jl
<img src="https://i.postimg.cc/ncXfgjXd/dancing-circles.gif" width="300"/>
-->

![](https://travis-ci.com/vavrines/Kinetic.jl.svg?branch=master)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://xiaotianbai.com/Kinetic.jl/dev/)
![](https://zenodo.org/badge/243490351.svg)

**Kinetic.jl** is a lightweight Julialang toolbox for kinetic theory and scientific machine learning. The ecosystem consists of the following elements. Currently,the main module here contains [KitBase.jl](https://github.com/vavrines/KitBase.jl) with basic physics and [KitML.jl](https://github.com/vavrines/KitML.jl) with neural differential equations. The high-performance Fortran library [KitFort.jl](https://github.com/vavrines/KitFort.jl) can be manually imported when needed.

| Kinetic | KitBase | KitML | KitFort |
| ----------   | --------- | ---------------- | ------ |
| ![CI](https://github.com/vavrines/Kinetic.jl/workflows/CI/badge.svg) | ![CI](https://github.com/vavrines/KitBase.jl/workflows/CI/badge.svg) | ![CI](https://github.com/vavrines/KitML.jl/workflows/CI/badge.svg) | ![CI](https://github.com/vavrines/KitFort.jl/workflows/CI/badge.svg) |
| [![codecov](https://codecov.io/gh/vavrines/Kinetic.jl/branch/master/graph/badge.svg?token=mMtuTG3qMo)](https://codecov.io/gh/vavrines/Kinetic.jl) | [![codecov](https://codecov.io/gh/vavrines/KitBase.jl/branch/main/graph/badge.svg?token=vGgQhyGJ6L)](https://codecov.io/gh/vavrines/KitBase.jl) | [![codecov](https://codecov.io/gh/vavrines/KitML.jl/branch/main/graph/badge.svg?token=OnazyqLA4K)](https://codecov.io/gh/vavrines/KitML.jl) | [![codecov](https://codecov.io/gh/vavrines/KitFort.jl/branch/main/graph/badge.svg?token=67tfVc3AtW)](https://codecov.io/gh/vavrines/KitFort.jl) |

## Theory

This Kinetic.jl serves as a ecosystem for theoretical and numerical studies of the kinetic theory of gases, photons, plasmas, neutrons, etc.
It employs the finite volume method (FVM) to conduct 1-3 dimensional numerical simulations on CPUs and GPUs, which solve the Boltzmann

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;f}{\partial&space;t}&plus;&space;\mathbf&space;u&space;\cdot&space;\nabla_\mathbf&space;x&space;f&space;&plus;&space;\mathbf&space;a&space;\cdot&space;\nabla_\mathbf&space;u&space;f&space;=&space;\int_{\mathcal&space;R^3}&space;\int_{\mathcal&space;S^2}&space;B&space;(f(\mathbf&space;u_*')f(\mathbf&space;u')-f(\mathbf&space;u_*)&space;f(\mathbf&space;u))&space;d\Omega&space;d\mathbf&space;u_*" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;f}{\partial&space;t}&plus;&space;\mathbf&space;u&space;\cdot&space;\nabla_\mathbf&space;x&space;f&space;&plus;&space;\mathbf&space;a&space;\cdot&space;\nabla_\mathbf&space;u&space;f&space;=&space;\int_{\mathcal&space;R^3}&space;\int_{\mathcal&space;S^2}&space;B&space;(f(\mathbf&space;u_*')f(\mathbf&space;u')-f(\mathbf&space;u_*)&space;f(\mathbf&space;u))&space;d\Omega&space;d\mathbf&space;u_*" title="\frac{\partial f}{\partial t}+ \mathbf u \cdot \nabla_\mathbf x f + \mathbf a \cdot \nabla_\mathbf u f = \int_{\mathcal R^3} \int_{\mathcal S^2} B (f(\mathbf u_*')f(\mathbf u')-f(\mathbf u_*) f(\mathbf u)) d\Omega d\mathbf u_*" /></a>

and its related deep neural models, or their upscaling moment system

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;\mathbf&space;W}{\partial&space;t}&space;&plus;&space;\nabla_\mathbf&space;x&space;\cdot&space;\mathbf&space;F&space;=&space;\mathbf&space;S" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;\mathbf&space;W}{\partial&space;t}&space;&plus;&space;\nabla_\mathbf&space;x&space;\cdot&space;\mathbf&space;F&space;=&space;\mathbf&space;S" title="\frac{\partial \mathbf W}{\partial t} + \nabla_\mathbf x \cdot \mathbf F = \mathbf S" /></a>

## Documentation

For the detailed information on using the package, please
[check the documentation](https://xiaotianbai.com/Kinetic.jl/dev/).

## Contributing

If you have further questions regarding Kinetic.jl or have got an idea on improving it, please feel free to get in touch. Open an issue or pull request if you'd like to work on a new feature or even if you're new to open-source and want to find a cool little project or issue to work on that fits your interests. We're more than happy to help along the way.
