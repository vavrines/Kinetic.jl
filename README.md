# Kinetic.jl

![](https://travis-ci.com/vavrines/Kinetic.jl.svg?branch=master)
![CI](https://github.com/vavrines/Kinetic.jl/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/vavrines/Kinetic.jl/branch/master/graph/badge.svg?token=mMtuTG3qMo)](https://codecov.io/gh/vavrines/Kinetic.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://xiaotianbai.com/Kinetic.jl/dev/)
![](https://zenodo.org/badge/243490351.svg)
<!--
[![Coverage Status](https://coveralls.io/repos/github/vavrines/Kinetic.jl/badge.svg?branch=master)](https://coveralls.io/github/vavrines/Kinetic.jl?branch=master)
-->

<img src="https://i.postimg.cc/ncXfgjXd/dancing-circles.gif" width="300"/>

This Julia package serves as a toolbox for theoretical and numerical studies towards the kinetic theory of gases, photons, plasmas and neutrons.
It employs the finite volume method (FVM) to conduct 1-3 dimensional numerical simulations on CPUs and GPUs, which solve the Boltzmann and its related model equations

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;f}{\partial&space;t}&plus;&space;\mathbf&space;u&space;\cdot&space;\nabla_\mathbf&space;x&space;f&space;&plus;&space;\mathbf&space;a&space;\cdot&space;\nabla_\mathbf&space;u&space;f&space;=&space;\int_{\mathcal&space;R^3}&space;\int_{\mathcal&space;S^2}&space;B&space;(f(\mathbf&space;u_*')f(\mathbf&space;u')-f(\mathbf&space;u_*)&space;f(\mathbf&space;u))&space;d\Omega&space;d\mathbf&space;u_*" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;f}{\partial&space;t}&plus;&space;\mathbf&space;u&space;\cdot&space;\nabla_\mathbf&space;x&space;f&space;&plus;&space;\mathbf&space;a&space;\cdot&space;\nabla_\mathbf&space;u&space;f&space;=&space;\int_{\mathcal&space;R^3}&space;\int_{\mathcal&space;S^2}&space;B&space;(f(\mathbf&space;u_*')f(\mathbf&space;u')-f(\mathbf&space;u_*)&space;f(\mathbf&space;u))&space;d\Omega&space;d\mathbf&space;u_*" title="\frac{\partial f}{\partial t}+ \mathbf u \cdot \nabla_\mathbf x f + \mathbf a \cdot \nabla_\mathbf u f = \int_{\mathcal R^3} \int_{\mathcal S^2} B (f(\mathbf u_*')f(\mathbf u')-f(\mathbf u_*) f(\mathbf u)) d\Omega d\mathbf u_*" /></a>

or their upscaling moment system

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;\mathbf&space;W}{\partial&space;t}&space;&plus;&space;\nabla_\mathbf&space;x&space;\cdot&space;\mathbf&space;F&space;=&space;\mathbf&space;S" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;\mathbf&space;W}{\partial&space;t}&space;&plus;&space;\nabla_\mathbf&space;x&space;\cdot&space;\mathbf&space;F&space;=&space;\mathbf&space;S" title="\frac{\partial \mathbf W}{\partial t} + \nabla_\mathbf x \cdot \mathbf F = \mathbf S" /></a>

## Documentation

For the detailed information on using the package,
[see the documentation](https://xiaotianbai.com/Kinetic.jl/dev/).

## Contributing

If you have further questions regarding Kinetic.jl or have got an idea on improving it, please feel free to get in touch. Open an issue or pull request if you'd like to work on a new feature or if you're new to open-source and want to find a cool little project or issue to work on that fits your interests. We're more than happy to help along the way.
