# Kinetic.jl

![](https://travis-ci.com/vavrines/Kinetic.jl.svg?branch=master)
<!--
[![Coverage Status](https://coveralls.io/repos/github/vavrines/Kinetic.jl/badge.svg?branch=master)](https://coveralls.io/github/vavrines/Kinetic.jl?branch=master)
-->

This Julia package serves as a tool box for theoretical and numerical studies towards the kinetic theory of gases, photons, plasmas and neutrons. 
It can be used to solve either Boltzmann and related model equations

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;f}{\partial&space;t}&plus;&space;\mathbf&space;u&space;\cdot&space;\nabla_\mathbf&space;x&space;f&space;&plus;&space;\mathbf&space;a&space;\cdot&space;\nabla_\mathbf&space;u&space;f&space;=&space;\int_{\mathcal&space;R^3}&space;\int_{\mathcal&space;S^2}&space;B&space;(f(\mathbf&space;u_*')f(\mathbf&space;u')-f(\mathbf&space;u_*)&space;f(\mathbf&space;u))&space;d\Omega&space;d\mathbf&space;u_*" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;f}{\partial&space;t}&plus;&space;\mathbf&space;u&space;\cdot&space;\nabla_\mathbf&space;x&space;f&space;&plus;&space;\mathbf&space;a&space;\cdot&space;\nabla_\mathbf&space;u&space;f&space;=&space;\int_{\mathcal&space;R^3}&space;\int_{\mathcal&space;S^2}&space;B&space;(f(\mathbf&space;u_*')f(\mathbf&space;u')-f(\mathbf&space;u_*)&space;f(\mathbf&space;u))&space;d\Omega&space;d\mathbf&space;u_*" title="\frac{\partial f}{\partial t}+ \mathbf u \cdot \nabla_\mathbf x f + \mathbf a \cdot \nabla_\mathbf u f = \int_{\mathcal R^3} \int_{\mathcal S^2} B (f(\mathbf u_*')f(\mathbf u')-f(\mathbf u_*) f(\mathbf u)) d\Omega d\mathbf u_*" /></a>

or their upscaling moment system

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;\mathbf&space;W}{\partial&space;t}&space;&plus;&space;\nabla_\mathbf&space;x&space;\cdot&space;\mathbf&space;F&space;=&space;\mathbf&space;S" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;\mathbf&space;W}{\partial&space;t}&space;&plus;&space;\nabla_\mathbf&space;x&space;\cdot&space;\mathbf&space;F&space;=&space;\mathbf&space;S" title="\frac{\partial \mathbf W}{\partial t} + \nabla_\mathbf x \cdot \mathbf F = \mathbf S" /></a>

The finite volume method (FVM) is employed in all cases. 

The package is compatible with Julia 1.3 or newer version. 
To make use of it, execute `Julia` and type
```julia
julia> ]
(v1.3) pkg> add Kinetic
```
This will install Kinetic and all its dependencies.
After that, load the package,
```julia
julia> using Kinetic
```
