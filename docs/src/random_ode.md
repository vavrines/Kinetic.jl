```@setup mysetup
using PolyChaos, DifferentialEquations
x0 = 2.0
μ, σ = -0.5, 0.05
tend, Δt = 3.0, 0.01
using PolyChaos
L, Nrec = 6, 40
opq = GaussOrthoPoly(L; Nrec=Nrec, addQuadrature=true)
using DifferentialEquations

a = [ convert2affinePCE(μ, σ, opq); zeros(Float64,L-1) ] # PCE coefficients of a
xinit = [ x0; zeros(Float64,L) ] # PCE coefficients of initial condition

t2 = Tensor(2, opq); # \langle \phi_i, \phi_j \rangle
t3 = Tensor(3, opq); # \langle \phi_i \phi_j, \phi_k \rangle

# Galerkin-projected random differential equation
function ODEgalerkin(du,u,p,t)
   du[:] = [ sum( p[j+1]*u[k+1]*t3.get([j,k,m])/t2.get([m,m]) for j=0:L for k=0:L) for m=0:L ]
end

probgalerkin = ODEProblem(ODEgalerkin,xinit,(0,tend),a)
solgalerkin = solve(probgalerkin;saveat=0:Δt:tend)
t, x = solgalerkin.t, solgalerkin.u;
# an advantage of PCE is that moments can be computed from the PCE coefficients alone; no sampling required
mean_pce = [ mean(x_, opq) for x_ in x]  
std_pce = [ std(x_, opq) for x_ in x]
using Statistics
Nsmpl = 5000
ξ = sampleMeasure(Nsmpl,opq)     # sample from Gaussian measure; effectively randn() here    
asmpl = evaluatePCE(a,ξ,opq)     # sample random variable with PCE coefficients a; effectively μ + σ*randn() here
# or: asmpl = samplePCE(Nsmpl,a,opq)
xmc = [ solve(ODEProblem((u,p,t)->aa*u,x0,(0,tend));saveat=0:Δt:tend).u for aa in asmpl]
xmc = hcat(xmc...);
[ mean(xmc,dims=2)-mean_pce std(xmc,dims=2)-std_pce]
logx_pce = [ log.(evaluatePCE(x[i],ξ,opq)) for i=1:length(t)]
[mean.(logx_pce)-(log(x0) .+ μ*t) std.(logx_pce)-σ*t ]
```
# Galerkin-based Solution of Random Differential Equation

This tutorial demonstrates how random differential equations can be solved using polynomial chaos expansions (PCE).

## Theory

A random differential equation is an ordinary differential equation that has random parameters, hence its solution is itself a (time-varying) random variable.
Perhaps the simplest non-trivial example is the following scalar, linear ordinary differential equation
```math
\dot{x}(t) = a x(t), \quad x(0) = x_{0},
```
where $a$ is the realization of a Gaussian random variable $\mathsf{a} \sim \mathcal{N}(\mu, \sigma^2)$ with mean $\mu$ and variance $\sigma^2$.
Arguably, for every realization $a$ we can solve the differential equation and obtain
```math
x(t) = x_0 \mathrm{e}^{a t},
```
from which we find that
```math
\ln (x(t)) = \ln (x_0) + at \sim \mathcal{N}(\ln(x_0) + \mu t, (\sigma t)^2).
```
In other words, the logarithm of the solution is normally distributed (so-called [log-normal distribution](https://en.wikipedia.org/wiki/Log-normal_distribution)).

We'd like to obtain this result numerically with the help of PCE.
The first step is to define the (truncated) PCE for the random variable $\mathsf{a}$
```math
\mathsf{a} = \sum_{i=0}^{L} a_i \phi_i,
```
where $a_i$ are the so-called PCE coefficients, and $\phi_i$ are the orthogonal basis polynomials.
As the solution to the random differential equation is itself a random variable, we treat $x(t)$ as the realization of the random variable $\mathsf{x}(t)$, and define its PCE
```math
\mathsf{x}(t) = \sum_{i=0}^{L} x_i(t) \phi_i.
```
The question is how to obtain the unknown PCE coefficients $x_i(t)$ from the known PCE coefficients $a_i$ relative to the orthogonal basis polynomials $\phi_i$.
This can be done using Galerkin projection, which is nothing else than projecting onto the orthogonal basis.
Think of a three-dimensional space, in which you have placed some three-dimensional object.
If you know project the silhouett of the object onto every axis of the three-dimensional space, then you are doing a Galerkin projection.
With PCE the concept is equivalent, but the imagination has a harder time.
The first step for Galerkin projection is to insert the PCEs
```math
\sum_{i=0}^{L} \dot{x}_i(t) \phi_i = \sum_{j=0}^{L} a_j \phi_j \sum_{k=0}^{L} x_k(t) \phi_k;
```
the second step is to project onto every basis polynomial $\phi_m$ for $m = 0, 1, \dots, L$, and to exploit orthogonality of the basis.
This gives
```math
\dot{x}_m(t) \langle \phi_m, \phi_m \rangle = \sum_{j=0}^{L} \sum_{k=0}^{L} a_j x_k(t) \langle \phi_l \phi_k, \phi_m \rangle \quad m = 0, 1, \dots, L.
```
Of course, the initial condition must not be forgotten:
```math
x_0(0) = x_0, \quad x_m(0) = 0 \quad m = 1, \dots, L.
```
If we can solve this enlarged system of ordinary random differential equations, we can reconstruct the analytic solution.

## Practice
We begin by defining the random differential equation


```@example mysetup
x0 = 2.0
μ, σ = -0.5, 0.05
tend, Δt = 3.0, 0.01
```

Next, we define an orthogonal basis (and its quadrature rule) relative to the Gaussian measure using `PolyChaos`.
We choose a maximum degree of `L`.


```@example mysetup
using PolyChaos
L, Nrec = 6, 40
opq = GaussOrthoPoly(L; Nrec=Nrec, addQuadrature=true)
```

Now we can define the PCE for $\mathsf{a}$ and solve the Galerkin-projected ordinary differential equation using `DifferentialEquations.jl`.


```@example mysetup
using DifferentialEquations

a = [ convert2affinePCE(μ, σ, opq); zeros(Float64,L-1) ] # PCE coefficients of a
xinit = [ x0; zeros(Float64,L) ] # PCE coefficients of initial condition

t2 = Tensor(2, opq); # \langle \phi_i, \phi_j \rangle
t3 = Tensor(3, opq); # \langle \phi_i \phi_j, \phi_k \rangle

# Galerkin-projected random differential equation
function ODEgalerkin(du,u,p,t)
   du[:] = [ sum( p[j+1]*u[k+1]*t3.get([j,k,m])/t2.get([m,m]) for j=0:L for k=0:L) for m=0:L ]
end

probgalerkin = ODEProblem(ODEgalerkin,xinit,(0,tend),a)
solgalerkin = solve(probgalerkin;saveat=0:Δt:tend)
t, x = solgalerkin.t, solgalerkin.u;
```

For later purposes we compute the expected value and the standard deviation at all time instants using PCE.


```@example mysetup
# an advantage of PCE is that moments can be computed from the PCE coefficients alone; no sampling required
mean_pce = [ mean(x_, opq) for x_ in x]  
std_pce = [ std(x_, opq) for x_ in x]
```

We compare the solution from PCE to a Monte-Carlo-based solution.
That means to solve the ordinary differential equation for many samples of $\mathsf{a}$.
We first sample from the measure using `sampleMeasure`, and then generate samples of $\mathsf{a}$ using `evaluatePCE`.
After that we solve the ODE and store the results in `xmc`.


```@example mysetup
using Statistics
Nsmpl = 5000
ξ = sampleMeasure(Nsmpl, opq)     # sample from Gaussian measure; effectively randn() here    
asmpl = evaluatePCE(a, ξ, opq)     # sample random variable with PCE coefficients a; effectively μ + σ*randn() here
# or: asmpl = samplePCE(Nsmpl,a,opq)
xmc = [ solve(ODEProblem((u,p,t)->aa*u,x0,(0,tend));saveat=0:Δt:tend).u for aa in asmpl]
xmc = hcat(xmc...);
```

Now we can compare the Monte Carlo mean and standard deviation to the expression from PCE for every time instant.


```@example mysetup
[ mean(xmc,dims=2)-mean_pce std(xmc,dims=2)-std_pce]
```

Clearly, the accuracy of PCE deteriorates over time.
Possible remedies are to increase the dimension of PCE, and to tweak the tolerances of the integrator.

Finally, we compare whether the samples follow a log-normal distribution, and compare the result to the analytic mean and standard deviation.


```@example mysetup
logx_pce = [ log.(evaluatePCE(x_,ξ,opq)) for x_ in x]
[ mean.(logx_pce)-(log(x0) .+ μ*t) std.(logx_pce)-σ*t ]
```

