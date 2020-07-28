```@setup mysetup
using PolyChaos
```
# Overview
PolyChaos is a collection of numerical routines for orthogonal polynomials written in the [Julia](https://julialang.org/) programming language.
Starting from some non-negative weight (aka an absolutely continuous nonnegative measure), PolyChaos allows
- to compute the coefficients for the monic three-term recurrence relation,
- to evaluate the orthogonal polynomials at arbitrary points,
- to compute the quadrature rule,
- to compute tensors of scalar products,
- to do all of the above in a multivariate setting (aka product measures).

If the weight function is a probability density function, PolyChaos further provides routines to compute [polynomial chaos expansions](https://en.wikipedia.org/wiki/Polynomial_chaos) (PCEs) of random variables with this very density function.
These routines allow
- to compute affine PCE coefficients for arbitrary densities,
- to compute moments,
- to compute the tensors of scalar products.

PolyChaos contains several *canonical* orthogonal polynomials such as Jacobi or Hermite polynomials.
For these, closed-form expressions and state-of-the art quadrature rules are used whenever possible.
However, a cornerstone principle of PolyChaos is to provide all the functionality for user-specific, arbitrary weights.

!!! note
    What PolyChaos is not (at least currently):
    - a self-contained introduction to orthogonal polynomials, quadrature rules and/or polynomial chaos expansions. We assume the user brings some experience to the table. However, over time we will focus on strengthening the tutorial charater of the package.
    - a symbolic toolbox
    - a replacement for [FastGaussQuadrature.jl](https://github.com/ajt60gaibb/FastGaussQuadrature.jl)

## Installation
The package requires `Julia 1.0` or newer.
In `Julia` switch to the package manager
```julia
julia> ]
(v1.0) pkg> add PolyChaos
```
This will install PolyChaos and its dependencies.
Once that is done, load the package:
```julia
julia> using PolyChaos
```
That's it.

Let's take a look at a simple example.
We would like to solve the integral
```math
\int_0^1 6 x^5 \mathrm{d}x.
```
Exploiting the underlying uniform measure, the integration can be done exactly with a 3-point quadrature rule.
```jldoctest
julia> using PolyChaos

julia> opq = Uniform01OrthoPoly(3, addQuadrature = true)
Uniform01OrthoPoly{Array{Float64,1},Uniform01Measure,Quad{Float64,Array{Float64,1}}}(3, [0.5, 0.5, 0.5, 0.5], [1.0, 0.08333333333333333, 0.06666666666666667, 0.06428571428571428], Uniform01Measure(PolyChaos.w_uniform01, (0.0, 1.0), true), Quad{Float64,Array{Float64,1}}("golubwelsch", 3, [0.11270166537925838, 0.49999999999999994, 0.8872983346207417], [0.2777777777777777, 0.4444444444444444, 0.27777777777777757]))

julia> integrate(x -> 6x^5, opq)
0.9999999999999993

julia> show(opq)

Univariate orthogonal polynomials
degree:         3
#coeffs:        4
α =             [0.5, 0.5, 0.5, 0.5]
β =             [1.0, 0.08333333333333333, 0.06666666666666667, 0.06428571428571428]

Measure dλ(t)=w(t)dt
w:      w_uniform01
dom:    (0.0, 1.0)
symmetric:      true
```

To get going with PolyChaos check out the tutorials such as the one on [numerical integration](@ref NumericalIntegration).
In case you are unfamiliar with orthogonal polynomials, perhaps [this background information](@ref MathematicalBackground) is of help.

## References
The code base of `PolyChaos` is partially based on Walter Gautschi's [Matlab suite of programs for generating orthogonal polynomials and related quadrature rules](https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html), with much of the theory presented in his book *Orthogonal Polynomials: Computation and Approximation* published in 2004 by the Oxford University Press.

For the theory of polynomial chaos expansion we mainly consulted T. J. Sullivan. *Introduction to Uncertainty Quantification*. Springer International Publishing Switzerland. 2015.

## Contributing
We are always looking for contributors.
If you are interested, just get in touch: tillmann [dot] muehlpfordt [at] kit [dot] edu.

Or just fork and/or star the repository:

Julia's package manager works nicely with Github: simply install the hosted package via `Pkg.clone` and the [repository's URL](https://help.github.com/articles/which-remote-url-should-i-use/).
A fork is created with

```julia
Pkg.clone("https://github.com/timueh/PolyChaos.jl")
```
The fork will replace the original package.

Call
```julia
Pkg.dir("PolyChaos")
```
to figure out where the package was cloned to.
Go to that location and figure out what branch you are on via `git branch`.

## Citing
If you found the software useful and applied it to your own research, we'd appreciate a citation.
Add the following to your BibTeX file

```tex
@misc{Muehlpfordt18,
  author       = {T. Mühlpfordt and F. Zahn and F. Becker and T. Faulwasser and V. Hagenmeyer},
  title        = {PolyChaos.jl -- A Julia package for orthogonal polynomials, quadrature, and polynomial chaos expansion},
  year         = 2019,
  doi          = {10.5281/zenodo.2592031},
  url          = {https://doi.org/10.5281/zenodo.2592031}
}
```

Of course you are more than welcome to partake in GitHub's gamification aka :star: and :fork_and_knife:.
