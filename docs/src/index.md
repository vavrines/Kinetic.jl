# Overview

Kinetic.jl is a collection of numerical routines for orthogonal polynomials written in the [Julia](https://julialang.org/) programming language.
Starting from some non-negative weight (aka an absolutely continuous nonnegative measure), PolyChaos allows
- to compute the coefficients for the monic three-term recurrence relation,
- to evaluate the orthogonal polynomials at arbitrary points,
- to compute the quadrature rule,
- to compute tensors of scalar products,
- to do all of the above in a multivariate setting (aka product measures).

```@docs
pdf_slope
maxwellian
```