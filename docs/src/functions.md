# Functions

!!! note
    The core interface of all essential functions are *not* dependent on specialized types such as `AbstractOrthoPoly`.
    Having said that, for exactly those essential functions there exist overloaded functions that accept specialized types such as `AbstractOrthoPoly` as arguments.

    Too abstract?
    For example, the function `evaluate` that evaluates a polynomial of degree `n` at points `x` has the core interface
    ```
        evaluate(n::Int,x::Array{<:Real},a::Vector{<:Real},b::Vector{<:Real})
    ```
    where `a` and `b` are the vectors of recurrence coefficients.
    For simplicity, there also exists the interface
    ```
        evaluate(n::Int64,x::Vector{<:Real},op::AbstractOrthoPoly)
    ```
    So fret not upon the encounter of multiply-dispatched versions of the same thing. It's there to simplify your life.

    The idea of this approach is to make it simpler for others to copy and paste code snippets and use them in their own work.

List of all functions in `PolyChaos`.

```@index
```

## Recurrence Coefficients for Monic Orthogonal Polynomials
The functions below provide analytic expressions for the recurrence coefficients of common orthogonal polynomials.
All of these provide *monic orthogonal polynomials* relative to the weights.

!!! note
    The number `N` of recurrence coefficients has to be positive for all functions below.

```@docs
r_scale
rm_compute
rm_logistic
rm_hermite
rm_hermite_prob
rm_laguerre
rm_legendre
rm_legendre01
rm_jacobi
rm_jacobi01
rm_meixner_pollaczek
stieltjes
lanczos
mcdiscretization
```

## Show Orthogonal Polynomials

To get a human-readable output of the orthognoal polynomials there is the function `showpoly`

```@docs
showpoly
```

In case you want to see the entire basis, just use `showbasis`
```@docs
showbasis
```

Both of these functions make excessive use of
```@docs
rec2coeff
```

## Evaluate Orthogonal Polynomials
```@docs
evaluate
```

## Scalar Products
```@docs
computeSP2
computeSP
```

## Quadrature Rules
```@docs
fejer
fejer2
clenshaw_curtis
quadgp
gauss
radau
lobatto
```

## Polynomial Chaos
```@docs
mean
var
std
sampleMeasure
evaluatePCE
samplePCE
calculateAffinePCE
convert2affinePCE
```

## Auxiliary Functions
```@docs
nw
coeffs
integrate
PolyChaos.issymmetric
```
