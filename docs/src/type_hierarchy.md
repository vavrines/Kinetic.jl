!!! note
    If you are unfamiliar with the mathematical background of orthogonal polynomials, check out [this tutorial](@ref MathematicalBackground).
# Type Hierarchy

Let's look at the types `PolyChaos` provides.
There are four `AbstractTypes`: `AbstractMeasure`, `AbstractOrthoPoly`, `AbstractQuad`, and `AbstractTensor`.
`AbstractMeasure` is the core on which `AbstractOrthoPoly` builds, on which `AbstractQuad` builds, which is then used by `AbstractTensor`.

## AbstractMeasure

The type tree for `AbstractMeasure` looks as follows

```julia
julia> using AbstractTrees, PolyChaos
julia> AbstractTrees.children(x::Type) = subtypes(x)
julia> print_tree(AbstractMeasure)
AbstractMeasure
├─ AbstractCanonicalMeasure
│  ├─ Beta01Measure
│  ├─ GammaMeasure
│  ├─ GaussMeasure
│  ├─ HermiteMeasure
│  ├─ JacobiMeasure
│  ├─ LaguerreMeasure
│  ├─ LegendreMeasure
│  ├─ LogisticMeasure
│  ├─ MeixnerPollaczekMeasure
│  ├─ Uniform01Measure
│  ├─ Uniform_11Measure
│  ├─ genHermiteMeasure
│  └─ genLaguerreMeasure
├─ Measure
└─ ProductMeasure
```



There are several canonical measures that `PolyChaos` provides, all gathered in as subtypes of `AbstractCanonicalMeasure`.
The `Measure` type is a generic measure, and `ProductMeasure` has an obvious meaning.

What are the relevant fields?

### Measure
It all begins with a measure, more specifically absolutely continuous measures.
What are the fields of such a type `Measure`?

| Field | Meaning |
| --- | --- |
| `name::String` | Name of measure|
| `w::Function` | Weight function $w: \Omega \rightarrow \mathbb{R}$ |
| `dom::Tuple{Real,Real}` | Domain $ \Omega$|
| `symmetric::Bool` | Is $w$ symmetric relative to some $m \in \Omega$, hence $w(m-x) = w(m+x)$ for all $x \in \Omega$? |
| `pars::Dict` | Additional parameters (e.g. shape parameters for Beta distribution) |

They are a `name`, a weight function $w: \Omega \rightarrow \mathbb{R}$ with domain $\Omega$ (`dom`).
If the weight function is symmetric relative to some $m \in \Omega$, the field `symmetric` should be set to `true`.
Symmetry relative to $m$ means that
```math
\forall x \in \Omega: \quad w(m-x) = w(m+x).
```
For example, the Gaussian probability density
```math
w(x) = \frac{1}{\sqrt{2\pi}} \mathrm{e}^{-x^2/2}
```
is symmetric relative to the origin $m=0$.
If the weight function has any parameters, then they are stored in the dictionary `pars`.
For example, the probability density of the Beta distribution on $\Omega = [0,1]$ has two positive shape parameters $\alpha, \beta > 0$
```math
w(x) = \frac{1}{B(\alpha,\beta)} x^{\alpha-1} (1-x)^{\beta-1}.
```
[This tutorial shows the above in action.](@ref UnivariateMonicOrthogonalPolynomials)

### ProductMeasure
So far, everything was univariate, the weight of the measure was mapping real numbers to real numbers.
`PolyChaos` can handle product measures too.
Let's assume the weight function is a product of two independent Gaussian densities
```math
w: \mathbb{R} \times \mathbb{R} \rightarrow \mathbb{R}, \quad w(x) = \frac{1}{\sqrt{2\pi}} \mathrm{e}^{x_1^2/2} \frac{1}{\sqrt{2\pi}} \mathrm{e}^{x_2^2/2}.
```

The type `ProductMeasure` serves this purpose, with its straightforward fields

| Field | Meaning |
| --- | --- |
| `w::Function`  | Weight function |
| `measures::Vector{<:AbstractMeasure}`  | Vector of univariate measures |

[This tutorial shows the above in action.](@ref MultivariateMonicOrthogonalPolynomials)

### Canonical Measures
Canonical measures are special, because we know their orthogonal polynomials.
That is why several canonical measures are pre-defined in `PolyChaos`.
Some of them may require additional parameters.
(alphabetical order)

#### Beta01Measure
| Field | Meaning |
| --- | --- |
|`w::Function` |  $\frac{1}{B(\alpha,\beta)}  t^{\alpha-1} (1-t)^{\beta-1}$ | 
| `dom::Tuple{<:Real,<:Real}` | $(0, 1)$  |
| `symmetric::Bool` | true if $\alpha = \beta$ |
| `ashapeParameter::Real` | $\alpha > 0$ |
| `bshapeParameter::Real` | $\beta > 0$ |

#### GammaMeasure

| Field | Meaning |
| --- | --- |
| `w::Function` | $\frac{\beta^\alpha}{\Gamma(\alpha)} t^{\alpha-1} \exp(-\beta t)$ |
| `dom::Tuple{<:Real,<:Real}` | $(0, \infty)$|
| `symmetric::Bool` | false|
| `shapeParameter::Real` | $\alpha > 0$|
| `rateParameter::Real` | $1$|

#### GaussMeasure
| Field | Meaning |
| --- | --- |
| `w::Function` |  $\frac{1}{\sqrt{2 \pi}} \, \exp \left( - \frac{t^2}{2} \right)$|
| `dom::Tuple{<:Real,<:Real}` | $(-\infty, \infty)$|
| `symmetric::Bool` | true |

#### HermiteMeasure

| Field | Meaning |
| --- | --- |
| `w::Function` | $ \exp\left( - t^2 \right)$ |
| `dom::Tuple{<:Real,<:Real}` |$(-\infty, \infty)$ |
| `symmetric::Bool` | true |

#### JacobiMeasure

| Field | Meaning |
| --- | --- |
| `dom::Tuple{<:Real,<:Real}` | $(-1, 1)$|
| `symmetric::Bool` | true if $\alpha = \beta$ |
| `ashapeParameter::Real` | $\alpha > -1$ |
| `bshapeParameter::Real` | $\beta > -1$ |

#### LaguerreMeasure

| Field | Meaning |
| --- | --- |
| `w::Function` | $\exp(-t)$ |
| `dom::Tuple{<:Real,<:Real}` | $(0, \infty)$ |
| `symmetric::Bool` | true |

#### LegendreMeasure

| Field | Meaning |
| --- | --- |
| `w::Function` |  $1$|
| `dom::Tuple{<:Real,<:Real}` | $(-1, 1)$|
| `symmetric::Bool` | true |

#### LogisticMeasure

| Field | Meaning |
| --- | --- |
| `w::Function` | $\frac{\exp(-t)}{(1+\exp(-t))^2}$ |
| `dom::Tuple{<:Real,<:Real}` | $(-\infty, \infty)$|
| `symmetric::Bool` | true |

#### MeixnerPollaczekMeasure

| Field | Meaning |
| --- | --- |
| `w::Function` | $\frac{1}{2 \pi} \exp((2\phi-\pi)t) \lvert\Gamma(\lambda + \mathrm{i}t)\rvert^2$   |
| `dom::Tuple{<:Real,<:Real}` | $(-\infty,\infty)$|
| `symmetric::Bool` | false |
| `λParameter::Real` | $\lambda > 0$|
| `ϕParameter::Real` | $0 < \phi < \pi$|

#### Uniform01Measure

| Field | Meaning |
| --- | --- |
| `w::Function` |  $1$|
| `dom::Tuple{<:Real,<:Real}` |$(0, 1)$ |
| `symmetric::Bool` | true |

#### Uniform_11Measure

| Field | Meaning |
| --- | --- |
| `w::Function` |  $0.5$|
| `dom::Tuple{<:Real,<:Real}` |$(-1, 1)$ |
| `symmetric::Bool` | true |

#### genHermiteMeasure

| Field | Meaning |
| --- | --- |
| `w::Function` | $ \lvert t \rvert^{2 \mu}\exp \left( - t^2 \right)$ |
| `dom::Tuple{<:Real,<:Real}` | $(-\infty, \infty)$ |
| `symmetric::Bool` | true |
| `muParameter::Real` | $\mu > -0.5$ |

#### genLaguerreMeasure

| Field | Meaning |
| --- | --- |
| `w::Function` |  $t^{\alpha}\exp(-t)$  |
| `dom::Tuple{<:Real,<:Real}` | $(0,\infty)$|
| `symmetric::Bool` | false |
| `shapeParameter::Bool` | $\alpha>-1$ |



## AbstractOrthoPoly
Orthogonal polynomials are at the heart of `PolyChaos`.
The type tree for `AbstractOrthoPoly` looks as follows

```julia
julia> print_tree(AbstractOrthoPoly)
AbstractOrthoPoly
├─ AbstractCanonicalOrthoPoly
│  ├─ Beta01OrthoPoly
│  ├─ GammaOrthoPoly
│  ├─ GaussOrthoPoly
│  ├─ HermiteOrthoPoly
│  ├─ JacobiOrthoPoly
│  ├─ LaguerreOrthoPoly
│  ├─ LegendreOrthoPoly
│  ├─ LogisticOrthoPoly
│  ├─ MeixnerPollaczekOrthoPoly
│  ├─ Uniform01OrthoPoly
│  ├─ Uniform_11OrthoPoly
│  ├─ genHermiteOrthoPoly
│  └─ genLaguerreOrthoPoly
├─ MultiOrthoPoly
└─ OrthoPoly
```

It mirrors the type tree from `AbstractMeasure`: there is a generica (univariate) type `OrthoPoly`, a multivariate extension `MultiOrthoPoly` for product measures, and several univariate canonical orthogonal polynomials.

## OrthoPoly
Given an absolutely continuous measure we are wondering what are the monic polynomials $\phi_i: \Omega \rightarrow \mathbb{R}$ that are orthogonal relative to this very measure?
Mathematically this reads
```math
\langle \phi_i, \phi_j \rangle = \int_{\Omega} \phi_i(t) \phi_j(t) w(t) \mathrm{d}t =
\begin{cases}
> 0, & i=j \\
= 0, & i\neq j.
\end{cases}
```
They can be constructed using the type `OrthoPoly`, which has the fields

| Name | Meaning |
| --- | --- |
| `name::String`| Name|
| `deg::Int`         | Maximum degree|
| `α::Vector{<:Real}` | Vector of recurrence coefficients α|
| `β::Vector{<:Real}` | Vector of recurrence coefficients β |
| `meas::AbstractMeasure` | Underlying measure|


The purpose of `name` is obvious.
The integer `deg` stands for the maxium degree of the polynomials.
Rather than storing the polynomials $\phi_i$ themselves we store the recurrence coefficients `α`, `β` that characterize the system of orthogonal polynomials.
These recurrence coefficients are the single most important piece of information for the orthogonal polynomials.
For several common measures, there exist analytic formulae.
These are built-in to `PolyChaos` and should be used whenever possible.

[This tutorial shows the above in action.](@ref UnivariateMonicOrthogonalPolynomials)

## MultiOrthoPoly
Just as we did in the univariate case, we use `ProductMeasure` as a building block for multivariate orthogonal polynomials.
The type `MultiOrthoPoly` combines product measures with the respective orthogonal polynomials and their quadrature rules.
Its fields are

| Name | Meaning |
| --- | --- |
| `name::Vector{String}` | Vector of names|
| `deg::Int` | Maximum degree|
| `dim::Int` | Dimension|
| `ind::Matrix{<:Int}` | Array of multi-indices|
| `measure::ProductMeasure` | Underlying product measure|

The names of the univariate bases are stored in `names`; the maximum degree of the basis is `deg`; the overall dimension of the multivariate basis is `dim`; the multi-index `ind` maps the $j$-th multivariate basis to the elements of the univariate bases; the product measure is stored in `meas`; finally, all univariate bases are collected in `uni`.

[This tutorial shows the above in action.](@ref MultivariateMonicOrthogonalPolynomials)

## AbstractCanonicalOrthoPoly
These are the bread-and-butter polynomials: polynomials for which we know analytic formulae for the recursion coefficients.
The following canonical orthogonal polynomials are implemented

```julia
julia> print_tree(AbstractCanonicalOrthoPoly)
AbstractCanonicalOrthoPoly
├─ Beta01OrthoPoly
├─ GammaOrthoPoly
├─ GaussOrthoPoly
├─ HermiteOrthoPoly
├─ JacobiOrthoPoly
├─ LaguerreOrthoPoly
├─ LegendreOrthoPoly
├─ LogisticOrthoPoly
├─ MeixnerPollaczekOrthoPoly
├─ Uniform01OrthoPoly
├─ Uniform_11OrthoPoly
├─ genHermiteOrthoPoly
└─ genLaguerreOrthoPoly
```

Their fields follow

| Name | Meaning |
| --- | --- |
| `deg::Int` | Maximum degree | 
| `α::Vector{<:Real}` | Vector of recurrence coefficients| 
| `β::Vector{<:Real}` | Vector of recurrence coefficients| 
| `measure::CanonicalMeasure`|  Underlying canonical measure | 
| `quad::AbstractQuad` | Quadrature rule | 




## Quad
Quadrature rules are intricately related to orthogonal polynomials.
An $n$-point quadrature rule is a pair of so-called nodes $t_k$ and weights $w_k$ for $k=1,\dots,n$ that allow to solve integrals relative to the measure
```math
\int_\Omega f(t) w(t) \mathrm{d} t \approx \sum_{k=1}^n w_k f(t_k).
```
If the integrand $f$ is polynomial, then the specific Gauss quadrature rules possess the remarkable property that an $n$-point quadrature rule can integrate polynomial integrands $f$ of degree at most $2n-1$ *exactly*; no approximation error is made.

The fields of `Quad` are

| Name | Meaning |
| --- | --- |
| `name::String`        | Name |
| `Nquad::Int`          | Number $n$ of quadrature points |
| `nodes::Vector{<:Real}` | Nodes |
| `weights::Vector{<:Real}` | Weights |

with obvious meanings.

`PolyChaos` provides the type `EmptyQuad` that is added in case no quadrature rule is desired.

[This tutorial shows the above in action.](@ref NumericalIntegration)






## Tensor
The last type we need to address is `Tensor`.
It is used to store the results of scalar products.
Its fields are

| Name | Meaning |
| --- | --- |
| `dim:`| *Dimension* $m$ of tensor $\langle \phi_{i_1} \phi_{i_2} \cdots \phi_{i_{m-1}}, \phi_{i_m} \rangle$ |
| `T::SparseVector{Float64,Int}`| Entries of tensor |
| `get::Function`| Function to get entries from `T` |
| `op::AbstractOrthoPoly`| Underlying univariate orthogonal polynomials |

The *dimension* $m$ of the tensor is the number of terms that appear in the scalar product.
Let's assume we set $m = 3$, hence have $\langle \phi_{i} \phi_{j}, \phi_{k} \rangle$, then the concrete entry is obtained as `Tensor.get([i,j,k])`.

[This tutorial shows the above in action.](@ref ComputationOfScalarProducts)
