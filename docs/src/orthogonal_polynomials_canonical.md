# [Univariate Monic Orthogonal Polynomials](@id UnivariateMonicOrthogonalPolynomials)
Univariate monic orthogonal polynomials make up the core building block of the package.
These are real polynomials $\{ \pi_k \}_{k \geq 0}$, which are univariate $\pi_k: \mathbb{R} \rightarrow \mathbb{R}$ and orthogonal relative to a nonnegative weight function $w: \mathbb{R} \rightarrow \mathbb{R}_{\geq 0}$, and which have a leading coefficient equal to one:
```math
\begin{aligned}
\pi_k(t) &= t^k + a_{k-1} t^{k-1} + \dots + a_1 t + a_0 \quad \forall k = 0, 1, \dots \\
\langle \pi_k, \pi_l \rangle &= \int_{\mathbb{R}} \pi_k(t) \pi_l(t) w(t) \mathrm{d}t =
\begin{cases}
0 & k \neq l, \text{ and }k,l \geq 0 \\
\| \pi_k \|^2 > 0 & k = l \geq 0
\end{cases}
\end{aligned}
```

These univariate monic orthogonal polynomials satisfy the paramount three-term recurrence relation
```math
\begin{aligned}
\pi_{k+1}(t) &= (t - \alpha_k) \pi_k(t) - \beta_k \pi_{k-1}(t), \quad k= 0, 1, \dots, \\
\pi_o(t) &= 1, \\
\pi_{-1}(t) &= 0.
\end{aligned}
```

Hence, every system of $n$ univariate monic orthogonal polynomials $\{ \pi_k \}_{k=0}^n$ is isomorphic to its recurrence coefficients $\{ \alpha_k, \beta_k \}_{k=0}^n$.


## Canonical Orthogonal Polynomials

The so-called *classical* or *canonical* orthogonal polynomials are polynomials named after famous mathematicians who each discovered a special family of orthogonal polynomials, for example [Hermite polynomials](https://en.wikipedia.org/wiki/Hermite_polynomials) or [Jacobi polynomials](https://en.wikipedia.org/wiki/Jacobi_polynomials).
For *classical* orthogonal polynomials there exist closed-form expressions of---among others---the recurrence coefficients.
Also quadrature rules for *classical* orthogonal polynomials are well-studied (with dedicated packages such as [FastGaussQuadrature.jl](https://github.com/ajt60gaibb/FastGaussQuadrature.jl).
However, more often than not these *classical* orthogonal polynomials are neither monic nor orthogonal, hence not normalized in any sense.
For example, there is a distinction between the [*probabilists'* Hermite polynomials](https://en.wikipedia.org/wiki/Hermite_polynomials#Definition) and the [*physicists'* Hermite polynomials](https://en.wikipedia.org/wiki/Hermite_polynomials#Definition).
The difference is in the weight function $w(t)$ relative to which the polynomials are orthogonal:
```math
\begin{aligned}
&\text{Probabilists':} &&&w(t) = \frac{1}{\sqrt{2 \pi}} \, \exp \left( - \frac{t^2}{2} \right) \\
&\text{Physicists':} &&&w(t) =  \exp \left( - t^2 \right).
\end{aligned}
```

To streamline the computations, all *classical* orthogonal polynomials are converted to __monic__ orthogonal polynomials (for which, of course, the closed-form expressions persist).
Currently, the following weight functions (hence *classical* orthogonal polynomials) are supported:

| Name | Weight $w(t)$| Parameters | Support| *Classical* polynomial |
| --- | --- | --- | --- | --- |
| `hermite` | $ \exp \left( - t^2 \right)$ | - | $(-\infty, \infty)$ | Hermite |
| `genhermite` | $ \lvert t \rvert^{2 \mu}\exp \left( - t^2 \right)$ | $\mu > -\frac{1}{2}$ | $(-\infty, \infty)$ | Generalized Hermite |
| `legendre` | $1$ | - | $(-1,1)$ | Legendre
| `jacobi` | $(1-t)^{\alpha} (1+t)^{\beta}$ | $\alpha, \beta > -1$ | $(-1,1)$ | Jacobi |
| `laguerre` | $\exp(-t)$ | - | $(0,\infty)$ | Laguerre |
| `genlaguerre` | $t^{\alpha}\exp(-t)$ | $\alpha>-1$ | $(0,\infty)$ | Generalized Laguerre |
| `meixnerpollaczek` | $\frac{1}{2 \pi} \exp((2\phi-\pi)t) \lvert\Gamma(\lambda + \mathrm{i}t)\rvert^2$ |$\lambda > 0, 0<\phi<\pi$ | $(-\infty,\infty)$ | Meixner-Pollaczek


Additionally, the following weight functions that are equivalent to probability density functions are supported:

| Name | Weight $w(t)$| Parameters | Support| *Classical* polynomial |
| --- | --- | --- | --- | --- |
| `gaussian` | $\frac{1}{\sqrt{2 \pi}} \, \exp \left( - \frac{t^2}{2} \right)$ | - | $(-\infty, \infty)$ | Probabilists' Hermite |
| `uniform01` | $1$ | - | $(0,1)$ |  Legendre
| `beta01` | $\frac{1}{B(\alpha,\beta)} \, t^{\alpha-1} (1-t)^{\beta-1}$ |$\alpha, \beta > 0$ | $(0,1)$ | Jacobi |
| `gamma` | $\frac{\beta^\alpha}{\Gamma(\alpha)} t^{\alpha-1} \exp(-\beta t)$ | $\alpha, \beta > 0$ | $(0,\infty)$ | Laguerre |
| `logistic` | $\frac{\exp(-t)}{(1+\exp(-t))^2}$ | - | $(-\infty,\infty)$ | -

To generate the orthogonal polynomials up to maximum degree `deg`, simply call


```jldoctest mylabel
julia> using PolyChaos

julia> deg = 4
4

julia> op = GaussOrthoPoly(deg)
GaussOrthoPoly{Array{Float64,1},GaussMeasure,Quad{Float64,Array{Float64,1}}}(4, [0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 2.0, 3.0, 4.0], GaussMeasure(PolyChaos.w_gaussian, (-Inf, Inf), true), Quad{Float64,Array{Float64,1}}("golubwelsch", 4, [-2.3344142183389778, -0.7419637843027257, 0.7419637843027258, 2.3344142183389778], [0.04587585476806844, 0.45412414523193134, 0.45412414523193106, 0.04587585476806852]))

julia> show(op)

Univariate orthogonal polynomials
degree:         4
#coeffs:        5
α =             [0.0, 0.0, 0.0, 0.0, 0.0]
β =             [1.0, 1.0, 2.0, 3.0, 4.0]

Measure dλ(t)=w(t)dt
w:      w_gaussian
dom:    (-Inf, Inf)
symmetric:      true

```

This generates `op`as a `GaussOrthoPoly` type with the underlying Gaussian measure `op.measure`.
The recurrence coefficients are accessible via `coeffs()`.


```jldoctest mylabel
julia> coeffs(op)
5×2 Array{Float64,2}:
 0.0  1.0
 0.0  1.0
 0.0  2.0
 0.0  3.0
 0.0  4.0
```

By default, the constructor for `OrthoPoly` generates `deg+1` recurrence coefficients.
Sometimes, some other number `Nrec` may be required.
This is why `Nrec` is a keyword for the constructor `OrthoPoly`.


```jldoctest mylabel
julia> N = 100
100

julia> opLogistic = LogisticOrthoPoly(deg; Nrec=N)
LogisticOrthoPoly{Array{Float64,1},LogisticMeasure,Quad{Float64,Array{Float64,1}}}(4, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 3.289868133696453, 10.527578027828648, 22.841084471092515, 40.10505915363294, 62.30810859273584, 89.4476035231595, 121.52266752315666, 158.53293971318436, 200.47824915030117  …  19986.565781520196, 20433.165380253333, 20884.69978120051, 21341.168984361153, 21802.572989734712, 22268.91179732067, 22740.185407118537, 23216.393819127847, 23697.53703334815, 24183.61504977904], LogisticMeasure(PolyChaos.w_logistic, (-Inf, Inf), true), Quad{Float64,Array{Float64,1}}("golubwelsch", 99, [-285.97091675697385, -266.56611354854135, -251.01698966393153, -237.53179686807928, -225.4187633699017, -214.31820469129195, -204.0126795649811, -194.35793540921836, -185.25200558110012, -176.61940782973926  …  176.61940782973895, 185.25200558110018, 194.35793540921847, 204.01267956498108, 214.31820469129212, 225.4187633699016, 237.53179686807948, 251.01698966393138, 266.56611354854135, 285.9709167569736], [1.4541663108207099e-123, 2.897917000559268e-115, 1.3858976222735606e-108, 8.826460482953542e-103, 1.4618715331286334e-97, 8.935651454381735e-93, 2.49282531464423e-88, 3.6557113389197252e-84, 3.1147999002113552e-80, 1.660700338355251e-76  …  1.6607003383554774e-76, 3.1147999002111335e-80, 3.6557113389195227e-84, 2.492825314644278e-88, 8.935651454380596e-93, 1.461871533128785e-97, 8.826460482953113e-103, 1.3858976222735651e-108, 2.8979170005595435e-115, 1.4541663108207404e-123]))

julia> show(opLogistic)

Univariate orthogonal polynomials
degree:         4
#coeffs:        100
α =             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]...
β =             [1.0, 3.289868133696453, 10.527578027828648, 22.841084471092515, 40.10505915363294, 62.30810859273584, 89.4476035231595]...

Measure dλ(t)=w(t)dt
w:      w_logistic
dom:    (-Inf, Inf)
symmetric:      true

```

Let's check whether we truly have more coefficients:


```jldoctest mylabel
julia> size(coeffs(opLogistic),1) == N
true
```

## Arbitrary Weights
If you are given a weight function $w$ that does not belong to the Table above, it is still possible to generate the respective univariate monic orthogonal polynomials.
First, we define the measure by specifying a name, the weight, the support, symmetry, and parameters


```jldoctest mylabel
julia> supp = (-1, 1)
(-1, 1)

julia> w(t) = 1 + t
w (generic function with 1 method)

julia> my_meas = Measure("my_meas", w, supp, false, Dict())
Measure("my_meas", w, (-1.0, 1.0), false, Dict{Any,Any}())
```

Notice: it is advisable to define the weight such that an error is thrown for arguments outside of the support.

Now, we want to construct the univariate monic orthogonal polynomials up to degree `deg` relative to `my_meas`.
The constructor is


```jldoctest mylabel
julia> my_op = OrthoPoly("my_op", deg, my_meas; Nquad=200);

julia> show(my_op)

Univariate orthogonal polynomials
degree:         4
#coeffs:        5
α =             [0.3333333333333335, 0.06666666666666644, 0.028571428571428848, 0.015873015873015657, 0.010101010101010171]
β =             [2.0, 0.2222222222222223, 0.23999999999999996, 0.24489795918367344, 0.2469135802469136]

Measure dλ(t)=w(t)dt
name:   my_meas
w:      w
dom:    (-1.0, 1.0)
symmetric:      false
pars:   Dict{Any,Any}()

```

By default, the recurrence coefficients are computed using the [Stieltjes procuedure](https://warwick.ac.uk/fac/sci/maths/research/grants/equip/grouplunch/1985Gautschi.pdf) with [Clenshaw-Curtis](https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature) quadrature (with `Nquad` nodes and weights).
Hence, the choice of `Nquad` influences accuracy.


## [Multivariate Monic Orthogonal Polynomials](@id MultivariateMonicOrthogonalPolynomials)
Suppose we have $p$ systems of univariate monic orthogonal polynomials,
```math
\{ \pi_k^{(1)} \}_{k\geq 0}, \: \{ \pi_k^{(2)} \}_{k\geq 0}, \dots, \{ \pi_k^{(p)} \}_{k\geq 0},
```
each system being orthogonal relative to the weights $w^{(1)}, w^{(2)}, \dots, w^{(p)}$ with supports $\mathcal{W}^{(1)}, \mathcal{W}^{(2)}, \dots, \mathcal{W}^{(p)}$.
Also, let $d^{(i)}$ be the maximum degree of the $i$-th system of univariate orthogonal polynomials.
We would like to construct a $p$-variate monic basis $\{ \psi_k \}_{k \geq 0}$ with $\psi: \mathbb{R}^p \rightarrow \mathbb{R}$ of degree at most $0 \leq d \leq \min_{i=1,\dots,k}\{ d^{(i)}\}$.
Further, this basis shall be orthogonal relative to the product measure $w: \mathcal{W} = \mathcal{W}^{(1)} \otimes \mathcal{W}^{(2)} \mathcal{W}^{(1)} \cdots \otimes \mathcal{W}^{(p)} \rightarrow \mathbb{R}_{\geq 0}$ given by
```math
w(t) = \prod_{i=1}^{p} w^{(i)}(t_i),
```
hence satisfies
```math
\langle \psi_k, \psi_l \rangle = \int_{\mathcal{W}} \psi_k(t) \psi_l(t) w(t) \mathrm{d} t =
\begin{cases}
0 & k \neq l, \text{ and }k,l \geq 0 \\
\| \psi_k \|^2 > 0 & k = l \geq 0
\end{cases}
```

For this, there exists the composite struct `MultiOrthoPoly`.
Let's consider an example where we mix *classical* orthogonal polynomials with an arbitrary weight.


```jldoctest mylabel
julia> deg = [3, 5, 6, 4]
4-element Array{Int64,1}:
 3
 5
 6
 4

julia> d = minimum(deg)
3

julia> op1 = GaussOrthoPoly(deg[1]);

julia> op2 = Uniform01OrthoPoly(deg[2]);

julia> op3 = Beta01OrthoPoly(deg[3], 2, 1.2);

julia> ops = [op1, op2, op3, my_op];

julia> mop = MultiOrthoPoly(ops, d);

julia> show(mop)

4-variate orthogonal polynomials
name:           GaussOrthoPoly{Array{Float64,1},GaussMeasure,Quad{Float64,Array{Float64,1}}}
                Uniform01OrthoPoly{Array{Float64,1},Uniform01Measure,Quad{Float64,Array{Float64,1}}}
                Beta01OrthoPoly{Array{Float64,1},Beta01Measure,Quad{Float64,Array{Float64,1}}}
                my_op
deg:            3
dim:            35
ind:            [0, 0, 0, 0]
                [1, 0, 0, 0]
                [0, 1, 0, 0]
                [0, 0, 1, 0]
                [0, 0, 0, 1]
                [2, 0, 0, 0]
                [1, 1, 0, 0]
                ...

false

```

The total number of  basis polynomials is stored in the field `dim`.
The univariate basis polynomials making up the multivariate basis are stored in the field `uni`.
The field `ind` contains the multi-index, i.e. row $i$ stores what combination of univariate polynomials makes up the $i$-th multivariate polynomial.
For example,


```jldoctest mylabel
julia> i = 11;

julia> mop.ind[i+1, :]
4-element Array{Int64,1}:
 0
 1
 0
 1
```

translates mathematically to
```math
\psi_{11}(t) = \pi_0^{(1)}(t_1) \pi_1^{(2)}(t_2) \pi_0^{(3)}(t_3) \pi_1^{(4)}(t_4).
```

Notice that there is an offset by one, because the basis counting starts at 0, but Julia is 1-indexed.
The underlying measure of `mop` is now of type `ProductMeasure`, and stored in the field `measure`
The weight $w$ can be evaluated as one would expect.
