```@setup mysetup
using PolyChaos
d = 6
op_gauss = GaussOrthoPoly(d, addQuadrature = true)
points, degrees = randn(10), 0:2:d
[ evaluate(degree, points, op_gauss) for degree in degrees ]
μ, σ = 2., 0.2
pce_gauss = convert2affinePCE(μ, σ, op_gauss)
a, b = -0.3, 1.2
op_uniform = Uniform01OrthoPoly(d)
convert2affinePCE(a, b, op_uniform)
pce_uniform = convert2affinePCE(μ, σ, op_uniform; kind="μσ")
α, β = 1.2, 3.4
op_beta = Beta01OrthoPoly(d, α, β)
convert2affinePCE(a, b, op_beta)
pce_beta = convert2affinePCE(μ, σ, op_beta; kind="μσ")
a1, a2 = μ, sqrt(3)*σ/pi
op_logistic = LogisticOrthoPoly(d)
pce_logistic = convert2affinePCE(a1, a2, op_logistic)
mean(pce_gauss, op_gauss), std(pce_gauss, op_gauss)
mean(pce_uniform, op_uniform), std(pce_uniform, op_uniform)
mean(pce_beta, op_beta), std(pce_beta, op_beta)
mean(pce_logistic, op_beta), std(pce_logistic, op_beta)
using Statistics
N = 1000
ξ_gauss = sampleMeasure(N, op_gauss)
samples_gauss = evaluatePCE(pce_gauss, ξ_gauss, op_gauss)
# samplePCE(N,pce_gauss,myops["gauss"])
ξ_uniform = sampleMeasure(N, op_uniform)
samples_uniform = evaluatePCE(pce_uniform, ξ_uniform, op_uniform)
# samples_uniform = samplePCE(N,pce_uniform, op_uniform)
ξ_beta = sampleMeasure(N, op_beta)
samples_beta = evaluatePCE(pce_beta, ξ_beta, op_beta)
# samples_beta = samplePCE(N, pce_beta, op_beta)
ξ_logistic = sampleMeasure(N, op_logistic)
samples_logistic = evaluatePCE(pce_logistic, ξ_logistic, op_logistic)
# samples_logistic = samplePCE(N, pce_logistic, op_logistic)
```

# [Common Random Variables](@id CommonRandomVariables)
Polynomial chaos expansion (PCE) is a Hilbert space technique for random variables with finite variance.
Mathematically equivalent to Fourier series expansions for periodic signals, PCE allows to characterize a random variable in terms of its PCE coefficients (aka Fourier coefficients).
That is, the PCE of a random variable $\mathsf{x}$ is given by
```math
\mathsf{x} = \sum_{i=0}^L x_i \phi_i,
```
where $x_i$ are the so-called PCE coefficients, and $\phi_i$ are the orthogonal polynomials that are orthogonal relative to the probability density function of $\mathsf{x}$.

This tutorial walks you through the PCE of common random variables, namely Gaussian, Beta, Uniform, Logistic, and shows how they are implemented in `PolyChaos`.

## Construction of Basis
Let's begin with the usual suspect, a Gaussian random variable.
We construct the respective orthogonal polynomials of degree at most `d`
```@example mysetup
using PolyChaos
d = 6
op_gauss = GaussOrthoPoly(d, addQuadrature = true)
```
A quadrature rule is added by setting they keyword `addQuadrature = true` (which it is by default).

Let's evaluate the Gaussian basis polynomials at some points for varying degrees.
```@example mysetup
points, degrees = randn(10), 0:2:d
[ evaluate(degree, points, op_gauss) for degree in degrees ]
```

## Finding PCE Coefficients
Having constructed the orthogonal bases, the question remains how to find the PCE coefficients for the common random variables.
Every random variable can be characterized exactly by two PCE coefficients.
For a Gaussian random variable, this is familiar: the mean and the variance suffice to describe a Gaussian random variable entirely.
The same is true for any random variable of finite variance given the right basis.
The function `convert2affinePCE` provides the first two PCE coefficients (hence the name affine) for the common random variables.

### Gaussian
Given the Gaussian random variable $\mathsf{x} \sim \mathcal{N}(\mu, \sigma^2)$ with $\sigma > 0$, the affine PCE coefficients are


```@example mysetup
μ, σ = 2., 0.2
pce_gauss = convert2affinePCE(μ, σ, op_gauss)
```

### Uniform
Given the uniform random variable $\mathsf{x} \sim \mathcal{U}(a, b)$ with finite support $a<b$, the affine PCE coefficients are

```@example mysetup
a, b = -0.3, 1.2
op_uniform = Uniform01OrthoPoly(d)
convert2affinePCE(a, b, op_uniform)
```

Instead, if the expected value and standard deviation are known, the affine PCE coefficients of the uniform random variable are

```@example mysetup
pce_uniform = convert2affinePCE(μ, σ, op_uniform; kind="μσ")
```
Notice that the zero-order coefficient *is* equivalent to the expected value μ.

### Beta
Given the Beta random variable $\mathsf{x} \sim \mathcal{B}(a, b, \alpha, \beta)$ with finite support $a<b$ and shape parameters $\alpha, \beta > 0$, the affine PCE coefficients are


```@example mysetup
α, β = 1.2, 3.4
op_beta = Beta01OrthoPoly(d, α, β)
convert2affinePCE(a, b, op_beta)
```

Instead, if the expected value and standard deviation are known, the affine PCE coefficients of the uniform random variable are


```@example mysetup
pce_beta = convert2affinePCE(μ, σ, op_beta; kind="μσ")
```

### Logistic

Given the logstic random variable $\mathsf{x} \sim \mathcal{L}(a_1,a_2)$ where $a_2>0$ with the probability density function
```math
\rho(t) = \frac{1}{4 a_2} \, \operatorname{sech}^2 \left(\frac{t-a_1}{2a_2}\right)
```
the affine PCE coefficients of the uniform random variable are


```@example mysetup
a1, a2 = μ, sqrt(3)*σ/pi
op_logistic = LogisticOrthoPoly(d)
pce_logistic = convert2affinePCE(a1, a2, op_logistic)
```

## Moments
It is a key feature of PCE to compute moments from the PCE coefficients alone; no sampling is required.

### Gaussian

```@example mysetup
mean(pce_gauss, op_gauss), std(pce_gauss, op_gauss)
```

### Uniform

```@example mysetup
mean(pce_uniform, op_uniform), std(pce_uniform, op_uniform)
```

### Beta

```@example mysetup
mean(pce_beta, op_beta), std(pce_beta, op_beta)
```

### Logistic

```@example mysetup
mean(pce_logistic, op_beta), std(pce_logistic, op_beta)
```

## Sampling
Having found the PCE coefficients, it may be useful to sample the random variables.
That means, find $N$ realizations of the random variable that obey the random variable's probability density function.
This is done in two steps:
1. Draw $N$ samples from the measure (`sampleMeasure()`), and then
2. Evaluate the basis polynomials and multiply times the PCE coefficients, i.e. $\sum_{i=0}^L x_i \phi_i(\xi_j)$ where $\xi_j$ is the $j$-th sample from the measure (`evaluatePCE()`).

Both steps are combined in the function `samplePCE()`.

### Gaussian


```@example mysetup
N = 1000
ξ_gauss = sampleMeasure(N, op_gauss)
samples_gauss = evaluatePCE(pce_gauss, ξ_gauss, op_gauss);
# samplePCE(N,pce_gauss,myops["gauss"])
```

### Uniform

```@example mysetup
ξ_uniform = sampleMeasure(N, op_uniform)
samples_uniform = evaluatePCE(pce_uniform, ξ_uniform, op_uniform);
# samples_uniform = samplePCE(N,pce_uniform, op_uniform)
```

### Beta

```@example mysetup
ξ_beta = sampleMeasure(N, op_beta)
samples_beta = evaluatePCE(pce_beta, ξ_beta, op_beta);
# samples_beta = samplePCE(N, pce_beta, op_beta)
```

### Logistic

```@example mysetup
ξ_logistic = sampleMeasure(N, op_logistic)
samples_logistic = evaluatePCE(pce_logistic, ξ_logistic, op_logistic);
# samples_logistic = samplePCE(N, pce_logistic, op_logistic)
```
