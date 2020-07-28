```@setup mysetup
using Plots
function f(x,μ,σ)
    1/sqrt(2 *π*σ^2) * exp(-(x - μ)^2 / (2σ^2))
end
μ, σ = [1., 1.7], [0.2, 0.3]
ρ(x) = 0.5*f(x,μ[1],σ[1]) + 0.5*f(x,μ[2],σ[2])
x = 0:0.01:3
plot(x,ρ.(x))
xlabel!("x")
ylabel!("rho(x)")
using PolyChaos
deg = 4
meas = Measure("my_GaussMixture", ρ, (-Inf,Inf), false, Dict(:μ=>μ, :σ=>σ)) # build measure
op = OrthoPoly("my_op", deg, meas; Nquad = 100,Nrec = 2*deg) # construct orthogonal polynomial
showbasis(op, digits=2) # in case you wondered
T2 = Tensor(2,op) # compute scalar products
T2num_1 = [ T2.get([i,j]) for i in 0:deg, j in 0:deg]
using QuadGK
T2num_2 = [quadgk(x -> evaluate(i,x,op)*evaluate(j,x,op)*ρ(x),-Inf,Inf)[1] for i in 0:deg, j in 0:deg ]
T2num_1 - T2num_2
```
# Gaussian Mixture Models
Gaussian mixture models are popular for clustering data.
Generally speaking, they are continuous random variables with a special probability density, namely
```math
\rho(x) = \sum_{i = 1}^{n} \frac{w_i}{\sqrt{2 \pi \sigma_i^2}} \exp \left( \frac{(x - \mu_i)^2}{2 \sigma_i^2} \right) \quad \text{with} \quad \sum_{i = 1}^n w_i = 1,
```
where the pairs of means and standard deviations $(\mu_i, \sigma_i)$, and the weights $w_i$ for all $i \in \{ 1, \dots, n \}$ are given.
Let's consider a simple example.


```@example mysetup
using Plots
f(x,μ,σ) = 1 / sqrt(2*π*σ^2) * exp(-(x - μ)^2 / (2σ^2))
μs, σs, ws = [1., 1.7], [0.2, 0.3], [0.5, 0.5]
ρ(x) = sum(w*f(x, μ, σ) for (μ, σ, w) in zip(μs, σs, ws))
x = 0:0.01:3;
plot(x, ρ.(x))
xlabel!("x"); ylabel!("rho(x)");
```

This looks nice!

What are now the polynomials that are orthogonal relative to this specific density?

```@example mysetup
using PolyChaos
deg = 4
meas = Measure("my_GaussMixture", ρ, (-Inf,Inf), false, Dict(:μ=>μ, :σ=>σ)) # build measure
op = OrthoPoly("my_op", deg, meas; Nquad = 100, Nrec = 2*deg) # construct orthogonal polynomial
showbasis(op, digits=2) # in case you wondered
```

Let's add the quadrature rule and compute the square norms of the basis polynomials.


```@example mysetup
T2 = Tensor(2,op) # compute scalar products
[ T2.get([i,j]) for i in 0:deg, j in 0:deg ]
```

Great!

