```@setup mysetup
using PolyChaos
k = 1
deg, Nrec = 2, 20
opq = GaussOrthoPoly(deg; Nrec=Nrec, addQuadrature=true);
showbasis(opq; sym="ξ") # works for `op` too!
showpoly(0:2:deg,opq)
L = dim(opq)
mu, sig = 0., 1.
x = [ convert2affinePCE(mu, sig, opq); zeros(Float64,L-2) ]
t2 = Tensor(2, opq);
t3 = Tensor(3, opq)
y = [ sum( x[i]*x[j]*t3.get([i-1,j-1,m-1])/t2.get([m-1,m-1])  for i=1:L, j=1:L ) for m=1:L ]
moms_analytic(k) = [k, sqrt(2k), sqrt(8/k)]
function myskew(y)
   e3 = sum( y[i]*y[j]*y[k]*t3.get([i-1,j-1,k-1]) for i=1:L,j=1:L,k=1:L )
   μ = y[1]
   σ = std(y,opq)
   (e3-3*μ*σ^2-μ^3)/(σ^3)
end
print("Expected value:\t\t$(moms_analytic(k)[1]) = $(mean(y,opq))\n")
print("\t\t\terror = $(abs(mean(y,opq)-moms_analytic(k)[1]))\n")
print("Standard deviation:\t$(moms_analytic(k)[2]) = $(std(y,opq))\n")
print("\t\t\terror = $(moms_analytic(k)[2]-std(y,opq))\n")
print("Skewness:\t\t$(moms_analytic(k)[3]) = $(myskew(y))\n")
print("\t\t\terror = $(moms_analytic(k)[3]-myskew(y))\n")
using Plots
Nsmpl = 10000
ysmpl = samplePCE(Nsmpl, y, opq)
histogram(ysmpl;normalize=true,xlabel="t",ylabel="rho(t)")
import SpecialFunctions: gamma
ρ(t) = 1/(sqrt(2)*gamma(0.5))*1/sqrt(t)*exp(-0.5*t)
t = range(0.1; stop=maximum(ysmpl), length=100)
plot!(t, ρ.(t), w=4)
```

# Chi-squared Distribution ($k=1$)


## Theory
Given a standard random variable $X \sim \mathcal{N}(0,1)$ we would like to find the random variable $Y = X^2$.
The analytic solution is known: $Y$ follows a chi-squared distribution with $k=1$ degree of freedom.

Using polynomial chaos expansion (PCE), the problem can be solved using Galerkin projection.
Let $\{\phi_k \}_{k=0}^{n}$ be the monic orthogonal basis relative to the probability density of $X$, namely
```math
f_X(x) = \frac{1}{\sqrt{2 \pi}} \exp \left( - \frac{x^2}{2} \right).
```
Then, the PCE of $X$ is given by
```math
X = \sum_{k=0}^n x_k \phi_k,
```
with
```math
x_0 = 0, \quad x_1 = 1, \quad x_i = 0 \quad \forall i =2,\dots,n.
```
To find the PCE coefficients $y_k$ for $Y = \sum_{k=}^n y_k \phi_k$, we apply Galerkin projection, which leads to
```math
y_m \langle \phi_m, \phi_m \rangle = \sum_{i=0}^n \sum_{j=0}^n x_i x_j \langle \phi_i \phi_j, \phi_m \rangle \quad \forall m = 0, \dots, n.
```
Hence, knowing the scalars $\langle \phi_m, \phi_m \rangle$, and $\langle \phi_i \phi_j, \phi_m \rangle$, the PCE coefficients $y_k$ can be obtained immediately.
From the PCE coefficients we can get the moments and compare them to the closed-form expressions.

__Notice:__ A maximum degree of 2 suffices to get the *exact* solution with PCE.
In other words, increasing the maximum degree to values greater than 2 introduces nothing but computational overhead (and numerical errors, possibly).


## Practice
First, we create a orthogonal basis relative to $f_X(x)$ of degree at most $d=2$ (`deg` below).

Notice that we consider a total of `Nrec` recursion coefficients, and that we also add a quadrature rule by setting `addQuadrature = true`.


```@example mysetup
using PolyChaos
k = 1
deg, Nrec = 2, 20
opq = GaussOrthoPoly(deg; Nrec=Nrec, addQuadrature=true);
```

What are the basis polynomials?

```@example mysetup
showbasis(opq; sym="ξ")
```

Note that the command `showbasis` is based on the more general `showpoly`:

```@example mysetup
showpoly(0:2:deg,opq)
```

Next, we define the PCE for $X$.


```@example mysetup
L = dim(opq)
mu, sig = 0., 1.
x = [ convert2affinePCE(mu, sig, opq); zeros(Float64,L-2) ]
```

With the orthogonal basis and the quadrature at hand, we can compute the tensors `t2` and `t3` that store the entries $\langle \phi_m, \phi_m \rangle$, and $\langle \phi_i \phi_j, \phi_m \rangle$, respectively.


```@example mysetup
t2 = Tensor(2, opq);
t3 = Tensor(3, opq)
```

With the tensors at hand, we can compute the Galerkin projection.


```@example mysetup
y = [ sum( x[i]*x[j]*t3.get([i-1,j-1,m-1])/t2.get([m-1,m-1])  for i=1:L, j=1:L ) for m=1:L ]
```

Let's compare the moments via PCE to the closed-form expressions.

```@example mysetup
moms_analytic(k) = [k, sqrt(2k), sqrt(8/k)]
function myskew(y)
   e3 = sum( y[i]*y[j]*y[k]*t3.get([i-1,j-1,k-1]) for i=1:L, j=1:L, k=1:L )
   μ = y[1]
   σ = std(y,opq)
   (e3-3*μ*σ^2-μ^3)/(σ^3)
end

print("Expected value:\t\t$(moms_analytic(k)[1]) = $(mean(y,opq))\n")
print("\t\t\terror = $(abs(mean(y,opq)-moms_analytic(k)[1]))\n")
print("Standard deviation:\t$(moms_analytic(k)[2]) = $(std(y,opq))\n")
print("\t\t\terror = $(moms_analytic(k)[2]-std(y,opq))\n")
print("Skewness:\t\t$(moms_analytic(k)[3]) = $(myskew(y))\n")
print("\t\t\terror = $(moms_analytic(k)[3]-myskew(y))\n")

```

Let's plot the probability density function to compare results.
We first draw samples from the measure with the help of `sampleMeasure()`, and then evaluate the basis at these samples and multiply times the PCE coefficients.
The latter stop is done using `evaluatePCE()`.
Finally, we compare the result agains the analytical PDF $\rho(t) = \frac{\mathrm{e}^{-0.5t}}{\sqrt{2 t} \, \Gamma(0.5)}$ of the chi-squared distribution with one degree of freedom.


```@example mysetup
using Plots
Nsmpl = 10000
# long way: ξ = sampleMeasure(Nsmpl,opq), ysmpl = evaluatePCE(y,ξ,opq)
ysmpl = samplePCE(Nsmpl, y, opq)
histogram(ysmpl; normalize=true, xlabel="t", ylabel="rho(t)")

import SpecialFunctions: gamma
ρ(t) = 1/(sqrt(2)*gamma(0.5))*1/sqrt(t)*exp(-0.5*t)
t = range(0.1; stop=maximum(ysmpl), length=100)
plot!(t, ρ.(t), w=4)
```
