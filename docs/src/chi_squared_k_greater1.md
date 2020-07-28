```@setup mysetup
k = 12
using PolyChaos
degree, Nrec = 2, 20
opq = GaussOrthoPoly(degree; Nrec=Nrec, addQuadrature = true);
mop = MultiOrthoPoly([opq for i in 1:k], degree)
L = dim(mop)
mu, sig = 0., 1.
x = [ assign2multi(convert2affinePCE(mu, sig, opq),i,mop.ind) for i in 1:k ]
t2 = Tensor(2,mop)
t3 = Tensor(3,mop)
y = [ sum( x[i][j1]*x[i][j2]*t3.get([j1-1,j2-1,m-1])/t2.get([m-1,m-1])  for i=1:k, j1=1:L, j2=1:L ) for m=1:L ]
moms_analytic(k) = [k, sqrt(2k), sqrt(8/k)]
function myskew(y)
   e3 = sum( y[i]*y[j]*y[k]*t3.get([i-1,j-1,k-1]) for i=1:L,j=1:L,k=1:L )
   μ = y[1]
   σ = std(y,mop)
   (e3-3*μ*σ^2-μ^3)/(σ^3)
end
print("Expected value:\t\t$(moms_analytic(k)[1]) = $(mean(y,mop))\n")
print("\t\t\terror = $(abs(mean(y,mop)-moms_analytic(k)[1]))\n")
print("Standard deviation:\t$(moms_analytic(k)[2]) = $(std(y,mop))\n")
print("\t\t\terror = $(moms_analytic(k)[2]-std(y,mop))\n")
print("Skewness:\t\t$(moms_analytic(k)[3]) = $(myskew(y))\n")
print("\t\t\terror = $(moms_analytic(k)[3]-myskew(y))\n")
using Plots
Nsmpl = 10000
ysmpl = samplePCE(Nsmpl, y, mop)
histogram(ysmpl;normalize=true, xlabel="t",ylabel="rho(t)")
import SpecialFunctions: gamma
ρ(t) = 1/(2^(0.5*k)*gamma(0.5*k))*t^(0.5*k-1)*exp(-0.5*t)
t = range(0.1; stop=maximum(ysmpl), length=100)
plot!(t,ρ.(t),w=4)
```

# Chi-squared Distribution ($k>1$)


## Theory
Given $k$ standard random variables $X_i \sim \mathcal{N}(0,1)$ for $i=1,\dots,k$ we would like to find the random variable $Y = \sum_{i=1}^k X_i^2$.
The analytic solution is known: $Y$ follows a chi-squared distribution with $k$ degrees of freedom.

Using polynomial chaos expansion (PCE), the problem can be solved using Galerkin projection.
Let $\{\phi_k \}_{k=0}^{n}$ be the monic orthogonal basis relative to the probability density of $X = [X_1, \dots, X_k]$, namely
```math
f_X(x) =  \prod_{i=1}^k \frac{1}{\sqrt{2 \pi}} \, \exp \left( - \frac{x_i^2}{2} \right).
```
Then, the PCE of $X_i$ is given by
```math
X_i = \sum_{k=0}^n x_{i,k} \phi_k,
```
with
```math
x_{i,0} = 0, \quad x_{i,i+1} = 1, \quad x_{i,l} = 0 \quad \forall l \in \{1,\dots,n\} \setminus \{i+1\}.
```
To find the PCE coefficients $y_k$ for $Y = \sum_{k=}^n y_k \phi_k$, we apply Galerkin projection, which leads to
```math
y_m \langle \phi_m, \phi_m \rangle = \sum_{i=1}^k \sum_{j_1=0}^n \sum_{j_2=0}^n x_{i,j_1} x_{i,j_2} \langle \phi_{j_1} \phi_{j_2}, \phi_m \rangle \quad \forall m = 0, \dots, n.
```
Hence, knowing the scalars $\langle \phi_m, \phi_m \rangle$, and $\langle \phi_{j_1} \phi_{j_2}, \phi_m \rangle$, the PCE coefficients $y_k$ can be obtained immediately.
From the PCE coefficients we can get the moments and compare them to the closed-form expressions.

__Notice:__ A maximum degree of 2 suffices to get the *exact* solution with PCE.
In other words, increasing the maximum degree to values greater than 2 introduces nothing but computational overhead (and numerical errors, possibly).


## Practice
First, we create a orthogonal basis relative to $f_X(x)$ of degree at most $d=2$ (`degree` below).

Notice that we consider a total of `Nrec` recursion coefficients, and that we also add a quadrature rule by setting `addQuadrature = true`.


```@example mysetup
k = 12
using PolyChaos
degree, Nrec = 2, 20
opq = GaussOrthoPoly(degree; Nrec=Nrec, addQuadrature = true);
```

Now let's define a multivariate basis

```@example mysetup
mop = MultiOrthoPoly([opq for i in 1:k], degree)
```

Next, we define the PCE for all $X_i$ with $i = 1, \dots, k$.


```@example mysetup
L = dim(mop)
mu, sig = 0., 1.
x = [ assign2multi(convert2affinePCE(mu, sig, opq), i, mop.ind) for i in 1:k ]
```

With the orthogonal basis and the quadrature at hand, we can compute the tensors `t2` and `t3` that store the entries $\langle \phi_m, \phi_m \rangle$, and $\langle \phi_{j_1} \phi_{j_2}, \phi_m \rangle$, respectively.


```@example mysetup
t2 = Tensor(2,mop)
t3 = Tensor(3,mop)
```

With the tensors at hand, we can compute the Galerkin projection.

Notice: there are more efficient ways to do this, but let's keep it simple.


```@example mysetup
y = [ sum( x[i][j1]*x[i][j2]*t3.get([j1-1,j2-1,m-1])/t2.get([m-1,m-1])  for i=1:k, j1=1:L, j2=1:L ) for m=1:L ]
```

Let's compare the moments via PCE to the closed-form expressions.


```@example mysetup
moms_analytic(k) = [k, sqrt(2k), sqrt(8/k)]
function myskew(y)
   e3 = sum( y[i]*y[j]*y[k]*t3.get([i-1,j-1,k-1]) for i=1:L,j=1:L,k=1:L )
   μ = y[1]
   σ = std(y,mop)
   (e3-3*μ*σ^2-μ^3)/(σ^3)
end

print("Expected value:\t\t$(moms_analytic(k)[1]) = $(mean(y,mop))\n")
print("\t\t\terror = $(abs(mean(y,mop)-moms_analytic(k)[1]))\n")
print("Standard deviation:\t$(moms_analytic(k)[2]) = $(std(y,mop))\n")
print("\t\t\terror = $(moms_analytic(k)[2]-std(y,mop))\n")
print("Skewness:\t\t$(moms_analytic(k)[3]) = $(myskew(y))\n")
print("\t\t\terror = $(moms_analytic(k)[3]-myskew(y))\n")

```

Let's plot the probability density function to compare results.
We first draw samples from the measure with the help of `sampleMeasure()`, and then evaluate the basis at these samples and multiply times the PCE coefficients.
The latter stop is done using `evaluatePCE()`.
Both steps are combined in the function `samplePCE()`.
Finally, we compare the result agains the analytical PDF $\rho(t) = \frac{t^{t/2-1}\mathrm{e}^{-t/2}}{2^{k/2} \, \Gamma(k/2)}$ of the chi-squared distribution with one degree of freedom.


```@example mysetup
using Plots
Nsmpl = 10000
# long way: ξ = sampleMeasure(Nsmpl,mop), ysmpl = evaluatePCE(y,ξ,mop)
ysmpl = samplePCE(Nsmpl, y, mop)
histogram(ysmpl;normalize=true, xlabel="t",ylabel="rho(t)")

import SpecialFunctions: gamma
ρ(t) = 1/(2^(0.5*k)*gamma(0.5*k))*t^(0.5*k-1)*exp(-0.5*t)
t = range(0.1; stop=maximum(ysmpl), length=100)
plot!(t, ρ.(t), w=4)
```
