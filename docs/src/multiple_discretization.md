```@setup mysetup
using PolyChaos, LinearAlgebra, Plots
γ = 0.5;
int_exact = 1+pi/2; # exact value of the integral
function my_w(t, γ)
    γ + (1 - γ) * 1 / sqrt(1 - t^2)
end
N = 1000;
n,w = fejer(N);
int_fejer = dot(w,my_w.(n,γ))
print("Fejer error:\t$(abs(int_exact-int_fejer))\twith $N nodes")
function quad_gaussleg(N,γ)
    a, b = rm_legendre(N)
    n, w = golubwelsch(a,b)
end
n, w = quad_gaussleg(N+1, γ)
int_gaussleg = dot(w,γ .+ (1-γ)/sqrt.(1 .- n.^2))
print("Gauss-Legendre error:\t$(abs(int_exact-int_gaussleg))\twith $N nodes")
function quad_gausscheb(N,γ)
    a, b = rm_chebyshev1(N)
    n, w = golubwelsch(a, b)
end
n, w = quad_gausscheb(N+1,γ)
int_gausscheb = dot(w,γ .+ (1-γ)*sqrt.(1 .- n.^2))
print("Gauss-Chebyshev error:\t$(abs(int_exact-int_gausscheb))\twith $(length(n)) nodes")
function quad_gaussleg_mod(N::Int,γ::Float64)
    n,w = quad_gaussleg(N+1,γ)
    return n,γ*w
end
function quad_gausscheb_mod(N::Int,γ::Float64)
            n,w = quad_gausscheb(N+1,γ)
    return n,(1-γ)*w
end
N = 8
a,b = mcdiscretization(N,[n->quad_gaussleg_mod(n,γ); n->quad_gausscheb_mod(n,γ)])
n,w = golubwelsch(a,b)
int_mc = sum(w)
print("Discretization error:\t$(abs(int_exact-int_mc))\twith $(length(n)) nodes")
Γ = 0:0.1:1;
ab = [ mcdiscretization(N,[n->quad_gaussleg_mod(n,gam); n->quad_gausscheb_mod(n,gam)]) for gam in Γ ];
bb = hcat([ ab[i][2] for i=1:length(Γ)]...);
b_leg = rm_legendre(N)[2];
b_cheb = rm_chebyshev1(N)[2]
bb[:,1]-b_cheb
bb[:,end]-b_leg
using Plots
plot(Γ,bb',yaxis=:log10, w=3, legend=false)
zs, os = zeros(N), ones(N)
scatter!(zs,b_cheb,marker=:x)
scatter!(os,b_leg,marker=:circle)
```
# Multiple Discretization

This tutorial shows how to compute recurrence coefficients for non-trivial weight functions, and how they are being used for quadrature.
The method we use is called *multiple discretization*, and follows W. Gautschi's book "Orthogonal Polynomials: Computation and Approximation", specifically Section 2.2.4, and Example 2.38.

Suppose we have the weight function
```math
\forall t \in [-1,1], \gamma \in [0,1]: \quad w(t;\gamma) = \gamma + (1-\gamma) \frac{1}{\sqrt{1-t^2}},
```
and we would like to solve
```math
\int_{-1}^{1} f(t) w(t;c) \mathrm{d}t = \sum_{\nu=1}^{N} f(\tau_\nu) w_\nu
```
by some quadrature rule.
We will see that ad-hoc quadrature rules will fail to solve the integral even for the simplest choice $f \equiv 1$.
However, finding the recurrence coefficients of the underlying orthogonal polynomials, and then finding the quadrature rule will do just fine.

Let us first try to solve the integral for $f \equiv 1$ by Féjer's rule.


```@example mysetup
using PolyChaos, LinearAlgebra
γ = 0.5;
int_exact = 1 + pi / 2; # exact value of the integral
function my_w(t, γ)
    γ + (1 - γ) * 1 / sqrt(1 - t^2)
end

N = 1000;
nodes, weights = fejer(N);
int_fejer = dot(weights, my_w.(nodes, γ))
print("Fejer error:\t$(abs(int_exact - int_fejer))\twith $N nodes")
```

Clearly, that is not satisfying.
Well, the term $\gamma$ of the weight $w$ makes us think of Gauss-Legendre integration, so let's try it instead.


```@example mysetup
function quad_gaussleg(N, γ)
    a, b = rm_legendre(N)
    nodes, weights = golubwelsch(a, b)
end
nodes, weights = quad_gaussleg(N+1, γ)
int_gaussleg = dot(weights, γ .+ (1-γ)/sqrt.(1 .- nodes.^2))
print("Gauss-Legendre error:\t$(abs(int_exact-int_gaussleg))\twith $N nodes")
```

Even worse!
Well, we can factor out $\frac{1}{\sqrt{1-t^2}}$, making the integral amenable to a Gauss-Chebyshev rule.
So, let's give it anothery try.


```@example mysetup
function quad_gausscheb(N, γ)
    a, b = rm_chebyshev1(N)
    nodes, weights = golubwelsch(a, b)
end
nodes, weights = quad_gausscheb(N+1, γ)
int_gausscheb = dot(weights, γ .+ (1-γ)*sqrt.(1 .- nodes.^2))
# int=sum(xw(:,2).*(1+sqrt(1-xw(:,1).^2)))
print("Gauss-Chebyshev error:\t$(abs(int_exact - int_gausscheb))\twith $(length(n)) nodes")
```

Okay, that's better, but it took us a lot of nodes to get this result.
Is there a different way?
Indeed, there is.
As we have noticed, the weight $w$ has a lot in common with Gauss-Legendre *and* Gauss-Chebyshev.
We can decompose the integral as follows
```math
\int_{-1}^1 f(t) w(t) \mathrm{d}t = \sum_{i=1}^{m} \int_{-1}^{1} f(t) w_i(t) \mathrm{d} t,
```
with
```math
\begin{align*}
w_1(t) &= \gamma \\
w_2(t) &= (1-\gamma) \frac{1}{\sqrt{1-t^2}}.
\end{align*}
```
To the weight $w_1$ we can apply Gauss-Legendre quadrature, to the weight $w_2$ we can apply Gauss-Chebyshev quadrature (with tiny modifications).
This *discretization* of the measure can be used in our favor.
The function `mcdiscretization()` takes the $m$ discretization rules as an input


```@example mysetup
function quad_gaussleg_mod(N, γ)
    nodes, weights = quad_gaussleg(N + 1, γ)
    nodes, γ*weights
end
function quad_gausscheb_mod(N, γ)
            nodes, weights = quad_gausscheb(N + 1,γ)
    return nodes, (1-γ)*weights
end

N = 8
a, b = mcdiscretization(N, [n -> quad_gaussleg_mod(n, γ); n -> quad_gausscheb_mod(n, γ)])
nodes, weights = golubwelsch(a, b)
int_mc = sum(w)
print("Discretization error:\t$(abs(int_exact-int_mc))\twith $(length(n)) nodes")
```

Et voilà, no error with fewer nodes.
(For this example, we'd need in fact just a single node.)

The function `mcdiscretization()` is able to construct the recurrence coefficients of the orthogonal polynomials relative to the weight $w$.
Let's inspect the values of the recurrence coefficients a little more.
For $\gamma = 0$, we are in the world of Chebyshev polynomials, for $\gamma = 1$, we enter the realm of Legendre polynomials. And in between?
That's exactly where the weight $w$ comes in: it can be thought of as an interpolatory weight, interpolating Legendre polynomials and Chebyshev polynomials.
Let's verify this by plotting the recurrence coefficients for several values of $\gamma$.




```@example mysetup
Γ = 0:0.1:1;
ab = [ mcdiscretization(N, [n -> quad_gaussleg_mod(n, gam); n -> quad_gausscheb_mod(n, gam)]) for gam in Γ ];
bb = hcat([ab[i][2] for i in 1:length(Γ)]...);
b_leg = rm_legendre(N)[2];
b_cheb = rm_chebyshev1(N)[2]
bb[:,1]-b_cheb
```


```@example mysetup
bb[:,end] - b_leg
```

Let's plot these values to get a better feeling.


```@example mysetup
using Plots
plot(Γ, bb', yaxis=:log10, w=3, legend=false)
zs, os = zeros(N), ones(N)
scatter!(zs, b_cheb, marker=:x)
scatter!(os, b_leg, marker=:circle)

xlabel!("Gamma")
ylabel!("Beta")
```

The crosses denote the values of the β recursion coefficients for Chebyshev polynomials; the circles the β recursion coefficients for Legendre polynomials.
The interpolating line in between stands for the β recursion coefficients of $w(t; \gamma)$.
