```@setup mysetup
using PolyChaos, LinearAlgebra
my_f(t) = t^2
a, b = 1.23, 3.45 # shape parameters of Jacobi weight
int_exact = 0.353897; # reference value 
N = 4
α, β = rm_jacobi(N+1,a,b)
n_gauss, w_gauss = gauss(N,α,β)
int_gauss = dot(w_gauss, my_f.(n_gauss))
print("first point:\t $(n_gauss[1])\n")
print("end point:\t $(n_gauss[end])\n")
print("error Gauss:\t $(int_gauss - int_exact)\n")
n_radau, w_radau = radau(N-1, α, β, 1.)
int_radau = dot(w_radau, my_f.(n_radau))
print("first point:\t $(n_radau[1])\n")
print("end point:\t $(n_radau[end])\n")
print("error Radau:\t $(int_radau - int_exact)")
n_lob, w_lob = lobatto(N-2, α, β, -1., 1.)
int_lob = dot(w_lob, my_f.(n_lob))
print("first point:\t $(n_lob[1])\n")
print("end point:\t $(n_lob[end])\n")
print("error Lobatto:\t $(int_lob - int_exact)")
n_fej, w_fej = fejer(N)
int_fej = dot(w_fej, my_f.(n_fej).*(1 .- n_fej).^a.*(1 .+ n_fej).^b)
print("first point:\t $(n_fej[1])\n")
print("end point:\t $(n_fej[end])\n")
print("error Fejer:\t $(int_fej - int_exact)")
n_fej2, w_fej2 = fejer2(N)
int_fej2 = dot(w_fej2, my_f.(n_fej2).*(1 .- n_fej2).^a.*(1 .+ n_fej2).^b)
print("first point:\t $(n_fej2[1])\n")
print("end point:\t $(n_fej2[end])\n")
print("error Fejer2:\t $(int_fej2 - int_exact)")
n_cc, w_cc = clenshaw_curtis(N)
int_cc = dot(w_cc, my_f.(n_cc).*(1 .- n_cc).^a.*(1 .+ n_cc).^b)
print("first point:\t\t $(n_cc[1])\n")
print("end point:\t\t $(n_cc[end])\n")
print("error Clenshaw-Curtis:\t $(int_cc - int_exact)")
```

# [Quadrature Rules](@id QuadratureRules)
In this tutorial we investigate how recurrence coefficients of orthogonal polynomials lead to quadrature rules.

We want to solve the integral
```math
I = \int_{-1}^{1} f(t) w(t) \mathrm{d} t,
```
with the weight function
```math
w(t) = (1-t)^a (1+t)^b
```
for all $t \in [-1, 1]$ and $a, b > -1$.
For the function $f$ we choose
```math
f(t) = t^2.
```
To solve the integral we do the following:
    
0. Choose number of nodes $N$;
1. Generate recurrence coefficients;
2. Generate quadrature rule from those recurrence coefficients.

We will compare Gauss quadrature to Gauss-Radau quadrature and Gauss-Lobatto quadrature.

Make sure to check out [this tutorial](@ref NumericalIntegration) too.

Let's begin:


```@example mysetup
using PolyChaos, LinearAlgebra
my_f(t) = t^2
a, b = 1.23, 3.45 # shape parameters of Jacobi weight
int_exact = 0.353897; # reference value 
```

Now we compute $N$ recurrence coefficients.


```@example mysetup
N = 4
α, β = rm_jacobi(N+1, a, b)
```

## Gauss
The first quadrature rule is Gauss quadrature.
This method goes back to [Golub and Welsch](https://en.wikipedia.org/wiki/Gaussian_quadrature#The_Golub-Welsch_algorithm).


```@example mysetup
n_gauss, w_gauss = gauss(N, α, β)
int_gauss = dot(w_gauss, my_f.(n_gauss))
print("first point:\t $(n_gauss[1])\n")
print("end point:\t $(n_gauss[end])\n")
print("error Gauss:\t $(int_gauss - int_exact)\n")
```

Since Gauss quadrature has a degree of exactness of $2N-1$, the value of the integral is exact.
## Gauss-Radau
Gauss-Radau quadrature is a variant of Gauss quadrature that allows to specify a value of a node that *has to be included*.
We choose to include the right end point $t = 1.0$.


```@example mysetup
n_radau, w_radau = radau(N-1, α, β, 1.)
int_radau = dot(w_radau, my_f.(n_radau))
print("first point:\t $(n_radau[1])\n")
print("end point:\t $(n_radau[end])\n")
print("error Radau:\t $(int_radau - int_exact)")
```

## Gauss-Lobatto
Next, we look at Gauss-Lobatto quadrature, which allows to include two points.
We choose to include the left and end point of the interval, which are $t \in [-1.0, 1.0]$.


```@example mysetup
n_lob, w_lob = lobatto(N-2, α, β, -1., 1.)
int_lob = dot(w_lob, my_f.(n_lob))
print("first point:\t $(n_lob[1])\n")
print("end point:\t $(n_lob[end])\n")
print("error Lobatto:\t $(int_lob - int_exact)")
```

There are other quadratures that we subsume as *all-purpose* quadrature rules.
These include Fejér's first and second rule, and Clenshaw-Curtis quadrature.
## Fejér's First Rule
Fejér's first rule does *not* include the end points of the interval.


```@example mysetup
n_fej, w_fej = fejer(N)
int_fej = dot(w_fej, my_f.(n_fej).*(1 .- n_fej).^a.*(1 .+ n_fej).^b)
print("first point:\t $(n_fej[1])\n")
print("end point:\t $(n_fej[end])\n")
print("error Fejer:\t $(int_fej-int_exact)")
```

## Fejér's Second Rule
[Fejér's second](https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature) rule does include the end points of the interval.


```@example mysetup
n_fej2, w_fej2 = fejer2(N)
int_fej2 = dot(w_fej2, my_f.(n_fej2).*(1 .- n_fej2).^a.*(1 .+ n_fej2).^b)
print("first point:\t $(n_fej2[1])\n")
print("end point:\t $(n_fej2[end])\n")
print("error Fejer2:\t $(int_fej2 - int_exact)")
```

## Clenshaw-Curtis
[Clenshaw-Curtis quadrature](https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature) is similar to Féjer's second rule, as in it includes the end points of the integration interval. For the same number of nodes it is also more accurate than Féjer's rules, generally speaking.


```@example mysetup
n_cc, w_cc = clenshaw_curtis(N)
int_cc = dot(w_cc, my_f.(n_cc).*(1 .- n_cc).^a.*(1 .+ n_cc).^b)
print("first point:\t\t $(n_cc[1])\n")
print("end point:\t\t $(n_cc[end])\n")
print("error Clenshaw-Curtis:\t $(int_cc - int_exact)")
```

As we can see, for the same number of nodes $N$, the quadrature rules based on the recurrence coefficients can greatly outperform the all-purpose quadratures.
So, whenever possible, use quadrature rules based on recurrence coefficients of the orthogonal polynomials relative to the underlying measure.
Make sure to check out [this tutorial](@ref NumericalIntegration) too.
