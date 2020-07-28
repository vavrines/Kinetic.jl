```@setup mysetup
using PolyChaos, LinearAlgebra
deg, n = 4, 20
s_α, s_β = 2.1, 3.2
opq = Beta01OrthoPoly(deg, s_α, s_β; Nrec=n, addQuadrature=true)
normsq = computeSP2(opq)
m = 3
t = Tensor(3,opq)
t.get([1,2,3])
T = [ t.get([i1,i2,i3]) for i1=0:dim(opq)-1, i2=0:dim(opq)-1, i3=0:dim(opq)-1]
#@show normsq == diag(T[:,:,1])
#@show normsq == diag(T[:,1,:])
#@show normsq == diag(T[1,:,:])
t2 = Tensor(2, opq)
@show normsq == [ t2.get([i, i]) for i in 0:dim(opq)-1]
using SpecialFunctions
supp = (0, 1)
w(t) = (t^(s_α-1)*(1-t)^(s_β-1)/SpecialFunctions.beta(s_α,s_β))
my_meas = Measure("my_meas", w, supp, false)
my_opq = OrthoPoly("my_op", deg, my_meas; Nrec=n, addQuadrature = true)
my_normsq = computeSP2(my_opq)
my_t = Tensor(m, my_opq)
my_T = [ my_t.get([i1,i2,i3]) for i1=0:dim(opq)-1,i2=0:dim(opq)-1,i3=0:dim(opq)-1]
@show abs.(normsq-my_normsq)
@show norm(T-my_T)
mop = MultiOrthoPoly([opq, my_opq], deg)
mt2 = Tensor(2,mop)
mt3 = Tensor(3,mop)
mT2 = [ mt2.get([i,i]) for i=0:dim(mop)-1 ]
mop.ind
ind_opq = findUnivariateIndices(1,mop.ind)
ind_my_opq = findUnivariateIndices(2,mop.ind)
@show mT2[ind_opq] - normsq
@show mT2[ind_my_opq] - my_normsq;
```
# [Computation of Scalar Products](@id ComputationOfScalarProducts)
By now, we are able to construct orthogonal polynomials, and to construct quadrature rules for a given nonnegative weight function, respectively.
Now we combine both ideas to solve integrals involving the orthogonal polynomials
```math
\langle \phi_{i_1} \phi_{i_2} \cdots \phi_{i_{m-1}}, \phi_{i_m} \rangle
= \int \phi_{i_1}(t) \phi_{i_2}(t) \cdots \phi_{i_{m-1}}(t) \phi_{i_m}(t) w(t) \mathrm{d} t,
```
both for the univariate and multivariate case.
The integrand is a polynomial (possibly multivariate) that can be solved exactly with the appropriate Gauss quadrature rules.

!!! note
    To simplify notation we drop the integration interval.
    It is clear from the context.


## Univariate Polynomials
### Classical Polynomials
Let's begin with a univariate basis for some *classical* orthogonal polynomial


```@example mysetup
using PolyChaos
deg, n = 4, 20
s_α, s_β = 2.1, 3.2
opq = Beta01OrthoPoly(deg, s_α, s_β; Nrec=n, addQuadrature=true)
```

By setting `addQuadrature = true` (which is default), an $n$-point Gauss quadrature rule is create relative to the underlying measure `opq.measure`, where $n$ is the number of recurrence coefficients stored in `opq.α` and `opq.β`.

To compute the squared norms
```math
\| \phi_k \|^2 = \langle \phi_k, \phi_k  \rangle
= \int \phi_k(t) \phi_k(t) w(t) \mathrm{d} t
```

of the basis we call `computeSP2()`


```@example mysetup
normsq = computeSP2(opq)
```

For the general case
```math
\langle \phi_{i_1} \phi_{i_2} \cdots \phi_{i_{m-1}}, \phi_{i_m} \rangle
= \int \phi_{i_1}(t) \phi_{i_2}(t) \cdots \phi_{i_{m-1}}(t) \phi_{i_m}(t) w(t) \mathrm{d} t,
```
there exists a type `Tensor` that requires only two arguments: the *dimension* $m \geq 1$, and an `AbstractOrthoPoly`


```@example mysetup
m = 3
t = Tensor(3,opq)
```

To get the desired entries, `Tensor` comes with a `get()` function that is called for some index $a \in \mathbb{N}_0^m$ that has the entries $a = [i_1, i_2, \dots, i_m]$.
For example



```@example mysetup
t.get([1,2,3])
```

Or using comprehension


```@example mysetup
T = [ t.get([i1,i2,i3]) for i1=0:dim(opq)-1, i2=0:dim(opq)-1, i3=0:dim(opq)-1]
```

Notice that we can cross-check the results.


```@example mysetup
using LinearAlgebra
normsq == diag(T[:,:,1]) == diag(T[:,1,:]) == diag(T[1,:,:])
```

Also, `normsq` can be computed analogously in `Tensor` format


```@example mysetup
t2 = Tensor(2, opq)
normsq == [ t2.get([i, i]) for i in 0:dim(opq)-1]
```

### Arbitrary Weights
Of course, the type `OrthoPoly` can be constructed for arbitrary weights $w(t)$.
In this case we have to compute the orthogonal basis and the respective quadrature rule.
Let's re-work the above example by hand.


```@example mysetup
using SpecialFunctions
supp = (0, 1)
w(t) = (t^(s_α-1)*(1-t)^(s_β-1)/SpecialFunctions.beta(s_α,s_β))
my_meas = Measure("my_meas", w, supp, false)
my_opq = OrthoPoly("my_op", deg, my_meas; Nrec=n, addQuadrature = true)
```

Now we can compute the squared norms $\| \phi_k \|^2$


```@example mysetup
my_normsq = computeSP2(my_opq)
```

And the tensor


```@example mysetup
my_t = Tensor(m, my_opq)
my_T = [ my_t.get([i1,i2,i3]) for i1=0:dim(opq)-1,i2=0:dim(opq)-1,i3=0:dim(opq)-1]
```

Let's compare the results:

```@example mysetup
abs.(normsq-my_normsq)
```

```@example mysetup
norm(T-my_T)
```

!!! note
    The possibility to create quadrature rules for arbitrary weights should be reserved to cases different from *classical* ones.

## Multivariate Polynomials
For multivariate polynomials the syntax for `Tensor` is very much alike, except that we are dealing with the type `MultiOrthoPoly` now.


```@example mysetup
mop = MultiOrthoPoly([opq, my_opq], deg)
```


```@example mysetup
mt2 = Tensor(2,mop)
mt3 = Tensor(3,mop)
mT2 = [ mt2.get([i,i]) for i=0:dim(mop)-1 ]
```

Notice that `mT2` carries the elements of the 2-dimensional tensors for the univariate bases `opq` and `my_opq`.
The encoding is given by the multi-index `mop.ind`


```@example mysetup
mop.ind
```

To cross-check the results we can distribute the multi-index back to its univariate indices with the help of `findUnivariateIndices`.


```@example mysetup
ind_opq = findUnivariateIndices(1,mop.ind)
ind_my_opq = findUnivariateIndices(2,mop.ind)
```

```@example mysetup
mT2[ind_opq] - normsq
```

```@example mysetup
mT2[ind_my_opq] - my_normsq
```
