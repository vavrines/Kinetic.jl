# [Mathematical Background](@id MathematicalBackground)

This section is heavily based on the book "Orthogonal Polynomials: Computation and Approximation" by Walter Gautschi (Oxford University Press)

## Orthogonal Polynomials

### Basic Theory
We always work with absolutely continuous measures for which we write $\mathrm{d} \lambda (t) = w(t) \mathrm{d}t$, where the so-called *weight function* $w$
- is a nonnegative integrable function on the real line $\mathbb{R}$, i.e. $$w: \mathcal{W} \subseteq \mathbb{R} \rightarrow \mathbb{R}_{\geq 0},$$
- has finite limits in case $\mathcal{W} = \mathbb{R}$, i.e.
```math
\lim_{t \to \pm \infty} w(t) < \infty,
```
- has finite moments of all orders
```math
\mu_r(\mathrm{d}\lambda) = \int_{\mathcal{W}} t^r \mathrm{d} \lambda (t), \quad r = 0, 1, 2, \dots \quad \text{with}\: \mu_0 > 0.
```

For any pair of integrable functions $u, v$, their scalar product relative to $\mathrm{d} \lambda$ is defined as
```math
\langle u, v \rangle_{\mathrm{d} \lambda} = \int_{\mathcal{W}} u(t) v(t) \mathrm{d} \lambda(t)
```
Let $\mathcal{P}$ be the set of real polynomials and $\mathcal{P}_d \subset \mathcal{P}$ be the set of real polynomials of degree at most $d$ on $\mathcal{W}$, respectively.
Monic real polynomials are real polynomials with leading coefficient equal to one, i.e. $\pi_k(t) = t^k + \dots$ for $k = 0, 1, \dots.$

The polynomials $u,v \in \mathcal{P}$ with $u \neq v$ are said to be orthogonal if
```math
\langle u, v \rangle_{\mathrm{d} \lambda} = \int_{\mathcal{W}} u(t) v(t) \mathrm{d} \lambda(t) = 0.
```
The norm of $u$ is given by
```math
\| u \|_ {\mathrm{d}\lambda} = \sqrt{\langle u, u \rangle}.
```
If the polynomials $u \in \mathcal{P}$ has unit length $\| u \|_ {\mathrm{d}\lambda} = 1$, it is called *orthonormal*.

*Monic orthogonal polynomials* are polynomials that are monic and orthogonal, hence satisfy
- $\pi_k(t) = \pi_k(t; \mathrm{d} \lambda) = t^k + \dots$ for $k = 0, 1, \dots$, and
- $\langle \pi_k, \pi_l \rangle_{\mathrm{d}\lambda} = 0$ for $k \neq l$ and $k, l = 0, 1, \dots$, and $\langle \pi_k, \pi_k \rangle_{\mathrm{d}\lambda} = \| \pi_k \|^2_ {\mathrm{d}\lambda} > 0$ for $k = 0, 1, \dots$.

!!! note
    The support $\mathcal{W}$ of $\mathrm{d} \lambda$ can be an interval (finite, half-finite, infinite), or a finite number of disjoint intervals.
    If the support consists of a finite or denumerably infinite number of distinct points $t_k$ at which $\lambda$ has positive jumps $w_k$, the measure is called a *discrete measure*.
    For a finite number $N$ of points, the discrete measure is denoted by $\mathrm{d}\lambda_N$, and it is characterized by its nodes and weights $\{ t_k, w_k \}_{k=1}^N$ according to
    ```math
    \mathrm{d} \lambda_N (t) = \sum_{k=1}^N w_k \delta(t-t_k),
    ```
    where $\delta$ is the $\delta$-function.

    The inner product associated with $\mathrm{d} \lambda_N$ is
    ```math
    \langle u, v \rangle_{\mathrm{d}\lambda_N} = \int_{\mathcal{W}} u(t) v(t) \mathrm{d} \lambda_N (t) = \sum_{k=1}^{N} w_k u(t_k) v(t_k).
    ```
    There exist only $N$ orthogonal polynomials $\{ \pi_k(\mathrm{d} \lambda_N) \}_{k=0}^{N-1}$ that are orthogonal relative to the discrete measure $\mathrm{d} \lambda_N$ in the sense
    ```math
    \langle \pi_k(\mathrm{d} \lambda_N) \pi_l(\mathrm{d} \lambda_N) \rangle_{\mathrm{d}\lambda_N} = \| \pi_k(\mathrm{d} \lambda_N) \|_{\mathrm{d} \lambda_N} \delta_{kl},
    ```
    where $\delta_{kl}$ is the Dirac-delta, for $k,l = 0, 1, \dots, N-1$.

### Properties

#### Symmetry
An absolutely continuous measure $\mathrm{d} \lambda(t) = w(t) \mathrm{d} t$ is symmetric (with respect to the origin) if its support is $\mathcal{W} = [-a,a]$ for some $0 < a \leq \infty$, and if $w(-t) = w(t)$ for all $t \in \mathcal{W}$.

Similarly, a discrete measure $\mathrm{d} \lambda_N (t) = \sum_{k=1}^N w_k \delta(t-t_k)$ is symmetric if $t_k = - t_{N+1-k}$, and $w_k = w_{N+1-k}$ for $k=1, 2, \dots, N$.

Theorem 1.17 states that: If $\mathrm{d} \lambda$ is symmetric, then
```math
\pi_k(-t; \mathrm{d} \lambda) = (-1)^k \pi_k(t; \mathrm{d} \lambda), \quad k=0,1, \dots,
```
hence the parity of $k$ decides whether $\pi_k$ is even/odd.

!!! note
    Symmetry is exploited in [`computeSP`](@ref), where symmetry need not be relative to the origin, but some arbitrary point of the support.

### Three-term Recurrence Relation
The fact that orthogonal polynomials can be represented in terms of a three-term recurrence formula is at the heart of all numerical methods of the package.
The importance of the three-term recurrence relation is difficult to overestimate. It provides
- efficient means of evaluating polynomials (and derivatives),
- zeros of orthogonal polynomials by means of a eigenvalues of a symmetric, tridiagonal matrix
- acces to quadrature rules,
- normalization constants to create orthonormal polynomials.

Theorem 1.27 states:

Let $\pi_k(\cdot) = \pi_k(\cdot; \mathrm{d}\lambda)$ for $k = 0, 1, \dots$ be the monic orthogonal polynomials with respect to the measure $\mathrm{d} \lambda$. Then
```math
\begin{aligned}
\pi_{k+1}(t) &= (t - \alpha_k) \pi_k(t) - \beta_k \pi_{k-1}(t), \quad k= 0, 1, \dots, \\
\pi_o(t) &= 1, \\
\pi_{-1}(t) &= 0,
\end{aligned}
```
where
```math
\begin{aligned}
\alpha = \alpha_k(\mathrm{d} \lambda) &= \frac{\langle t \pi_k, \pi_k \rangle_{\mathrm{d} \lambda}}{\langle \pi_k, \pi_k \rangle_{\mathrm{d} \lambda}},  &&k=0,1,2, \dots \\
\beta = \beta_k(\mathrm{d} \lambda) &= \frac{\langle \pi_k, \pi_k \rangle_{\mathrm{d} \lambda}}{\langle \pi_{k-1}, \pi_{k-1} \rangle_{\mathrm{d} \lambda}},  &&k=1,2,\dots
\end{aligned}
```

Let $\tilde{\pi}_k(\cdot) = \tilde{\pi}_k(\cdot; \mathrm{d} \lambda t)$ denote the orthonormal polynomials, then
```math
\begin{aligned}
\sqrt{\beta_{k+1}} \tilde{\pi}_k(t) &= (t - \alpha_k) \tilde{\pi}_{k}(t) - \sqrt{\beta_k} \tilde{\pi}_{k-1}(t), \quad k = 0, 1, \dots \\
\tilde{\pi}_0(t) &= 1 \\
\tilde{\pi}_{-1}(t) &= 0.
\end{aligned}
```

!!! note
    Within the package, the coefficients `(α,β)` are *the* building block to represent (monic) orthogonal polynomials.

Notice that $\beta_0$ is arbitrary.
Nevertheless, it is convenient to define it as
```math
\beta_0(\mathrm{d}\lambda) = \langle \pi_0, \pi_0 \rangle_{\mathrm{d} \lambda} = \int_{\mathcal{W}} \mathrm{d} \lambda (t),
```
because it allows to compute the norms of the polynomials based on $\beta_k$ alone
```math
\| \pi_n \|_{\mathrm{d} \lambda} = \beta_n(\mathrm{d} \lambda) \beta_{n-1}(\mathrm{d} \lambda) \cdots \beta_0(\mathrm{d} \lambda), \quad n = 0,1, \dots
```

Let the support be $\mathcal{W} = [a,b]$ for $0 < a,b < \infty$, then
```math
\begin{aligned}
& a < \alpha_k(\mathrm{d} \lambda) < b && k = 0,1,2, \dots \\
& 0 < \beta_k(\mathrm{d} \lambda) < \max(a^2, b^2) && k = 1, 2, \dots
\end{aligned}
```

## Quadrature Rules
An $n$-point quadrature rule for the measure $\mathrm{d} \lambda t$ is a formula of the form
```math
\int_{\mathcal{W}} f(t) \mathrm{d} \lambda(t) = \sum_{\nu = 1}^{n} w_\nu f(\tau_\nu) + R_n(f).
```
The quadrature rule $\{ (\tau_\nu, w_\nu) \}_{\nu=1}^n$ composed of (mutually distinct) nodes $\tau_\nu$ and weights $w_\nu$ provides an approximation to the integral.
The respective error is given by $R_n(f)$.
If, for polynomials $p \in \mathcal{P}_d$, the error $R_n(p)$ vanishes, the respective quadrature rule is said to have a degree of exactness $d$.
Gauss quadrature rule are special quadrature rules that have a degree of exactness $d = 2n - 1$.
That means, taking a $n =3$-point quadrature rule, polynomials up to degree 5 can be integrated *exactly*.
The nodes and weights for the Gauss quadrature rules have some remarkable properties:
- all Gauss nodes are mutually distinct and contained in the interior of the support of $\mathrm{d} \lambda$;
- the $n$ Gauss nodes are the zeros of $\pi_n$, the monic orthogonal polynomial of degree $n$ relative to the measure $\mathrm{d} \lambda$;
- all Gauss weights are positive.

The Gauss nodes and weights can be computed using the [Golub-Welsch algorithm](https://en.wikipedia.org/wiki/Gaussian_quadrature#The_Golub-Welsch_algorithm).
This means to solve an eigenvalue problem of a symmetric tridiagonal matrix.
