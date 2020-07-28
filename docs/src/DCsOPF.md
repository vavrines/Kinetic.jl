```@setup mysetup
using PolyChaos, JuMP, MosekTools, LinearAlgebra
A = [ -1 1 0 0; -1 0 1 0; -1 0 0 1 ; 0 1 -1 0; 0 0 -1 1] # incidence matrix
Nl, N = size(A,1), size(A,2)
Bbr = diagm(0 => -( 2 .+ 10*rand(Nl) )) # line parameters
Ψ = [ zeros(Nl)  -Bbr*A[:,2:end]*inv(A[:,2:end]'*Bbr*A[:,2:end]) ] # PTDF matrix
Cp, Cd = [1 0; 0 0; 0 0; 0 1], [0 0; 1 0; 0 1; 0 0 ] # book-keeping
Ng, Nd = size(Cp,2), size(Cd,2)
c = 4 .+ 10*rand(Ng) # cost function parameters
λp, λl = 1.6*ones(Ng), 1.6*ones(Nl) # lambdas for chance constraint reformulations
pmax, pmin = 10*ones(Ng), zeros(Ng) # engineering limits
plmax, plmin = 10*ones(Nl), -10*ones(Nl) # engineering limits
deg = 1
opq = [Uniform01OrthoPoly(deg; Nrec=5*deg), Uniform01OrthoPoly(deg; Nrec=5*deg)]
mop = MultiOrthoPoly(opq, deg)
Npce = mop.dim
d = zeros(Nd,Npce) # PCE coefficients of load
d[1,[1,2]] = convert2affinePCE(1., 0.1, mop.uni[1], kind="μσ")
d[2,[1,3]] = convert2affinePCE(2., 0.2, mop.uni[2], kind="μσ")
function buildSOC(x::Vector,mop::MultiOrthoPoly)
    t = [ sqrt(Tensor(2,mop).get([i,i])) for i in 0:mop.dim-1 ]
    (t.*x)[2:end]
end
model = Model(with_optimizer(Mosek.Optimizer))
@variable(model, p[i in 1:Ng,j in 1:Npce], base_name="p")
@constraint(model, energy_balance[j in 1:Npce], sum(p[i,j] for i in 1:Ng) - sum(d[i,j] for i in 1:Nd) == 0)
@constraint(model, con_pmax[i in 1:Ng], [1/λp[i]*(pmax[i] - mean(p[i,:],mop)); buildSOC(p[i,:],mop)] in SecondOrderCone())
@constraint(model, con_pmin[i in 1:Ng], [1/λp[i]*(mean(p[i,:],mop) - pmin[i]); buildSOC(p[i,:],mop)] in SecondOrderCone())
pl = Ψ*(Cp*p + Cd*d)
@constraint(model, con_plmax[i in 1:Nl], [1/λl[i]*(plmax[i] - mean(pl[1,:],mop)); buildSOC(pl[i,:],mop)] in SecondOrderCone())
@constraint(model, con_plmin[i in 1:Nl], [1/λl[i]*(mean(pl[1,:],mop) - plmin[i]); buildSOC(pl[i,:],mop)] in SecondOrderCone())
@objective(model, Min, sum( mean(p[i,:],mop)*c[i] for i in 1:Ng) )
optimize!(model) # here we go
@assert termination_status(model) == MOI.OPTIMAL "Model not solved to optimality."
psol, plsol, obj = value.(p), value.(pl), objective_value(model)
p_moments = [ [mean(psol[i,:],mop) var(psol[i,:],mop) ] for i in 1:Ng ]
pbr_moments = [ [mean(plsol[i,:],mop) var(plsol[i,:],mop) ] for i in 1:Nl ]
```

# Chance-Constrained DC Optimal Power Flow
The purpose of this tutorial is to show how polynomial chaos can be leveraged to solve optimization problems under uncertainty.
Specifically, we study chance-constrained DC optimal power flow as it is presented in [this paper](https://www.sciencedirect.com/science/article/pii/S235246771830105X).

We consider the following 4-bus system that has a total of two generators (buses 1 and 3) and two loads (buses 2 and 4):

![4-bus system](assets/DCsOPF_drawing.png)

We formalize the numbering of the generators (superscript $g$), loads (superscript $d$ for demand), and branches (superscript $br$) as follows
```math
\mathcal{N}^g = \{ 1, 3\}, \, \mathcal{N}^d = \{ 2, 4\}, \, \mathcal{N}^{br} = \{ 1, 2, 3, 4, 5 \}.
```
With each generator we associate a linear cost with cost coefficient $c_i$ for all $i \in \mathcal{N}^g$.
Each generator must adhere to its engineering limits given by $(\underline{p}_i^g , \overline{p}_i^g )$ for all $i \in \mathcal{N}^g$.
Also, each line is constrained by its limits $(\underline{p}_i^{br}, \overline{p}_i^{br})$ for all $i \in \mathcal{N}^{br}$.

We model the demand at the buses $i \in \mathcal{N}^d$ in terms of uniform distributions with known mean $\mu_i$ and standard deviation $\sigma_i$.
We concisely write
```math
\mathsf{p}_i^d \sim \mathsf{U}(\mu_i, \sigma_i) \quad \forall i \in \mathcal{N}^d.
```
For simplicity we consider DC conditions.
Hence, energy balance reads
```math
\sum_{i \in \mathcal{N}^g} \mathsf{p}_i^g - \sum_{i \in \mathcal{N}^d} \mathsf{p}_i^d = 0,
```
and the vector of branch flows is computed from the power transfer distribution factor (PTDF) matrix $\Psi$ via
```math
\mathsf{p}^{br} = \Psi (C^p\mathsf{p}^g + C^d\mathsf{p}^d).
```
The matrices $C^p$ and $C^d$ map the generators and the loads to the correct buses, respectively.

We want to solve a chance-constrained optimal power flow problem under DC conditions.
According to [this paper](https://www.sciencedirect.com/science/article/pii/S235246771830105X), we can formulate the problem as
$$\underset{\mathsf{p^{g}}}{\operatorname{min}} \, \sum_{i \in \mathcal{N}_g} c_i \mathbb{E}( \mathsf{p}_i^g)$$
subject to
```math
\sum_{i \in \mathcal{N}^g} \mathsf{p}_i^g - \sum_{i \in \mathcal{N}^d} \mathsf{p}_i^d = 0, \\
\underline{p}_i^g \leq \mathbb{E}(\mathsf{p}_i^g) \pm \lambda_i^g \sqrt{\mathbb{V}(\mathsf{p}_i^g)} \leq \overline{p}_i^g  \forall i \in \mathcal{N}^g,\\
\underline{p}_i^{br} \leq \mathbb{E}(\mathsf{p}_i^{br}) \pm \lambda_i^{br} \sqrt{\mathbb{V}(\mathsf{p}_i^{br})} \leq \overline{p}_i^{br} \forall i \in \mathcal{N}^{br},
```
which minimizes the total expected generation cost subject to the DC power flow equations and chance-constrained engineering limits.

Let's solve the problem using `PolyChaos` and `JuMP`, using `Mosek` as a solver.


```@example mysetup
using PolyChaos, JuMP, MosekTools, LinearAlgebra
```

Let's define system-specific quantities such as the incidence matrix and the branch flow parameters.
From these we can compute the PTDF matrix $\Psi$ (assuming the slack is at bus 1).


```@example mysetup
A = [ -1 1 0 0; -1 0 1 0; -1 0 0 1 ; 0 1 -1 0; 0 0 -1 1] # incidence matrix
Nl, N = size(A,1), size(A,2)
Bbr = diagm(0 => -( 2 .+ 10*rand(Nl) )) # line parameters
Ψ = [ zeros(Nl)  -Bbr*A[:,2:end]*inv(A[:,2:end]'*Bbr*A[:,2:end]) ] # PTDF matrix
```

Now we can continue the remaining ingredients that specify our systems:


```@example mysetup
Cp, Cd = [1 0; 0 0; 0 0; 0 1], [0 0; 1 0; 0 1; 0 0 ] # book-keeping
Ng, Nd = size(Cp,2), size(Cd,2)
c = 4 .+ 10*rand(Ng) # cost function parameters
λp, λl = 1.6*ones(Ng), 1.6*ones(Nl) # lambdas for chance constraint reformulations
pmax, pmin = 10*ones(Ng), zeros(Ng) # engineering limits
plmax, plmin = 10*ones(Nl), -10*ones(Nl) # engineering limits
```

We specify the uncertainty using `PolyChaos`:


```@example mysetup
deg = 1
opq = [Uniform01OrthoPoly(deg; Nrec=5*deg), Uniform01OrthoPoly(deg; Nrec=5*deg)]
mop = MultiOrthoPoly(opq, deg)
Npce = mop.dim
```

It remains to specify the PCE coefficients, for which we will use `convert2affine`.


```@example mysetup
d = zeros(Nd,Npce) # PCE coefficients of load
d[1,[1,2]] = convert2affinePCE(1., 0.1, mop.uni[1], kind="μσ")
d[2,[1,3]] = convert2affinePCE(2., 0.2, mop.uni[2], kind="μσ")
```

Now, let's put it all into an optimization problem, specifically a second-order cone program.
To build the second-order cone constraints we define a helper function `buildSOC`.


```@example mysetup
function buildSOC(x::Vector,mop::MultiOrthoPoly)
    t = [ sqrt(Tensor(2,mop).get([i,i])) for i in 0:mop.dim-1 ]
    (t.*x)[2:end]
end
```

Finally, let's use `JuMP` to formulate and then solve the problem.
We use `Mosek` to solve the problem; for academic use there are [free license](https://www.mosek.com/products/academic-licenses/).


```@example mysetup
model = Model(with_optimizer(Mosek.Optimizer))
@variable(model, p[i in 1:Ng,j in 1:Npce], base_name="p")
@constraint(model, energy_balance[j in 1:Npce], sum(p[i,j] for i in 1:Ng) - sum(d[i,j] for i in 1:Nd) == 0)
@constraint(model, con_pmax[i in 1:Ng], [1/λp[i]*(pmax[i] - mean(p[i,:],mop)); buildSOC(p[i,:],mop)] in SecondOrderCone())
@constraint(model, con_pmin[i in 1:Ng], [1/λp[i]*(mean(p[i,:],mop) - pmin[i]); buildSOC(p[i,:],mop)] in SecondOrderCone())
pl = Ψ*(Cp*p + Cd*d)
@constraint(model, con_plmax[i in 1:Nl], [1/λl[i]*(plmax[i] - mean(pl[1,:],mop)); buildSOC(pl[i,:],mop)] in SecondOrderCone())
@constraint(model, con_plmin[i in 1:Nl], [1/λl[i]*(mean(pl[1,:],mop) - plmin[i]); buildSOC(pl[i,:],mop)] in SecondOrderCone())
@objective(model, Min, sum( mean(p[i,:],mop)*c[i] for i in 1:Ng) )
optimize!(model) # here we go
```

Let's extract the numerical values of the optimal solution.

```@example mysetup
@assert termination_status(model) == MOI.OPTIMAL "Model not solved to optimality."
psol, plsol, obj = value.(p), value.(pl), objective_value(model)
```

Great, we've solved the problem.
How do we now make sense of the solution?
For instance, we can look at the moments of the generated power:

```@example mysetup
p_moments = [ [mean(psol[i,:],mop) var(psol[i,:],mop) ] for i in 1:Ng ]
```

Simiarly, we can study the moments for the branch flows:

```@example mysetup
pbr_moments = [ [mean(plsol[i,:],mop) var(plsol[i,:],mop) ] for i in 1:Nl ]
```
