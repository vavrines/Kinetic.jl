using KitBase, OrdinaryDiffEq, Plots

set = config_ntuple(; u0=-8, u1=8, nu=80, t1=8, nt=30, Kn=1)

tspan = (0, set.t1)
tsteps = linspace(tspan[1], tspan[2], set.nt)
γ = 3.0
vs = VSpace1D(set.u0, set.u1, set.nu)

f0 = @. 0.5 * (1 / π)^0.5 * (exp.(-(vs.u - 2)^2) + 0.5 * exp(-(vs.u + 2)^2))
prim0 = conserve_prim(moments_conserve(f0, vs.u, vs.weights), γ)
M0 = maxwellian(vs.u, prim0)

mu0 = ref_vhs_vis(set.Kn, set.α, set.ω)
τ0 = mu0 * 2.0 * prim0[end]^(0.5) / prim0[1]

# BGK
prob1 = ODEProblem(bgk_ode!, f0, tspan, [M0, τ0])
sol_bgk = solve(prob1, Tsit5(); saveat=tsteps)

# Shakhov
q = heat_flux(f0, prim0, vs.u, vs.weights)
S0 = shakhov(vs.u, M0, q, prim0, 2 / 3)
prob2 = ODEProblem(bgk_ode!, f0, tspan, [M0 .+ S0, τ0])
sol_shakhov = solve(prob2, Tsit5(); saveat=tsteps)

# ES-BGK
prob3 = ODEProblem(esbgk_ode!, f0, tspan, [vs.u, vs.weights, prim0, 2 / 3, τ0])
sol_es = solve(prob3, Tsit5(); saveat=tsteps)

idx = (rand() * set.nt |> round |> Int) + 1
begin
    plot(vs.u, sol_bgk[idx]; label="BGK")
    plot!(vs.u, sol_shakhov[idx]; label="Shakhov", line=:dash)
    scatter!(vs.u, sol_es[idx]; label="ES", alpha=0.5)
end
