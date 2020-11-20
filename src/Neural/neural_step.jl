"""
Collision solver of universal Boltzmann equation

    step_ube!(fwL, fhL, fbL, w, prim, h, b, fwR, fhR, fbR, u, weights, p; mode=:nn)

@notes: u & weights kept in args for multi-dispatch

"""
function step_ube!(fwL, fhL, fbL, w, prim, h, b, fwR, fhR, fbR, u, weights, p; mode=:nn)

    K, γ, μ, ω, Pr, dx, dt, RES, AVG = p[1:9]
    ann = p[10:end]

    #--- record W^{n} ---#
    w_old = deepcopy(w)
    h_old = deepcopy(h)
    b_old = deepcopy(b)

    H_old = maxwellian(u, prim)
    B_old = H_old .* K ./ (2. .* prim[end])
    τ_old = vhs_collision_time(prim, μ, ω)

    #--- update W^{n+1} ---#
    @. w += (fwL - fwR) / dx
    prim .= conserve_prim(w, γ)

    H = maxwellian(u, prim)
    B = H .* K ./ (2. .* prim[end])
    τ = vhs_collision_time(prim, μ, ω)

    #--- update f^{n+1} ---#
    for i in eachindex(h)
        h[i] += (fhL[i] - fhR[i]) / dx
        b[i] += (fbL[i] - fbR[i]) / dx
    end

    if mode == :bgk
        @. h = (h + H / τ * dt) / (1. + dt / τ)
        @. b = (b + B / τ * dt) / (1. + dt / τ)
    elseif mode == :shakhov
        qf = heat_flux(h, b, prim, u, weights)
        H1, B1 = shakhov(u, H, B, qf, prim, Pr, K)
        H .+= H1
        B .+= B1

        @. h = (h + H / τ * dt) / (1. + dt / τ)
        @. b = (b + B / τ * dt) / (1. + dt / τ)
    elseif mode == :nn
        df = ube_dfdt([h_old; b_old], ([H_old; B_old], τ_old, ann), dt)

        @. h += df[1:end÷2] * dt
        @. b += df[end÷2+1:end] *dt
    end

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

end