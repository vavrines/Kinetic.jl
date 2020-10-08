# ============================================================
# Hydrodynamic Fluxes
# http://www.cfdbooks.com/
# ============================================================

"""
Lax-Friedrichs flux

`flux_lax!(fw::AbstractArray{<:Real,1}, wL::AbstractArray{<:Real,1}, wR::AbstractArray{<:Real,1}, γ::Real, dt::Real, dx::Real)`

_P. D. Lax, Weak Solutions of Nonlinear Hyperbolic Equations and Their Numerical Computation, 
Commun. Pure and Applied Mathematics, 7, 159-193, 1954._

"""
function flux_lax!(
    fw::AbstractArray{<:Real,1},
    wL::AbstractArray{<:Real,1},
    wR::AbstractArray{<:Real,1},
    γ::Real,
    dt::Real,
    dx::Real,
)
    fw .= 0.5 * dt .* (euler_flux(wL, γ)[1] + euler_flux(wR, γ)[1] - dx / dt .* (wR - wL))
    return nothing
end


"""
HLL flux for the Euler equations

`flux_hll!(fw::AbstractArray{<:Real,1}, wL::AbstractArray{<:Real,1}, wR::AbstractArray{<:Real,1}, γ::Real, dt::Real)`

* @args: variables at left & right sides of interface
* @args: specific heat ratio

"""
function flux_hll!(
    fw::AbstractArray{<:Real,1},
    wL::AbstractArray{<:Real,1},
    wR::AbstractArray{<:Real,1},
    γ::Real,
    dt::Real,
)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    aL = sound_speed(primL, γ)
    aR = sound_speed(primR, γ)

    λmin = primL[2] - aL
    λmax = primR[2] + aR

    if λmin >= 0.0
        fw .= euler_flux(wL, γ)
    elseif λmax <= 0.0
        fw .= euler_flux(wR, γ)
    else
        factor = 1.0 / (λmax - λmin)

        flux1 = euler_flux(wL, γ)
        flux2 = euler_flux(wR, γ)

        fw .= factor * (λmax * flux1 - λmin * flux2 + λmax * λmin * (wR - wL))
    end

    fw .*= dt

    return nothing

end


"""
Roe's flux with entropy fix

`flux_roe!(fw::AbstractArray{<:Real,1}, wL::AbstractArray{<:Real,1}, wR::AbstractArray{<:Real,1},
    γ::Real, dt::Real, n = [1.0, 0.0]::AbstractArray{<:Real,1})`

_P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference Schemes, Journal of Computational Physics, 43, pp. 357-372._
(_cf. http://cfdbooks.com/cfdcodes.html_)

* @args primL[1:4] = left state (rhoL, uL, vL, pL)
* @args primR[1:4] = right state (rhoR, uR, vR, pR)
* @args γ: specific heat ratio
* @args n[2]: unit face normal (L -> R)

"""
function flux_roe!(
    fw::AbstractArray{<:Real,1},
    wL::AbstractArray{<:Real,1},
    wR::AbstractArray{<:Real,1},
    γ::Real,
    dt::Real,
    n = [1.0, 0.0]::AbstractArray{<:Real,1},
)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    if length(fw) == 3 # 1D

        # left state
        rhoL = primL[1]
        vL = primL[2]
        pL = 0.5 * primL[1] / primL[end]
        aL = sound_speed(primL, γ)
        HL = (wL[end] + pL) / rhoL

        # right state
        rhoR = primR[1]
        vR = primR[2]
        pR = 0.5 * primR[1] / primR[end]
        aR = sound_speed(primR, γ)
        HR = (wR[end] + pR) / rhoR

        # compute Roe averages
        RT = sqrt(rhoR / rhoL)
        rho = RT * rhoL
        v = (vL + RT * vR) / (1.0 + RT)
        H = (HL + RT * HR) / (1.0 + RT)
        a = sqrt((γ - 1.0) * (H - 0.5 * v^2))

        # differences in primitive variables
        drho = rhoR - rhoL
        du = vR - vL
        dP = pR - pL

        # wave strength (characteristic variables)
        dV = [
            0.5 * (dP - rho * a * du) / a^2,
            -(dP / a^2 - drho),
            0.5 * (dP + rho * a * du) / a^2,
        ]

        # absolute values of the wave speeds (eigenvalues)
        ws = [abs(v - a), abs(v), abs(v + a)]

        # modified wave speeds for nonlinear fields (to remove expansion shocks).
        Da = max(0.0, 4.0 * ((vR - aR) - (vL - aL)))
        if ws[1] < 0.5 * Da
            ws[1] = ws[1]^2 / Da + 0.25 * Da
        end
        Da = max(0.0, 4.0 * ((vR + aR) - (vL + aL)))
        if ws[3] < 0.5 * Da
            ws[3] = ws[3]^2 / Da + 0.25 * Da
        end

        # right eigenvectors
        R = zeros(3, 3)
        R[1, 1] = 1.0
        R[2, 1] = v - a
        R[3, 1] = H - v * a

        R[1, 2] = 1.0
        R[2, 2] = v
        R[3, 2] = 0.5 * v * v

        R[1, 3] = 1.0
        R[2, 3] = v + a
        R[3, 3] = H + v * a

        # compute average flux
        fw .= 0.5 .* (euler_flux(wL, γ)[1] + physical_flux(wR)[1])

        # add matrix dissipation term to complete Roe flux
        for j = 1:3, k = 1:3
            fw[j] -= 0.5 * ws[k] * dV[k] * R[j, k]
        end

    elseif length(fw) == 4 # 2D

        # normal vector
        nx = n[1]
        ny = n[2]

        # tangent vector
        mx = -ny
        my = nx

        #--- primitive and other variables ---#
        # left state
        rhoL = primL[1]
        uL = primL[2]
        vL = primL[3]
        unL = uL * nx + vL * ny
        umL = uL * mx + vL * my
        pL = 0.5 * primL[1] / primL[4]
        aL = sound_speed(primL[4], γ)
        HL = aL^2 / (γ - 1.0) + 0.5 * (uL^2 + vL^2)

        # right state
        rhoR = primR[1]
        uR = primR[2]
        vR = primR[3]
        unR = uR * nx + vR * ny
        umR = uR * mx + vR * my
        pR = 0.5 * primR[1] / primR[4]
        aR = sound_speed(primR[4], γ)
        HR = aR^2 / (γ - 1.0) + 0.5 * (uR^2 + vR^2)

        # Roe averages
        RT = sqrt(rhoR / rhoL)
        rho = RT * rhoL
        u = (uL + RT * uR) / (1.0 + RT)
        v = (vL + RT * vR) / (1.0 + RT)
        H = (HL + RT * HR) / (1.0 + RT)
        a = sqrt((γ - 1.0) * (H - 0.5 * (u^2 + v^2)))
        un = u * nx + v * ny
        um = u * mx + v * my

        # wave strengths
        drho = rhoR - rhoL
        dp = pR - pL
        dun = unR - unL
        dum = umR - umL

        LdU = [
            (dp - rho * a * dun) / (2.0 * a^2),
            rho * dum,
            drho - dp / (a^2),
            (dp + rho * a * dun) / (2.0 * a^2),
        ]

        # wave speed
        ws = abs.([un - a, un, un, un + a])

        # Harten's entropy fix JCP(1983), 49, pp357-393
        # only for the nonlinear fields.
        if ws[1] < 0.2
            ws[1] = 0.5 * (ws[1]^2 / 0.2 + 0.2)
        end
        if ws[4] < 0.2
            ws[4] = 0.5 * (ws[4]^2 / 0.2 + 0.2)
        end

        # right eigenvectors
        Rv = zeros(4, 4)

        Rv[1, 1] = 1.0
        Rv[2, 1] = u - a * nx
        Rv[3, 1] = v - a * ny
        Rv[4, 1] = H - un * a

        Rv[1, 2] = 0.0
        Rv[2, 2] = mx
        Rv[3, 2] = my
        Rv[4, 2] = um

        Rv[1, 3] = 1.0
        Rv[2, 3] = u
        Rv[3, 3] = v
        Rv[4, 3] = 0.5 * (u * u + v * v)

        Rv[1, 4] = 1.0
        Rv[2, 4] = u + a * nx
        Rv[3, 4] = v + a * ny
        Rv[4, 4] = H + un * a

        # dissipation term
        diss = zeros(4)
        for i = 1:4, j = 1:4
            diss[i] += ws[j] * LdU[j] * Rv[i, j]
        end

        # compute fluxes
        fL = [
            rhoL * unL,
            rhoL * unL * uL + pL * nx,
            rhoL * unL * vL + pL * ny,
            rhoL * unL * HL,
        ]
        fR = [
            rhoR * unR,
            rhoR * unR * uR + pR * nx,
            rhoR * unR * vR + pR * ny,
            rhoR * unR * HR,
        ]

        @. fw = 0.5 * dt * (fL + fR - diss)

    end

    return nothing
    
end
