# ============================================================
# Macroscopic Fluxes
# http://www.cfdbooks.com/
# ============================================================


"""
Lax-Friedrichs flux
* P. D. Lax, Weak Solutions of Nonlinear Hyperbolic Equations and Their Numerical Computation, Commun. Pure and Applied Mathematics, 7, 159-193, 1954.

"""
function flux_lax!(fw::AbstractArray{<:Real,1}, wL::AbstractArray{<:Real,1}, wR::AbstractArray{<:Real,1}, γ::Real, dt::Real, dx::Real)
    fw .= 0.5 * dt .* (euler_flux(wL, γ)[1] + euler_flux(wR, γ)[1] - dx / dt .* (wR - wL))
end


"""
HLL flux for the Euler equations

* @arg: variables at left & right sides of interface
* @arg: specific heat ratio

"""
function flux_hll!(fw::AbstractArray{<:Real,1}, wL::AbstractArray{<:Real,1}, wR::AbstractArray{<:Real,1}, γ::Real, dt::Real)

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

end

"""
Roe's flux with entropy fix

* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference Schemes, Journal of Computational Physics, 43, pp. 357-372.
* http://cfdbooks.com/cfdcodes.html

* @arg primL[1:4] = left state (rhoL, uL, vL, pL)
* @arg primR[1:4] = right state (rhoR, uR, vR, pR)
* @arg γ: specific heat ratio
* @arg n[2]: unit face normal (L -> R)

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
        RT = sqrt(rhoR/rhoL);
        rho = RT*rhoL
        v = (vL+RT*vR)/(1.0+RT)
        H = (HL+RT*HR)/(1.0+RT)
        a = sqrt( (γ-1.0)*(H-0.5*v^2) )

        # differences in primitive variables
        drho = rhoR - rhoL
        du =   vR - vL
        dP =   pR - pL

        # wave strength (characteristic variables)
        dV = [0.5*(dP-rho*a*du)/a^2, -( dP/a^2 - drho ), 0.5*(dP+rho*a*du)/a^2]

        # absolute values of the wave speeds (eigenvalues)
        ws = [abs(v-a), abs(v  ), abs(v+a)]

        # modified wave speeds for nonlinear fields (to remove expansion shocks).
        Da = max(0.0, 4.0*((vR-aR)-(vL-aL)) )
        if ws[1] < 0.5*Da
            ws[1] = ws[1]^2/Da + 0.25*Da
        end
        Da = max(0.0, 4.0*((vR+aR)-(vL+aL)) )
        if ws[3] < 0.5*Da
            ws[3] = ws[3]^2/Da + 0.25*Da
        end

        # right eigenvectors
        R = zeros(3, 3)
        R[1,1] = 1.0
        R[2,1] = v - a
        R[3,1] = H - v*a

        R[1,2] = 1.0
        R[2,2] = v
        R[3,2] = 0.5*v*v

        R[1,3] = 1.0
        R[2,3] = v + a
        R[3,3] = H + v*a

        # compute average flux
        fw .= 0.5.*( euler_flux(wL, γ)[1] + physical_flux(wR)[1] )

        # add matrix dissipation term to complete Roe flux
        for j in 1:3, k in 1:3
            fw[j] -= 0.5*ws[k]*dV[k]*R[j,k] 
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
        fL = [rhoL * unL, rhoL * unL * uL + pL * nx, rhoL * unL * vL + pL * ny, rhoL * unL * HL]
        fR = [rhoR * unR, rhoR * unR * uR + pR * nx, rhoR * unR * vR + pR * ny, rhoR * unR * HR]

        @. fw = 0.5 * dt * (fL + fR - diss)

    end

end


"""
Wave propagation method for the Maxwell's equations

* @param[in]: variables in left-left, left, right, and right-right cells
* @param[in]: eigenmatrix (A), eigenvalue (D)
* @param[in]: speed of light (sol)
* @param[in]: auxiliary parameters (χₑ, νᵦ)

"""
function flux_em!(
    femL::AbstractArray{<:Real,1},
    femR::AbstractArray{<:Real,1},
    ELL::AbstractArray{<:Real,1},
    BLL::AbstractArray{<:Real,1},
    EL::AbstractArray{<:Real,1},
    BL::AbstractArray{<:Real,1},
    ER::AbstractArray{<:Real,1},
    BR::AbstractArray{<:Real,1},
    ERR::AbstractArray{<:Real,1},
    BRR::AbstractArray{<:Real,1},
    ϕL::Real,
    ϕR::Real,
    ψL::Real,
    ψR::Real,
    dxL::Real,
    dxR::Real,
    A1p::AbstractArray{<:Real,2},
    A1n::AbstractArray{<:Real,2},
    D1::AbstractArray{<:Real,1},
    sol::Real,
    χ::Real,
    ν::Real,
    dt::Real,
)

    @assert length(femL) == length(femR) == 8

    slop = zeros(8, 8)
    slop[3, 1] = -0.5 * sol^2 * (BR[2] - BL[2]) + 0.5 * sol * (ER[3] - EL[3])
    slop[5, 1] = 0.5 * sol * (BR[2] - BL[2]) - 0.5 * (ER[3] - EL[3])
    slop[2, 2] = 0.5 * sol^2 * (BR[3] - BL[3]) + 0.5 * sol * (ER[2] - EL[2])
    slop[6, 2] = 0.5 * sol * (BR[3] - BL[3]) + 0.5 * (ER[2] - EL[2])
    slop[1, 3] = 0.5 * sol * χ * (ER[1] - EL[1])
    slop[7, 3] = 0.5 * χ * (ER[1] - EL[1])
    slop[4, 4] = 0.5 * sol * ν * (BR[1] - BR[1])
    slop[8, 4] = 0.5 * sol^2 * χ * (BR[1] - BR[1])
    slop[3, 5] = -0.5 * sol^2 * (BR[2] - BL[2]) - 0.5 * sol * (ER[3] - EL[3])
    slop[5, 5] = -0.5 * sol * (BR[2] - BL[2]) - 0.5 * (ER[3] - EL[3])
    slop[2, 6] = 0.5 * sol^2 * (BR[3] - BL[3]) - 0.5 * sol * (ER[2] - EL[2])
    slop[6, 6] = -0.5 * sol * (BR[3] - BL[3]) + 0.5 * (ER[2] - EL[2])
    slop[1, 7] = -0.5 * sol * χ * (ER[1] - EL[1])
    slop[7, 7] = 0.5 * χ * (ER[1] - EL[1])
    slop[4, 8] = -0.5 * sol * ν * (BR[1] - BR[1])
    slop[8, 8] = 0.5 * sol^2 * χ * (BR[1] - BR[1])

    limiter = zeros(8, 8)
    limiter[3, 1] = -0.5 * sol^2 * (BL[2] - BLL[2]) + 0.5 * sol * (EL[3] - ELL[3])
    limiter[5, 1] = 0.5 * sol * (BL[2] - BLL[2]) - 0.5 * (EL[3] - ELL[3])
    limiter[2, 2] = 0.5 * sol^2 * (BL[3] - BLL[3]) + 0.5 * sol * (EL[2] - ELL[2])
    limiter[6, 2] = 0.5 * sol * (BL[3] - BLL[3]) + 0.5 * (EL[2] - ELL[2])
    limiter[1, 3] = 0.5 * sol * χ * (EL[1] - ELL[1])
    limiter[7, 3] = 0.5 * χ * (EL[1] - ELL[1])
    limiter[4, 4] = 0.5 * sol * ν * (BL[1] - BL[1])
    limiter[8, 4] = 0.5 * sol^2 * χ * (BL[1] - BL[1])
    limiter[3, 5] = -0.5 * sol^2 * (BRR[2] - BR[2]) - 0.5 * sol * (ERR[3] - ER[3])
    limiter[5, 5] = -0.5 * sol * (BRR[2] - BR[2]) - 0.5 * (ERR[3] - ER[3])
    limiter[2, 6] = 0.5 * sol^2 * (BRR[3] - BR[3]) - 0.5 * sol * (ERR[2] - ER[2])
    limiter[6, 6] = -0.5 * sol * (BRR[3] - BR[3]) + 0.5 * (ERR[2] - ER[2])
    limiter[1, 7] = -0.5 * sol * χ * (ERR[1] - ER[1])
    limiter[7, 7] = 0.5 * χ * (ERR[1] - ER[1])
    limiter[4, 8] = -0.5 * sol * ν * (BRR[1] - BRR[1])
    limiter[8, 8] = 0.5 * sol^2 * χ * (BRR[1] - BRR[1])

    for i = 1:8
        limiter_theta = sum(slop[:, i] .* limiter[:, i]) / (sum(slop[:, i] .^ 2) + 1.e-7)
        slop[:, i] .*=
            max(0.0, min(min((1.0 + limiter_theta) / 2.0, 2.0), 2.0 * limiter_theta))
    end

    for i = 1:8
        femL[i] =
            sum(A1n[i, 1:3] .* (ER .- EL)) +
            sum(A1n[i, 4:6] .* (BR .- BL)) +
            A1n[i, 7] * (ϕR - ϕL) +
            A1n[i, 8] * (ψR - ψL) +
            0.5 * sum(
                fortsign.(1.0, D1) .* (1.0 .- dt ./ (0.5 * (dxL + dxR)) .* abs.(D1)) .*
                slop[i, :],
            )
        femR[i] =
            sum(A1p[i, 1:3] .* (ER .- EL)) +
            sum(A1p[i, 4:6] .* (BR .- BL)) +
            A1p[i, 7] * (ϕR - ϕL) +
            A1p[i, 8] * (ψR - ψL) -
            0.5 * sum(
                fortsign.(1.0, D1) .* (1.0 .- dt ./ (0.5 * (dxL + dxR)) .* abs.(D1)) .*
                slop[i, :],
            )
    end

end
