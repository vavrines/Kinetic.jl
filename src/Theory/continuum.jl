# ============================================================
# Continuum Theory
# ============================================================

"""
Transform primitive -> conservative variables

"""
function prim_conserve(prim::AbstractArray{<:Real,1}, γ::Real)

    if eltype(prim) <: Int
        W = zeros(axes(prim))
    else
        W = zeros(eltype(prim), axes(prim))
    end

    if length(prim) == 3 # 1D
        W[1] = prim[1]
        W[2] = prim[1] * prim[2]
        W[3] = 0.5 * prim[1] / prim[3] / (γ - 1.0) + 0.5 * prim[1] * prim[2]^2
    elseif length(prim) == 4 # 2D
        W[1] = prim[1]
        W[2] = prim[1] * prim[2]
        W[3] = prim[1] * prim[3]
        W[4] = 0.5 * prim[1] / prim[4] / (γ - 1.0) + 0.5 * prim[1] * (prim[2]^2 + prim[3]^2)
    elseif length(prim) == 5 # 3D
        W[1] = prim[1]
        W[2] = prim[1] * prim[2]
        W[3] = prim[1] * prim[3]
        W[4] = prim[1] * prim[4]
        W[5] =
            0.5 * prim[1] / prim[5] / (γ - 1.0) +
            0.5 * prim[1] * (prim[2]^2 + prim[3]^2 + prim[4]^2)
    else
        println("prim -> w : dimension error")
    end

    return W

end

prim_conserve(ρ::Real, U::Real, λ::Real, γ::Real) = prim_conserve([ρ, U, λ], γ)

prim_conserve(ρ::Real, U::Real, V::Real, λ::Real, γ::Real) = prim_conserve([ρ, U, V, λ], γ)


function mixture_prim_conserve(prim::AbstractArray{<:Real,2}, γ::Real)
    if eltype(prim) <: Int
        w = zeros(axes(prim))
    else
        w = zeros(eltype(prim), axes(prim))
    end

    for j in axes(w, 2)
        w[:, j] .= prim_conserve(prim[:, j], γ)
    end

    return w
end


"""
Transform conservative -> primitive variables

* scalar: pseudo primitive vector for scalar conservation laws
* vector: primitive vector for Euler, Navier-Stokes and extended equations

"""
conserve_prim(u::Real) = [u, 0.5 * u, 1.0]

conserve_prim(u::Real, a::Real) = [u, a, 1.0]

function conserve_prim(W::AbstractArray{<:Real,1}, γ::Real)

    if eltype(W) <: Int
        prim = zeros(axes(W))
    else
        prim = zeros(eltype(W), axes(W))
    end

    if length(W) == 3 # 1D
        prim[1] = W[1]
        prim[2] = W[2] / W[1]
        prim[3] = 0.5 * W[1] / (γ - 1.0) / (W[3] - 0.5 * W[2]^2 / W[1])
    elseif length(W) == 4 # 2D
        prim[1] = W[1]
        prim[2] = W[2] / W[1]
        prim[3] = W[3] / W[1]
        prim[4] = 0.5 * W[1] / (γ - 1.0) / (W[4] - 0.5 * (W[2]^2 + W[3]^2) / W[1])
    elseif length(W) == 5 # 3D
        prim[1] = W[1]
        prim[2] = W[2] / W[1]
        prim[3] = W[3] / W[1]
        prim[4] = W[4] / W[1]
        prim[5] = 0.5 * W[1] / (γ - 1.0) / (W[5] - 0.5 * (W[2]^2 + W[3]^2 + W[4]^2) / W[1])
    else
        throw("w -> prim : dimension dismatch")
    end

    return prim

end

conserve_prim(ρ::Real, M::Real, E::Real, gamma::Real) = conserve_prim([ρ, M, E], gamma)

conserve_prim(ρ::Real, MX::Real, MY::Real, E::Real, gamma::Real) =
    conserve_prim([ρ, MX, MY, E], gamma)


function mixture_conserve_prim(W::AbstractArray{<:Real,2}, γ::Real)
    if eltype(W) <: Int
        prim = zeros(axes(W))
    else
        prim = zeros(eltype(W), axes(W))
    end

    for j in axes(prim, 2)
        prim[:, j] .= conserve_prim(W[:, j], γ)
    end

    return prim
end


"""
Calculate electromagnetic coeffcients in hyperbolic Maxwell's equations

"""
function em_coefficients(
    prim::AbstractArray{<:Real,2},
    E::AbstractArray{<:Real,1},
    B::AbstractArray{<:Real,1},
    mr::Real,
    lD::Real,
    rL::Real,
    dt::Real,
)

    A = zeros(9, 9)
    A[1, 1] = -1.0 / (2.0 * rL)
    A[2, 2] = -1.0 / (2.0 * rL)
    A[3, 3] = -1.0 / (2.0 * rL)
    A[4, 1] = mr / (2.0 * rL)
    A[5, 2] = mr / (2.0 * rL)
    A[6, 3] = mr / (2.0 * rL)
    A[7, 1] = 1.0 / (dt)
    A[8, 2] = 1.0 / (dt)
    A[9, 3] = 1.0 / (dt)

    A[1, 4] = 1.0 / (dt)
    A[1, 5] = -B[3] / (2.0 * rL)
    A[1, 6] = B[2] / (2.0 * rL)
    A[2, 4] = B[3] / (2.0 * rL)
    A[2, 5] = 1.0 / (dt)
    A[2, 6] = -B[1] / (2.0 * rL)
    A[3, 4] = -B[2] / (2.0 * rL)
    A[3, 5] = B[1] / (2.0 * rL)
    A[3, 6] = 1.0 / (dt)

    A[4, 7] = 1.0 / (dt)
    A[4, 8] = mr * B[3] / (2.0 * rL)
    A[4, 9] = -mr * B[2] / (2.0 * rL)
    A[5, 7] = -mr * B[3] / (2.0 * rL)
    A[5, 8] = 1.0 / (dt)
    A[5, 9] = mr * B[1] / (2.0 * rL)
    A[6, 7] = mr * B[2] / (2.0 * rL)
    A[6, 8] = -mr * B[1] / (2.0 * rL)
    A[6, 9] = 1.0 / (dt)

    A[7, 4] = prim[1, 1] / (2.0 * rL * lD^2)
    A[8, 5] = prim[1, 1] / (2.0 * rL * lD^2)
    A[9, 6] = prim[1, 1] / (2.0 * rL * lD^2)
    A[7, 7] = -(prim[1, 2] * mr) / (2.0 * rL * lD^2)
    A[8, 8] = -(prim[1, 2] * mr) / (2.0 * rL * lD^2)
    A[9, 9] = -(prim[1, 2] * mr) / (2.0 * rL * lD^2)

    b = zeros(9)
    b[1] =
        prim[2, 1] / (dt) + E[1] / (2.0 * rL) - B[2] * prim[4, 1] / (2.0 * rL) +
        B[3] * prim[3, 1] / (2.0 * rL)
    b[2] =
        prim[3, 1] / (dt) + E[2] / (2.0 * rL) - B[3] * prim[2, 1] / (2.0 * rL) +
        B[1] * prim[4, 1] / (2.0 * rL)
    b[3] =
        prim[4, 1] / (dt) + E[3] / (2.0 * rL) - B[1] * prim[3, 1] / (2.0 * rL) +
        B[2] * prim[2, 1] / (2.0 * rL)
    b[4] =
        prim[2, 2] / (dt) - mr * E[1] / (2.0 * rL) + mr * B[2] * prim[4, 2] / (2.0 * rL) -
        mr * B[3] * prim[3, 2] / (2.0 * rL)
    b[5] =
        prim[3, 2] / (dt) - mr * E[2] / (2.0 * rL) + mr * B[3] * prim[2, 2] / (2.0 * rL) -
        mr * B[1] * prim[4, 2] / (2.0 * rL)
    b[6] =
        prim[4, 2] / (dt) - mr * E[3] / (2.0 * rL) + mr * B[1] * prim[3, 2] / (2.0 * rL) -
        mr * B[2] * prim[2, 2] / (2.0 * rL)
    b[7] =
        E[1] / (dt) - prim[1, 1] * prim[2, 1] / (2.0 * rL * lD^2) +
        prim[1, 2] * prim[2, 2] * mr / (2.0 * rL * lD^2)
    b[8] =
        E[2] / (dt) - prim[1, 1] * prim[3, 1] / (2.0 * rL * lD^2) +
        prim[1, 2] * prim[3, 2] * mr / (2.0 * rL * lD^2)
    b[9] =
        E[3] / (dt) - prim[1, 1] * prim[4, 1] / (2.0 * rL * lD^2) +
        prim[1, 2] * prim[4, 2] * mr / (2.0 * rL * lD^2)

    return A, b

end


"""
Theoretical flux of linear advection equation

"""
advection_flux(u::Real, a::Real) = a * u


"""
Theoretical flux of Burgers' equation

"""
burgers_flux(u::Real) = 0.5 * u^2


"""
Theoretical fluxes of Euler Equations

* @return: flux tuple

"""
function euler_flux(w::AbstractArray{<:Real,1}, γ::Real; frame = :cartesian::Symbol)

    prim = conserve_prim(w, γ)
    p = 0.5 * prim[1] / prim[end]

    if length(w) == 3
        F = zeros(axes(w))
        F[1] = w[2]
        F[2] = w[2]^2 / w[1] + p
        F[3] = (w[3] + p) * w[2] / w[1]

        return (F,)
    elseif length(w) == 4
        F = zeros(axes(w))
        F[1] = w[2]
        F[2] = w[2]^2 / w[1] + p
        F[3] = w[2] * w[3] / w[1]
        F[4] = (w[end] + p) * w[2] / w[1]

        G = zeros(axes(w))
        G[1] = w[3]
        G[2] = w[3] * w[2] / w[1]
        G[3] = w[3]^2 / w[1] + p
        G[4] = (w[end] + p) * w[3] / w[1]

        return F, G
    elseif length(w) == 5
        F = zeros(axes(w))
        F[1] = w[2]
        F[2] = w[2]^2 / w[1] + p
        F[3] = w[2] * w[3] / w[1]
        F[4] = w[2] * w[4] / w[1]
        F[5] = (w[end] + p) * w[2] / w[1]

        G = zeros(axes(w))
        G[1] = w[3]
        G[2] = w[3] * w[2] / w[1]
        G[3] = w[3]^2 / w[1] + p
        G[4] = w[3] * w[4] / w[1]
        G[5] = (w[end] + p) * w[3] / w[1]

        H = zeros(axes(w))
        H[1] = w[4]
        H[2] = w[4] * w[2] / w[1]
        H[3] = w[4] * w[3] / w[1]
        H[4] = w[4]^2 / w[1] + p
        H[5] = (w[end] + p) * w[4] / w[1]

        return F, G, H
    end

end


"""
Flux Jacobian of Euler Equations

"""
function euler_jacobi(w::AbstractArray{<:Real,1}, γ::Real)

    A = zeros(eltype(w), 3, 3)

    A[1, 1] = 0.0
    A[1, 2] = 1.0
    A[1, 3] = 0.0

    A[2, 1] = 0.5 * (γ - 3.0) * (w[2] / w[1])^2
    A[2, 2] = (3.0 - γ) * w[2] / w[1]
    A[2, 3] = γ - 1.0

    A[3, 1] = -γ * w[3] * w[2] / w[1]^2 + (γ - 1.0) * (w[2] / w[1])^3
    A[3, 2] = γ * w[3] / w[1] - 1.5 * (γ - 1.0) * (w[2] / w[1])^2
    A[3, 3] = γ * w[2] / w[1]

    return A

end
