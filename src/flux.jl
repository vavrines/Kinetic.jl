# ============================================================
# Numerical flux functions for finite volume method
# ============================================================


export flux_kfvs,
       flux_kcu,
       normal_moments,
       normal_moments_uv,
       normal_moments_au


"""
Kinetic flux vector splitting (KFVS) method

# >@param[in] : particle distribution functions and their slopes at left/right sides of interface
# >@param[in] : particle velocity quadrature points and weights
# >@param[in] : time step

# >@return : flux of particle distribution function and its velocity moments on conservative variables
"""

# ------------------------------------------------------------
# Pure 1D1F flux
# ------------------------------------------------------------
function flux_kfvs( fL::AbstractArray{Float64,1}, sfL::AbstractArray{Float64,1}, 
                    fR::AbstractArray{Float64,1}, sfR::AbstractArray{Float64,1},
                    u::AbstractArray{Float64,1}, ω::AbstractArray{Float64,1}, dt::Float64 )

    #--- upwind reconstruction ---#
    δ = heaviside.(u)

    f = @. fL * δ + fR * (1. - δ)
    sf = @. sfL * δ + sfR * (1. - δ)

    #--- calculate fluxes ---#
    fw = zeros(3); ff = similar(fL)

    fw[1] = dt * sum(ω .* u .* f) - 0.5 * dt^2 * sum(ω .* u.^2 .* sf)
    fw[2] = dt * sum(ω .* u.^2 .* f) - 0.5 * dt^2 * sum(ω .* u.^3 .* sf)
    fw[3] = dt * 0.5 * sum(ω .* u.^3 .* f) - 0.5 * dt^2 * 0.5 * sum(ω .* u.^4 .* sf)

    @. ff = dt * u * f - 0.5 * dt^2 * u^2 * sf

    return fw, ff

end

flux_kfvs(fL::AbstractArray{Float64,1}, fR::AbstractArray{Float64,1}, u::AbstractArray{Float64,1}, ω::AbstractArray{Float64,1}, dt::Float64) = 
flux_kfvs(fL, zeros(axes(fL)), fR, zeros(axes(fR)), u, ω, dt)


# ------------------------------------------------------------
# Reduced 1D2F flux
# ------------------------------------------------------------
function flux_kfvs( hL::AbstractArray{Float64,1}, bL::AbstractArray{Float64,1}, shL::AbstractArray{Float64,1}, sbL::AbstractArray{Float64,1}, 
                    hR::AbstractArray{Float64,1}, bR::AbstractArray{Float64,1}, shR::AbstractArray{Float64,1}, sbR::AbstractArray{Float64,1}, 
                    u::AbstractArray{Float64,1}, ω::AbstractArray{Float64,1}, dt::Float64 )

    #--- upwind reconstruction ---#
    δ = heaviside.(u)

    h = @. hL * δ + hR * (1. - δ)
    b = @. bL * δ + bR * (1. - δ)

    sh = @. shL * δ + shR * (1. - δ)
    sb = @. sbL * δ + sbR * (1. - δ)

    #--- calculate fluxes ---#
    fw = zeros(5); fh = similar(h); fb = similar(b)

    fw[1] = dt * sum(ω .* u .* h) - 0.5 * dt^2 * sum(ω .* u.^2 .* sh)
    fw[2] = dt * sum(ω .* u.^2 .* h) - 0.5 * dt^2 * sum(ω .* u.^3 .* sh)
    fw[3] = dt * 0.5 * (sum(ω .* u.^3 .* h) + sum(ω .* u .* b)) - 
            0.5 * dt^2 * 0.5 * (sum(ω .* u.^4 .* sh) + sum(ω .* u.^2 .* sb))

    @. fh = dt * u * h - 0.5 * dt^2 * u^2 * sh
    @. fb = dt * u * b - 0.5 * dt^2 * u^2 * sb

    return fw, fh, fb

end

flux_kfvs(hL::AbstractArray{Float64,1}, bL::AbstractArray{Float64,1}, hR::AbstractArray{Float64,1}, bR::AbstractArray{Float64,1},
u::AbstractArray{Float64,1}, ω::AbstractArray{Float64,1}, dt::Float64) = 
flux_kfvs(h0L, h1L, h2L, zeros(axes(h0L)), zeros(axes(h1L)), zeros(axes(h2L)), zeros(axes(h3L)), 
h0R, h1R, h2R, zeros(axes(h0R)), zeros(axes(h1R)), zeros(axes(h2R)), zeros(axes(h3R)), u, ω, dt)     


# ------------------------------------------------------------
# Reduced 1D4F flux
# ------------------------------------------------------------
function flux_kfvs( h0L::AbstractArray{Float64,1}, h1L::AbstractArray{Float64,1}, h2L :: AbstractArray{Float64,1}, h3L :: AbstractArray{Float64,1},
                    sh0L::AbstractArray{Float64,1}, sh1L::AbstractArray{Float64,1}, sh2L::AbstractArray{Float64,1}, sh3L::AbstractArray{Float64,1}, 
                    h0R::AbstractArray{Float64,1}, h1R::AbstractArray{Float64,1}, h2R :: AbstractArray{Float64,1}, h3R :: AbstractArray{Float64,1},
                    sh0R::AbstractArray{Float64,1}, sh1R::AbstractArray{Float64,1}, sh2R::AbstractArray{Float64,1}, sh3R::AbstractArray{Float64,1}, 
                    u::AbstractArray{Float64,1}, ω::AbstractArray{Float64,1}, dt::Float64 )

    #--- upwind reconstruction ---#
    δ = heaviside.(u)

    h0 = @. h0L * δ + h0R * (1. - δ)
    h1 = @. h1L * δ + h1R * (1. - δ)
    h2 = @. h2L * δ + h2R * (1. - δ)
    h3 = @. h3L * δ + h3R * (1. - δ)

    sh0 = @. sh0L * δ + sh0R * (1. - δ)
    sh1 = @. sh1L * δ + sh1R * (1. - δ)
    sh2 = @. sh2L * δ + sh2R * (1. - δ)
    sh3 = @. sh3L * δ + sh3R * (1. - δ)

    #--- calculate fluxes ---#
    fw = zeros(5); fh0 = similar(h0L); fh1 = similar(h1L); fh2 = similar(h2L); fh3 = similar(h3L)

    fw[1] = dt * sum(ω .* u .* h0) - 0.5 * dt^2 * sum(ω .* u.^2 .* sh0)
    fw[2] = dt * sum(ω .* u.^2 .* h0) - 0.5 * dt^2 * sum(ω .* u.^3 .* sh0)
    fw[3] = dt * sum(ω .* u .* h1) - 0.5 * dt^2 * sum(ω .* u.^2 .* sh1)
    fw[4] = dt * sum(ω .* u .* h2) - 0.5 * dt^2 * sum(ω .* u.^2 .* sh2)
    fw[5] = dt * 0.5 * (sum(ω .* u.^3 .* h0) + sum(ω .* u .* h3)) - 
            0.5 * dt^2 * 0.5 * (sum(ω .* u.^4 .* sh0) + sum(ω .* u.^2 .* sh3))

    @. fh0 = dt * u * h0 - 0.5 * dt^2 * u^2 * sh0
    @. fh1 = dt * u * h1 - 0.5 * dt^2 * u^2 * sh1
    @. fh2 = dt * u * h2 - 0.5 * dt^2 * u^2 * sh2
    @. fh3 = dt * u * h3 - 0.5 * dt^2 * u^2 * sh3

    return fw, fh0, fh1, fh2, fh3

end

flux_kfvs(h0L::AbstractArray{Float64,1}, h1L::AbstractArray{Float64,1}, h2L :: AbstractArray{Float64,1}, h3L :: AbstractArray{Float64,1},
          h0R::AbstractArray{Float64,1}, h1R::AbstractArray{Float64,1}, h2R :: AbstractArray{Float64,1}, h3R :: AbstractArray{Float64,1},
          u::AbstractArray{Float64,1}, ω::AbstractArray{Float64,1}, dt::Float64) = 
flux_kfvs(h0L, h1L, h2L, h3L, zeros(axes(h0L)), zeros(axes(h1L)), zeros(axes(h2L)), zeros(axes(h3L)), 
          h0R, h1R, h2R, h3R, zeros(axes(h0R)), zeros(axes(h1R)), zeros(axes(h2R)), zeros(axes(h3R)), u, ω, dt)     


# ------------------------------------------------------------
# Pure 2D1F flux
# ------------------------------------------------------------
function flux_kfvs( fL::AbstractArray{Float64,2}, sfL::AbstractArray{Float64,2}, 
                    fR::AbstractArray{Float64,2}, sfR::AbstractArray{Float64,2},
                    u::AbstractArray{Float64,2}, v::AbstractArray{Float64,2}, 
                    ω::AbstractArray{Float64,2}, dt::Float64, len::Float64 )
    
    #--- upwind reconstruction ---#
    δ = heaviside.(u)

    f = @. fL * δ + fR * (1. - δ)
    sf = @. sfL * δ + sfR * (1. - δ)

    #--- calculate fluxes ---#
    fw = zeros(4); ff = similar(fL)

    fw[1] = dt * sum(ω .* u .* f) - 0.5 * dt^2 * sum(ω .* u.^2 .* sf)
    fw[2] = dt * sum(ω .* u.^2 .* f) - 0.5 * dt^2 * sum(ω .* u.^3 .* sf)
    fw[3] = dt * sum(ω .* v .* u .* f) - 0.5 * dt^2 * sum(ω .* v .* u.^2 .* sf)
    fw[4] = dt * 0.5 * sum(ω .* u .* (u.^2 .+ v.^2) .* f) - 0.5 * dt^2 * 0.5 * sum(ω .* u.^2 .* (u.^2 .+ v.^2) .* sf)

    @. ff = dt * u * f - 0.5 * dt^2 * u^2 * sf

    return fw .* len, ff.* len

end


flux_kfvs(fL::AbstractArray{Float64,2}, fR::AbstractArray{Float64,2}, u::AbstractArray{Float64,2}, v::AbstractArray{Float64,2}, 
          ω::AbstractArray{Float64,2}, dt::Float64, len::Float64) = 
flux_kfvs(fL, zeros(axes(fL)), fR, zeros(axes(fR)), u, v, ω, dt, len)


"""
Kinetic central-upwind (KCU) method

# >@param[in] : particle distribution functions and their slopes at left/right sides of interface
# >@param[in] : particle velocity quadrature points and weights
# >@param[in] : time step

# >@return : flux of particle distribution function and its velocity moments on conservative variables
"""

# ------------------------------------------------------------
# Pure 1D1F flux
# ------------------------------------------------------------
function flux_kcu( wL::Array{Float64,1}, fL::AbstractArray{Float64,1}, 
                   wR::Array{Float64,1}, fR::AbstractArray{Float64,1},
                   u::AbstractArray{Float64,1}, ω::AbstractArray{Float64,1}, 
                   inK::Union{Int64,Float64}, γ::Float64, visRef::Float64, visIdx::Float64, pr::Float64, dt::Float64 )

    #--- upwind reconstruction ---#
    δ = heaviside.(u)
    f = @. fL * δ + fR * (1. - δ)

    primL = conserve_prim(wL, gam)
    primR = conserve_prim(wR, gam)

    #--- construct interface distribution ---#
    Mu1, Mxi1, MuL1, MuR1 = normal_moments(primL, inK)
    Muv1 = normal_moments_uv(MuL1, Mxi1, 0, 0)
    Mu2, Mxi2, MuL2, MuR2 = normal_moments(primR, inK)
    Muv2 = normal_moments_uv(MuR2, Mxi2, 0, 0)

    w = zeros(3)
    @. w = primL[1] * Muv1 + primR[1] * Muv2

    prim = conserve_primitive(w, γ)
    tau = collision_time(prim, visRef, visIdx)
    #tau = tau + abs(cellL.prim[1] / cellL.prim[end] - cellR.prim[1] / cellR.prim[end]) / 
    #       (cellL.prim[1] / cellL.prim[end] + cellR.prim[1] / cellR.prim[end]) * dt * 1.

    Mt = zeros(2)
    Mt[2] = tau * (1. - exp(-dt / tau)) # f0
    Mt[1] = dt - Mt[2] # M0

    #--- calculate fluxes ---#
    Mu, Mxi, MuL, MuR = normal_moments(prim, inK)

    # flux from M0
    Muv = normal_moments_uv(Mu, Mxi, 1, 0)
    fw = @. Mt[1] * prim[1] * Muv

    # flux from f0
    g = maxwellian(u, prim)

    fw[1] += Mt[2] * sum(ω .* u .* f)
    fw[2] += Mt[2] * sum(ω .* u.^2 .* f)
    fw[3] += Mt[2] * 0.5 * (sum(ω .* u.^3 .* f))

    ff = @. Mt[1] * u * g + Mt[2] * u * f

    return fw, ff

end


# ------------------------------------------------------------
# Pure 2D1F flux
# ------------------------------------------------------------
function flux_kcu( wL::Array{Float64,1}, fL::AbstractArray{Float64,2}, 
                   wR::Array{Float64,1}, fR::AbstractArray{Float64,2},
                   u::AbstractArray{Float64,2}, v::AbstractArray{Float64,2}, ω::AbstractArray{Float64,2}, 
                   inK::Union{Int64,Float64}, γ::Float64, visRef::Float64, visIdx::Float64, pr::Float64, 
                   dt::Float64, len::Float64 )
    
    #--- prepare ---#
    delta = heaviside.(u)

    #--- reconstruct initial distribution ---#
    δ = heaviside.(u)
    f = @. fL * δ + fR * (1. - δ)

    primL = conserve_prim(wL, gam)
    primR = conserve_prim(wR, gam)

    #--- construct interface distribution ---#
    Mu1, Mv1, Mxi1, MuL1, MuR1 = normal_moments(primL, inK)
    Muv1 = normal_moments_uv(MuL1, Mv1, Mxi1, 0, 0, 0)
    Mu2, Mv2, Mxi2, MuL2, MuR2 = normal_moments(primR, inK)
    Muv2 = normal_moments_uv(MuR2, Mv2, Mxi2, 0, 0, 0)

    w = @. primL[1] * Muv1 + primR[1] * Muv2
    prim = conserve_prim(w, γ)
    tau = collision_time(prim, visRef, visIdx)

    Mt = zeros(2)
    Mt[2] = tau * (1. - exp(-dt / tau)) # f0
    Mt[1] = dt - Mt[2] # M0

    #--- calculate interface flux ---#
    Mu, Mv, Mxi, MuL, MuR = normal_moments(prim, inK)

    # flux from M0
    Muv = normal_moments_uv(Mu, Mv, Mxi, 1, 0, 0)
    fw = @. Mt[1] * prim[1] * Muv

    # flux from f0
    g = maxwellian(u, v, prim)

    fw[1] += Mt[2] * sum(ω .* u .* f)
    fw[2] += Mt[2] * sum(ω .* u.^2 .* f)
    fw[3] += Mt[2] * sum(ω .* v .* u .* f)
    fw[4] += Mt[2] * 0.5 * (sum(ω .* u .* (u.^2 .+ v.^2) .* f))

    ff = @. Mt[1] * u * g + Mt[2] * u * f

    return fw .* len, ff.* len

end
