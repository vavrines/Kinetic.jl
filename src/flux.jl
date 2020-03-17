# ============================================================
# Numerical flux functions for finite volume method
# ============================================================


export flux_kfvs,
       flux_kcu


"""
Kinetic flux vector splitting (KFVS) method

# >@param[in] : particle distribution functions and their slopes at left/right sides of interface
# >@param[in] : particle velocity quadrature points and weights
# >@param[in] : time step

# >@return : flux of particle distribution function and its velocity moments on conservative variables
"""

# ------------------------------------------------------------
# 1D1F flux
# ------------------------------------------------------------
function flux_kfvs( fL::AbstractArray{Float64,1}, fR::AbstractArray{Float64,1}, 
                    u::AbstractArray{Float64,1}, ω::AbstractArray{Float64,1}, dt::Float64,
                    sfL = zeros(axes(fL))::AbstractArray{Float64,1}, sfR = zeros(axes(fR))::AbstractArray{Float64,1} )

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    f = @. fL * δ + fR * (1. - δ)
    sf = @. sfL * δ + sfR * (1. - δ)

    # --- calculate fluxes ---#
    fw = zeros(3); ff = similar(fL)

    fw[1] = dt * sum(ω .* u .* f) - 0.5 * dt^2 * sum(ω .* u.^2 .* sf)
    fw[2] = dt * sum(ω .* u.^2 .* f) - 0.5 * dt^2 * sum(ω .* u.^3 .* sf)
    fw[3] = dt * 0.5 * sum(ω .* u.^3 .* f) - 0.5 * dt^2 * 0.5 * sum(ω .* u.^4 .* sf)

    @. ff = dt * u * f - 0.5 * dt^2 * u^2 * sf

    return fw, ff

end


# ------------------------------------------------------------
# 1D2F flux
# ------------------------------------------------------------
function flux_kfvs( hL::AbstractArray{Float64,1}, bL::AbstractArray{Float64,1},  
                    hR::AbstractArray{Float64,1}, bR::AbstractArray{Float64,1}, 
                    u::AbstractArray{Float64,1}, ω::AbstractArray{Float64,1}, dt::Float64,
                    shL = zeros(axes(hL))::AbstractArray{Float64,1}, sbL = zeros(axes(bL))::AbstractArray{Float64,1},
                    shR = zeros(axes(hR))::AbstractArray{Float64,1}, sbR = zeros(axes(bR))::AbstractArray{Float64,1} )

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    h = @. hL * δ + hR * (1. - δ)
    b = @. bL * δ + bR * (1. - δ)

    sh = @. shL * δ + shR * (1. - δ)
    sb = @. sbL * δ + sbR * (1. - δ)

    # --- calculate fluxes ---#
    fw = zeros(5); fh = similar(h); fb = similar(b)

    fw[1] = dt * sum(ω .* u .* h) - 0.5 * dt^2 * sum(ω .* u.^2 .* sh)
    fw[2] = dt * sum(ω .* u.^2 .* h) - 0.5 * dt^2 * sum(ω .* u.^3 .* sh)
    fw[3] = dt * 0.5 * (sum(ω .* u.^3 .* h) + sum(ω .* u .* b)) - 
            0.5 * dt^2 * 0.5 * (sum(ω .* u.^4 .* sh) + sum(ω .* u.^2 .* sb))

    @. fh = dt * u * h - 0.5 * dt^2 * u^2 * sh
    @. fb = dt * u * b - 0.5 * dt^2 * u^2 * sb

    return fw, fh, fb

end


# ------------------------------------------------------------
# 1D4F flux
# ------------------------------------------------------------
function flux_kfvs( h0L::AbstractArray{Float64,1}, h1L::AbstractArray{Float64,1}, h2L::AbstractArray{Float64,1}, h3L::AbstractArray{Float64,1}, 
                    h0R::AbstractArray{Float64,1}, h1R::AbstractArray{Float64,1}, h2R::AbstractArray{Float64,1}, h3R::AbstractArray{Float64,1},   
                    u::AbstractArray{Float64,1}, ω::AbstractArray{Float64,1}, dt::Float64,
                    sh0L = zeros(axes(h0L))::AbstractArray{Float64,1}, sh1L = zeros(axes(h1L))::AbstractArray{Float64,1}, 
                    sh2L = zeros(axes(h2L))::AbstractArray{Float64,1}, sh3L = zeros(axes(h3L))::AbstractArray{Float64,1},
                    sh0R = zeros(axes(h0R))::AbstractArray{Float64,1}, sh1R = zeros(axes(h1R))::AbstractArray{Float64,1}, 
                    sh2R = zeros(axes(h2R))::AbstractArray{Float64,1}, sh3R = zeros(axes(h3R))::AbstractArray{Float64,1} )

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    h0 = @. h0L * δ + h0R * (1. - δ)
    h1 = @. h1L * δ + h1R * (1. - δ)
    h2 = @. h2L * δ + h2R * (1. - δ)
    h3 = @. h3L * δ + h3R * (1. - δ)

    sh0 = @. sh0L * δ + sh0R * (1. - δ)
    sh1 = @. sh1L * δ + sh1R * (1. - δ)
    sh2 = @. sh2L * δ + sh2R * (1. - δ)
    sh3 = @. sh3L * δ + sh3R * (1. - δ)

    # --- calculate fluxes ---#
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


# ------------------------------------------------------------
# 2D1F flux
# ------------------------------------------------------------
function flux_kfvs( fL::AbstractArray{Float64,2}, fR::AbstractArray{Float64,2},
                    u::AbstractArray{Float64,2}, v::AbstractArray{Float64,2}, 
                    ω::AbstractArray{Float64,2}, dt::Float64, len::Float64,
                    sfL = zeros(axes(fL))::AbstractArray{Float64,2}, sfR = zeros(axes(fR))::AbstractArray{Float64,2} )
    
    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    f = @. fL * δ + fR * (1. - δ)
    sf = @. sfL * δ + sfR * (1. - δ)

    # --- calculate fluxes ---#
    fw = zeros(4); ff = similar(fL)

    fw[1] = dt * sum(ω .* u .* f) - 0.5 * dt^2 * sum(ω .* u.^2 .* sf)
    fw[2] = dt * sum(ω .* u.^2 .* f) - 0.5 * dt^2 * sum(ω .* u.^3 .* sf)
    fw[3] = dt * sum(ω .* v .* u .* f) - 0.5 * dt^2 * sum(ω .* v .* u.^2 .* sf)
    fw[4] = dt * 0.5 * sum(ω .* u .* (u.^2 .+ v.^2) .* f) - 0.5 * dt^2 * 0.5 * sum(ω .* u.^2 .* (u.^2 .+ v.^2) .* sf)

    @. ff = dt * u * f - 0.5 * dt^2 * u^2 * sf

    return fw .* len, ff .* len

end


"""
Kinetic central-upwind (KCU) method

# >@param[in] : particle distribution functions and their slopes at left/right sides of interface
# >@param[in] : particle velocity quadrature points and weights
# >@param[in] : time step

# >@return : flux of particle distribution function and its velocity moments on conservative variables
"""

# ------------------------------------------------------------
# 1D1F flux
# ------------------------------------------------------------
function flux_kcu( wL::Array{Float64,1}, fL::AbstractArray{Float64,1}, 
                   wR::Array{Float64,1}, fR::AbstractArray{Float64,1},
                   u::AbstractArray{Float64,1}, ω::AbstractArray{Float64,1}, 
                   inK::Union{Int64,Float64}, γ::Float64, visRef::Float64, visIdx::Float64, pr::Float64, dt::Float64 )

    # --- upwind reconstruction ---#
    δ = heaviside.(u)
    f = @. fL * δ + fR * (1. - δ)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    # --- construct interface distribution ---#
    Mu1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Muv1 = moments_conserve(MuL1, Mxi1, 0, 0)
    Mu2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)
    Muv2 = moments_conserve(MuR2, Mxi2, 0, 0)

    w = zeros(3)
    @. w = primL[1] * Muv1 + primR[1] * Muv2

    prim = conserve_prim(w, γ)
    tau = vhs_collision_time(prim, visRef, visIdx)
    # tau = tau + abs(cellL.prim[1] / cellL.prim[end] - cellR.prim[1] / cellR.prim[end]) / 
    #       (cellL.prim[1] / cellL.prim[end] + cellR.prim[1] / cellR.prim[end]) * dt * 1.

    Mt = zeros(2)
    Mt[2] = tau * (1. - exp(-dt / tau)) # f0
    Mt[1] = dt - Mt[2] # M0

    # --- calculate fluxes ---#
    Mu, Mxi, MuL, MuR = gauss_moments(prim, inK)

    # flux from M0
    Muv = moments_conserve(Mu, Mxi, 1, 0)
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
# 1D2F flux
# ------------------------------------------------------------
function flux_kcu( wL::Array{Float64,1}, hL::AbstractArray{Float64,1}, bL::AbstractArray{Float64,1},
                   wR::Array{Float64,1}, hR::AbstractArray{Float64,1}, bR::AbstractArray{Float64,1},
                   u::AbstractArray{Float64,1}, ω::AbstractArray{Float64,1}, 
                   inK::Union{Int64,Float64}, γ::Float64, visRef::Float64, visIdx::Float64, pr::Float64, dt::Float64 )

    # --- upwind reconstruction ---#
    δ = heaviside.(u)
    h = @. hL * δ + hR * (1. - δ)
    b = @. bL * δ + bR * (1. - δ)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    # --- construct interface distribution ---#
    Mu1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Muv1 = moments_conserve(MuL1, Mxi1, 0, 0)
    Mu2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)
    Muv2 = moments_conserve(MuR2, Mxi2, 0, 0)

    w = zeros(3)
    @. w = primL[1] * Muv1 + primR[1] * Muv2

    prim = conserve_prim(w, γ)
    tau = vhs_collision_time(prim, visRef, visIdx)
    # tau = tau + abs(cellL.prim[1] / cellL.prim[end] - cellR.prim[1] / cellR.prim[end]) / 
    #       (cellL.prim[1] / cellL.prim[end] + cellR.prim[1] / cellR.prim[end]) * dt * 1.

    Mt = zeros(2)
    Mt[2] = tau * (1. - exp(-dt / tau)) # f0
    Mt[1] = dt - Mt[2] # M0

    # --- calculate fluxes ---#
    Mu, Mxi, MuL, MuR = gauss_moments(prim, inK)

    # flux from M0
    Muv = moments_conserve(Mu, Mxi, 1, 0)
    fw = @. Mt[1] * prim[1] * Muv

    # flux from f0
    Mh = maxwellian(u, prim)
    Mb = Mh .* inK ./ (2. * prim[end])

    fw[1] += Mt[2] * sum(ω .* u .* h)
    fw[2] += Mt[2] * sum(ω .* u.^2 .* h)
    fw[3] += Mt[2] * 0.5 * (sum(ω .* u.^3 .* h) + sum(ω .* u .* b))

    fh = @. Mt[1] * u * Mh + Mt[2] * u * h
    fb = @. Mt[1] * u * Mb + Mt[2] * u * b

    return fw, fh, fb

end


# ------------------------------------------------------------
# 2D1F flux
# ------------------------------------------------------------
function flux_kcu( wL::Array{Float64,1}, fL::AbstractArray{Float64,2}, 
                   wR::Array{Float64,1}, fR::AbstractArray{Float64,2},
                   u::AbstractArray{Float64,2}, v::AbstractArray{Float64,2}, ω::AbstractArray{Float64,2}, 
                   inK::Union{Int64,Float64}, γ::Float64, visRef::Float64, visIdx::Float64, pr::Float64, 
                   dt::Float64, len::Float64 )
    
    # --- prepare ---#
    delta = heaviside.(u)

    # --- reconstruct initial distribution ---#
    δ = heaviside.(u)
    f = @. fL * δ + fR * (1. - δ)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    # --- construct interface distribution ---#
    Mu1, Mv1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Muv1 = moments_conserve(MuL1, Mv1, Mxi1, 0, 0, 0)
    Mu2, Mv2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)
    Muv2 = moments_conserve(MuR2, Mv2, Mxi2, 0, 0, 0)

    w = @. primL[1] * Muv1 + primR[1] * Muv2
    prim = conserve_prim(w, γ)
    tau = vhs_collision_time(prim, visRef, visIdx)

    Mt = zeros(2)
    Mt[2] = tau * (1. - exp(-dt / tau)) # f0
    Mt[1] = dt - Mt[2] # M0

    # --- calculate interface flux ---#
    Mu, Mv, Mxi, MuL, MuR = gauss_moments(prim, inK)

    # flux from M0
    Muv = moments_conserve(Mu, Mv, Mxi, 1, 0, 0)
    fw = @. Mt[1] * prim[1] * Muv

    # flux from f0
    g = maxwellian(u, v, prim)

    fw[1] += Mt[2] * sum(ω .* u .* f)
    fw[2] += Mt[2] * sum(ω .* u.^2 .* f)
    fw[3] += Mt[2] * sum(ω .* v .* u .* f)
    fw[4] += Mt[2] * 0.5 * (sum(ω .* u .* (u.^2 .+ v.^2) .* f))

    ff = @. Mt[1] * u * g + Mt[2] * u * f

    return fw .* len, ff .* len

end


"""
Kinetic central-upwind (KCU) method for multi-component gas

# >@param[in] : particle distribution functions and their slopes at left/right sides of interface
# >@param[in] : particle velocity quadrature points and weights
# >@param[in] : time step

# >@return : flux of particle distribution function and its velocity moments on conservative variables
"""

# ------------------------------------------------------------
# 1D1F flux with AAP model
# ------------------------------------------------------------
function flux_kcu( wL::Array{Float64,2}, fL::AbstractArray{Float64,2},  
                   wR::Array{Float64,2}, fR::AbstractArray{Float64,2}, 
                   u::AbstractArray{Float64,2}, ω::AbstractArray{Float64,2}, inK::Union{Int64,Float64}, γ::Float64, 
                   mi::Float64, ni::Float64, me::Float64, ne::Float64, kn::Float64, dt::Float64 )
    
    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    f = @. fL * δ + fR * (1. - δ)

    primL = zeros(axes(wL)); primR = similar(primL)
    for j in 1:2
        primL[:,j] .= conserve_prim(wL[:,j], γ)
        primR{:,j} .= conserve_prim(wR[:,j], γ)
    end

    # --- construct interface distribution ---#
    Mu1 = OffsetArray{Float64}(undef, 0:6, 1:2); Mxi1 = similar(Mu1); MuL1 = similar(Mu1); MuR1 = similar(Mu1)
    Mu2 = similar(Mu1); Mxi2 = similar(Mu1); MuL2 = similar(Mu1); MuR2 = similar(Mu1)
    Muv1 = similar(wL); Muv2 = similar(wL)
    for j in 1:2
        Mu1[:,j], Mxi1[:,j], MuL1[:,j], MuR1[:,j] = gauss_moments(primL[:,j], inK)
        Muv1[:,j] = moments_conserve(MuL1[:,j], Mxi1[:,j], 0, 0)
        Mu2[:,j], Mxi2[:,j], MuL2[:,j], MuR2[:,j] = gauss_moments(primR[:,j], inK)
        Muv2[:,j] = moments_conserve(MuR2[:,j], Mxi2[:,j], 0, 0)
    end

    w = zeros(axes(wL)); prim = zeros(axes(wL))
    for j in 1:2
        @. w[:,j] = primL[1,j] * Muv1[:,j] + primR[1,j] * Muv2[:,j]
        prim[:,j] .= conserve_prim(w[:,j], γ)
    end

    tau = aap_hs_collision_time(prim, mi, ni, me, ne, kn)
    # tau .+= abs(cellL.prim[1,:] / cellL.prim[end,:] - cellR.prim[1,:] / cellR.prim[end,:]) / 
    #         (cellL.prim[1,:] / cellL.prim[end,:] + cellR.prim[1,:] / cellR.prim[end,:]) * dt * 2.
    prim = aap_hs_prim(prim, tau, mi, ni, me, ne, kn)

    Mt = zeros(2, 2)
    @. Mt[2,:] = tau * (1. - exp(-dt / tau)) # f0
    @. Mt[1,:] = dt - Mt[2,:] # M0

    # --- calculate fluxes ---#
    Mu = similar(Mu1); Mxi = similar(Mu1); MuL = similar(Mu1); MuR = similar(Mu1)
    Muv = similar(wL)
    for j in axes(Mu1, 2)
        Mu[:,j], Mxi[:,j], MuL[:,j], MuR[:,j] = gauss_moments(prim[:,j], inK)
        Muv[:,j] .= moments_conserve(MuL[:,j], Mxi[:,j], 1, 0)
    end

    # flux from M0
    fw = similar(wL)
    for j in 1:2
        @. fw[:,j] = Mt[1,j] * prim[1,j] * Muv[:,j]
    end

    # flux from f0
    g = similar(f)
    for j in 1:2
        g[:,j] .= maxwellian(u[:,j], prim[:,j])
    end

    ff = similar(f); 
    for j = 1:2
        fw[1,j] += Mt[2,j] * sum(ω[:,j] .* u[:,j] .* f[:,j])
        fw[2,j] += Mt[2,j] * sum(ω[:,j] .* u[:,j].^2 .* f[:,j])
        fw[3,j] += Mt[2,j] * 0.5 * sum(ω[:,j] .* u[:,j].^3 .* f[:,j])

        @. ff[:,j] = Mt[1,j] * u[:,j] * g[:,j] + Mt[2,j] * u[:,j] * f[:,j]
    end

    return fw, ff

end


# ------------------------------------------------------------
# 1D4F flux with AAP model
# ------------------------------------------------------------
function flux_kcu( wL::Array{Float64,2}, hL::AbstractArray{Float64,2}, bL::AbstractArray{Float64,2},  
                   wR::Array{Float64,2}, hR::AbstractArray{Float64,2}, bR::AbstractArray{Float64,2},  
                   u::AbstractArray{Float64,2}, ω::AbstractArray{Float64,2}, inK::Union{Int64,Float64}, γ::Float64, 
                   mi::Float64, ni::Float64, me::Float64, ne::Float64, kn::Float64, dt::Float64 )

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    h = @. hL * δ + hR * (1. - δ)
    b = @. bL * δ + bR * (1. - δ)

    primL = zeros(axes(wL)); primR = similar(primL)
    for j in 1:2
        primL[:,j] .= conserve_prim(wL[:,j], γ)
        primR{:,j} .= conserve_prim(wR[:,j], γ)
    end

    # --- construct interface distribution ---#
    Mu1 = OffsetArray{Float64}(undef, 0:6, 1:2); Mxi1 = similar(Mu1); MuL1 = similar(Mu1); MuR1 = similar(Mu1)
    Mu2 = similar(Mu1); Mxi2 = similar(Mu1); MuL2 = similar(Mu1); MuR2 = similar(Mu1)
    Muv1 = similar(wL); Muv2 = similar(wL)
    for j in 1:2
        Mu1[:,j], Mxi1[:,j], MuL1[:,j], MuR1[:,j] = gauss_moments(primL[:,j], inK)
        Muv1[:,j] = moments_conserve(MuL1[:,j], Mxi1[:,j], 0, 0)
        Mu2[:,j], Mxi2[:,j], MuL2[:,j], MuR2[:,j] = gauss_moments(primR[:,j], inK)
        Muv2[:,j] = moments_conserve(MuR2[:,j], Mxi2[:,j], 0, 0)
    end

    w = zeros(axes(wL)); prim = zeros(axes(wL))
    for j in 1:2
        @. w[:,j] = primL[1,j] * Muv1[:,j] + primR[1,j] * Muv2[:,j]
        prim[:,j] .= conserve_prim(w[:,j], γ)
    end

    tau = aap_hs_collision_time(prim, mi, ni, me, ne, kn)
    # tau .+= abs(cellL.prim[1,:] / cellL.prim[end,:] - cellR.prim[1,:] / cellR.prim[end,:]) / 
    #         (cellL.prim[1,:] / cellL.prim[end,:] + cellR.prim[1,:] / cellR.prim[end,:]) * dt * 2.
    prim = aap_hs_prim(prim, tau, mi, ni, me, ne, kn)

    Mt = zeros(2, 2)
    @. Mt[2,:] = tau * (1. - exp(-dt / tau)) # f0
    @. Mt[1,:] = dt - Mt[2,:] # M0

    # --- calculate fluxes ---#
    Mu = similar(Mu1); Mxi = similar(Mu1); MuL = similar(Mu1); MuR = similar(Mu1)
    Muv = similar(wL)
    for j in axes(Mu1, 2)
        Mu[:,j], Mxi[:,j], MuL[:,j], MuR[:,j] = gauss_moments(prim[:,j], inK)
        Muv[:,j] .= moments_conserve(MuL[:,j], Mxi[:,j], 1, 0)
    end

    # flux from M0
    fw = similar(wL)
    for j in 1:2
        @. fw[:,j] = Mt[1,j] * prim[1,j] * Muv[:,j]
    end

    # flux from f0
    g0 = similar(h); g1 = similar(b)
    for j in 1:2
        g0[:,j] .= maxwellian(u[:,j], prim[:,j])
        g1[:,j] .= g0[:,j] .* inK ./ (2. * prim[end,j])
    end

    fh = similar(h); fb = similar(b)
    for j = 1:2
        fw[1,j] += Mt[2,j] * sum(ω[:,j] .* u[:,j] .* h[:,j])
        fw[2,j] += Mt[2,j] * sum(ω[:,j] .* u[:,j].^2 .* h[:,j])
        fw[3,j] += Mt[2,j] * 0.5 * (sum(ω[:,j] .* u[:,j].^3 .* h[:,j]) + sum(ω[:,j] .* u[:,j] .* b[:,j]))

        @. fh[:,j] = Mt[1,j] * u[:,j] * g0[:,j] + Mt[2,j] * u[:,j] * h[:,j]
        @. fb[:,j] = Mt[1,j] * u[:,j] * g1[:,j] + Mt[2,j] * u[:,j] * b[:,j]
    end

    return fw, fh, fb

end


# ------------------------------------------------------------
# 1D4F flux with AAP model
# ------------------------------------------------------------
function flux_kcu( wL::Array{Float64,2}, h0L::AbstractArray{Float64,2}, h1L::AbstractArray{Float64,2}, h2L::AbstractArray{Float64,2}, h3L::AbstractArray{Float64,2}, 
                   wR::Array{Float64,2}, h0R::AbstractArray{Float64,2}, h1R::AbstractArray{Float64,2}, h2R::AbstractArray{Float64,2}, h3R::AbstractArray{Float64,2},   
                   u::AbstractArray{Float64,2}, ω::AbstractArray{Float64,2}, inK::Union{Int64,Float64}, γ::Float64, 
                   mi::Float64, ni::Float64, me::Float64, ne::Float64, kn::Float64, dt::Float64 )

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    h0 = @. h0L * δ + h0R * (1. - δ)
    h1 = @. h1L * δ + h1R * (1. - δ)
    h2 = @. h2L * δ + h2R * (1. - δ)
    h3 = @. h3L * δ + h3R * (1. - δ)

    primL = zeros(axes(wL)); primR = similar(primL)
    for j in 1:2
        primL[:,j] .= conserve_prim(wL[:,j], γ)
        primR{:,j} .= conserve_prim(wR[:,j], γ)
    end

    # --- construct interface distribution ---#
    Mu1 = OffsetArray{Float64}(undef, 0:6, 1:2); Mv1 = similar(Mu1); Mw1 = similar(Mu1); MuL1 = similar(Mu1); MuR1 = similar(Mu1)
    Mu2 = similar(Mu1); Mv2 = similar(Mu1); Mw2 = similar(Mu1); MuL2 = similar(Mu1); MuR2 = similar(Mu1)
    Muv1 = similar(wL); Muv2 = similar(wL)
    for j in 1:2
        Mu1[:,j], Mv1[:,j], Mw1[:,j], MuL1[:,j], MuR1[:,j] = gauss_moments(primL[:,j], inK)
        Muv1[:,j] = moments_conserve(MuL1[:,j], Mv1[:,j], Mw1[:,j], 0, 0, 0)
        Mu2[:,j], Mv2[:,j], Mw2[:,j], MuL2[:,j], MuR2[:,j] = gauss_moments(primR[:,j], inK)
        Muv2[:,j] = moments_conserve(MuR2[:,j], Mv2[:,j], Mw2[:,j], 0, 0, 0)
    end

    w = zeros(axes(wL)); prim = zeros(axes(wL))
    for j in 1:2
        @. w[:,j] = primL[1,j] * Muv1[:,j] + primR[1,j] * Muv2[:,j]
        prim[:,j] .= conserve_prim(w[:,j], γ)
    end

    tau = aap_hs_collision_time(prim, mi, ni, me, ne, kn)
    # tau .+= abs(cellL.prim[1,:] / cellL.prim[end,:] - cellR.prim[1,:] / cellR.prim[end,:]) / 
    #        (cellL.prim[1,:] / cellL.prim[end,:] + cellR.prim[1,:] / cellR.prim[end,:]) * dt * 2.
    prim = aap_hs_prim(prim, tau, mi, ni, me, ne, kn)

    Mt = zeros(2, 2)
    @. Mt[2,:] = tau * (1. - exp(-dt / tau)) # f0
    @. Mt[1,:] = dt - Mt[2,:] # M0

    # --- calculate fluxes ---#
    Mu = similar(Mu1); Mv = similar(Mu1); Mw = similar(Mu1); MuL = similar(Mu1); MuR = similar(Mu1)
    Muv = similar(wL)
    for j in axes(Mu1, 2)
        Mu[:,j], Mv[:,j], Mw[:,j], MuL[:,j], MuR[:,j] = gauss_moments(prim[:,j], inK)
        Muv[:,j] .= moments_conserve(MuL[:,j], Mv[:,j], Mw[:,j], 1, 0, 0)
    end

    # flux from M0
    fw = similar(wL)
    for j in 1:2
        @. fw[:,j] = Mt[1,j] * prim[1,j] * Muv[:,j]
    end

    # flux from f0
    g0 = similar(h0); g1 = similar(h0); g2 = similar(h0); g3 = similar(h0)
    for j in 1:2
        g0[:,j] .= maxwellian(u[:,j], prim[:,j])
        g1[:,j] .= Mv[1,j] * g0[:,j]
        g2[:,j] .= Mw[1,j] * g0[:,j]
        g3[:,j] .= (Mv[2,j] + Mw[2,j]) * g0[:,j]
    end

    fh0 = similar(h0); fh1 = similar(h0); fh2 = similar(h0); fh3 = similar(h0)
    for j = 1:2
        fw[1,j] += Mt[2,j] * sum(ω[:,j] .* u[:,j] .* h0[:,j])
        fw[2,j] += Mt[2,j] * sum(ω[:,j] .* u[:,j].^2 .* h0[:,j])
        fw[3,j] += Mt[2,j] * sum(ω[:,j] .* u[:,j] .* h1[:,j])
        fw[4,j] += Mt[2,j] * sum(ω[:,j] .* u[:,j] .* h2[:,j])
        fw[5,j] += Mt[2,j] * 0.5 * (sum(ω[:,j] .* u[:,j].^3 .* h0[:,j]) + sum(ω[:,j] .* u[:,j] .* h3[:,j]))

        @. fh0[:,j] = Mt[1,j] * u[:,j] * g0[:,j] + Mt[2,j] * u[:,j] * h0[:,j]
        @. fh1[:,j] = Mt[1,j] * u[:,j] * g1[:,j] + Mt[2,j] * u[:,j] * h1[:,j]
        @. fh2[:,j] = Mt[1,j] * u[:,j] * g2[:,j] + Mt[2,j] * u[:,j] * h2[:,j]
        @. fh3[:,j] = Mt[1,j] * u[:,j] * g3[:,j] + Mt[2,j] * u[:,j] * h3[:,j]
    end

    return fw, fh0, fh1, fh2, fh3

end