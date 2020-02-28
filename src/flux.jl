# ============================================================
# Numerical flux functions for finite volume method
# ============================================================


export flux_kfvs,
       flux_kcu,
       normal_moments,
       normal_moments_uv,
       normal_moments_au,
       heaviside


# ------------------------------------------------------------
# Kinetic flux vector splitting of particle distribution function
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


# ------------------------------------------------------------
# central-upwind flux of particle distribution function
# ------------------------------------------------------------
function flux_kcu( wL::Array{Float64,1}, fL::AbstractArray{Float64,1}, 
                   wR::Array{Float64,1}, fR::AbstractArray{Float64,1},
                   u::AbstractArray{Float64,1}, ω::AbstractArray{Float64,1}, 
                   inK::Union{Int,Float64}, γ::Float64, visRef::Float64, visIdx::Float64, pr::Float64, dt::Float64 )

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


function flux_kcu( wL::Array{Float64,1}, fL::AbstractArray{Float64,2}, 
                   wR::Array{Float64,1}, fR::AbstractArray{Float64,2},
                   u::AbstractArray{Float64,2}, v::AbstractArray{Float64,2}, ω::AbstractArray{Float64,2}, 
                   inK::Union{Int,Float64}, γ::Float64, visRef::Float64, visIdx::Float64, pr::Float64, 
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


# ------------------------------------------------------------
# Calculate moments of Maxwellian distribution function
# ------------------------------------------------------------
function normal_moments(prim::Array{Float64,1}, inK::Union{Int, Float64})

    MuL = zeros(7); MuR = zeros(7); Mu = zeros(7); Mxi = zeros(3)
    MuL = OffsetArray(MuL, 0:6); MuR = OffsetArray(MuR, 0:6); Mu = OffsetArray(Mu, 0:6); Mxi = OffsetArray(Mxi, 0:2)

    MuL[0] = 0.5 * erfc(-sqrt(prim[end]) * prim[2])
    MuL[1] = prim[2] * MuL[0] + 0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])
    MuR[0] = 0.5 * erfc(sqrt(prim[end]) * prim[2])
    MuR[1] = prim[2] * MuR[0] - 0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])

    for i=2:6
        MuL[i] = prim[2] * MuL[i-1] + 0.5 * (i-1) * MuL[i-2] / prim[end]
        MuR[i] = prim[2] * MuR[i-1] + 0.5 * (i-1) * MuR[i-2] / prim[end]
    end

    Mu = MuL .+ MuR

    Mxi[0] = 1.0
    Mxi[1] = 0.5 * inK / prim[end]
    Mxi[2] = (inK^2 + 2.0 * inK) / (4.0 * prim[end]^2)

    if length(prim) == 4
        Mv = zeros(7)
        Mv = OffsetArray(Mv, 0:6)
        
        Mv[0] = 1.0
        Mv[1] = prim[3]
        for i=2:6
            Mv[i] = prim[3] * Mv[i-1] + 0.5 * (i-1) * Mv[i-2] / prim[end]
        end

        return Mu, Mv, Mxi, MuL, MuR
    else
        return Mu, Mxi, MuL, MuR
    end

end


function normal_moments_uv(Mu::OffsetArray{Float64,1}, Mxi::OffsetArray{Float64,1}, alpha::Int64, delta::Int64)

    uv = zeros(3)

    uv[1] = Mu[alpha] * Mxi[Int64(delta/2)]
    uv[2] = Mu[alpha+1] * Mxi[Int64(delta/2)]
    uv[3] = 0.5 * (Mu[alpha+2] * Mxi[Int64(delta/2)] + Mu[alpha] * Mxi[Int64((delta+2)/2)])

    return uv

end


function normal_moments_uv(Mu::OffsetArray{Float64,1}, Mv::OffsetArray{Float64,1}, Mxi::OffsetArray{Float64,1}, alpha::Int64, beta::Int64, delta::Int64)

    uv = zeros(4)

    uv[1] = Mu[alpha] * Mv[beta] * Mxi[Int64(delta/2)]
    uv[2] = Mu[alpha+1] * Mv[beta] * Mxi[Int64(delta/2)]
    uv[3] = Mu[alpha] * Mv[beta+1] * Mxi[Int64(delta/2)]
    uv[4] = 0.5 * (Mu[alpha+2] * Mv[beta] * Mxi[delta÷2] + Mu[alpha] * Mv[beta+2] * Mxi[delta÷2] + Mu[alpha] * Mv[beta] * Mxi[(delta+2)÷2])

    return uv

end


function normal_moments_au(a::Array{Float64,1}, Mu::OffsetArray{Float64,1}, Mxi::OffsetArray{Float64,1}, alpha::Int64)

    au = zeros(3)

    au = a[1] .* get_moment_uv(Mu, Mxi, alpha+0, 0) +
         a[2] .* get_moment_uv(Mu, Mxi, alpha+1, 0) +
         0.5 .* a[3] .* get_moment_uv(Mu, Mxi, alpha+2, 0) +
         0.5 .* a[3] .* get_moment_uv(Mu, Mxi, alpha+0, 2)

    return au

end


function normal_moments_au(a::Array{Float64,1}, Mu::OffsetArray{Float64,1}, Mv::OffsetArray{Float64,1}, Mxi::OffsetArray{Float64,1}, alpha::Int64, beta::Int64)

    au = zeros(4)

    au = a[1] .* get_moment_uv(Mu, Mv, Mxi, alpha+0, beta+0, 0) +
         a[2] .* get_moment_uv(Mu, Mv, Mxi, alpha+1, beta+0, 0) +
         a[3] .* get_moment_uv(Mu, Mv, Mxi, alpha+0, beta+1, 0) +
         0.5 * a[4] * get_moment_uv(Mu, Mv, Mxi, alpha+2, beta+0, 0) +
         0.5 * a[4] * get_moment_uv(Mu, Mv, Mxi, alpha+0, beta+2, 0) +
         0.5 * a[4] * get_moment_uv(Mu, Mv, Mxi, alpha+0, beta+0, 2)

    return au

end


heaviside(x::Union{Int, AbstractFloat}) = ifelse(x >= 0, 1., 0.)