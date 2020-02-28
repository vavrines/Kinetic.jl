# ============================================================
# Numerical flux functions for finite volume method
# ============================================================


using SpecialFunctions


# ------------------------------------------------------------
# Calculate moments of Maxwellian distribution function
# ------------------------------------------------------------
function get_moments(prim::Array{Float64,1}, inK::Float64)

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


function get_moment_uv(Mu::OffsetArray{Float64,1}, Mxi::OffsetArray{Float64,1}, alpha::Int64, delta::Int64)

    uv = zeros(3)

    uv[1] = Mu[alpha] * Mxi[Int64(delta/2)]
    uv[2] = Mu[alpha+1] * Mxi[Int64(delta/2)]
    uv[3] = 0.5 * (Mu[alpha+2] * Mxi[Int64(delta/2)] + Mu[alpha] * Mxi[Int64((delta+2)/2)])

    return uv

end


function get_moment_uv(Mu::OffsetArray{Float64,1}, Mv::OffsetArray{Float64,1}, Mxi::OffsetArray{Float64,1}, alpha::Int64, beta::Int64, delta::Int64)

    uv = zeros(4)

    uv[1] = Mu[alpha] * Mv[beta] * Mxi[Int64(delta/2)]
    uv[2] = Mu[alpha+1] * Mv[beta] * Mxi[Int64(delta/2)]
    uv[3] = Mu[alpha] * Mv[beta+1] * Mxi[Int64(delta/2)]
    uv[4] = 0.5 * (Mu[alpha+2] * Mv[beta] * Mxi[delta÷2] + Mu[alpha] * Mv[beta+2] * Mxi[delta÷2] + Mu[alpha] * Mv[beta] * Mxi[(delta+2)÷2])

    return uv

end


function get_moment_au(a::Array{Float64,1}, Mu::OffsetArray{Float64,1}, Mxi::OffsetArray{Float64,1}, alpha::Int64)

    au = zeros(3)

    au = a[1] .* get_moment_uv(Mu, Mxi, alpha+0, 0) +
         a[2] .* get_moment_uv(Mu, Mxi, alpha+1, 0) +
         0.5 .* a[3] .* get_moment_uv(Mu, Mxi, alpha+2, 0) +
         0.5 .* a[3] .* get_moment_uv(Mu, Mxi, alpha+0, 2)

    return au

end


function get_moment_au(a::Array{Float64,1}, Mu::OffsetArray{Float64,1}, Mv::OffsetArray{Float64,1}, Mxi::OffsetArray{Float64,1}, alpha::Int64, beta::Int64)

    au = zeros(4)

    au = a[1] .* get_moment_uv(Mu, Mv, Mxi, alpha+0, beta+0, 0) +
         a[2] .* get_moment_uv(Mu, Mv, Mxi, alpha+1, beta+0, 0) +
         a[3] .* get_moment_uv(Mu, Mv, Mxi, alpha+0, beta+1, 0) +
         0.5 * a[4] * get_moment_uv(Mu, Mv, Mxi, alpha+2, beta+0, 0) +
         0.5 * a[4] * get_moment_uv(Mu, Mv, Mxi, alpha+0, beta+2, 0) +
         0.5 * a[4] * get_moment_uv(Mu, Mv, Mxi, alpha+0, beta+0, 2)

    return au

end