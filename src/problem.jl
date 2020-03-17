# ============================================================
# Initial and Boundary Conditions of Specific Problems
# ============================================================


export ib_rh,
       ib_briowu


# ------------------------------------------------------------
# Initialize Rankine-Hugoniot relation
# ------------------------------------------------------------
function ib_rh(MaL::Real, gam::Real, u::AbstractArray{Float64,1})

    #--- calculate Rankine-Hugoniot relation ---#
    primL = [1.0, MaL * sqrt(gam / 2.0), 1.0]

    MaR = sqrt((MaL^2 * (gam - 1.0) + 2.0) / (2.0 * gam * MaL^2 - (gam - 1.0)))
	ratioT = (1.0 + (gam - 1.0) / 2.0 * MaL^2) * (2.0 * gam / (gam - 1.0) * MaL^2 - 1.0) / 
		     (MaL^2 * (2.0 * gam / (gam - 1.0) + (gam - 1.0) / 2.0))

    primR = [ primL[1] * (gam + 1.0) * MaL^2 / ((gam - 1.0) * MaL^2 + 2.0),
              MaR * sqrt(gam / 2.0) * sqrt(ratioT),
              primL[3] / ratioT ]

    wL = prim_conserve(primL, gam)
    wR = prim_conserve(primR, gam)

    hL = maxwellian(u, primL)
    hR = maxwellian(u, primR)

    bcL = zeros(3)
    bcR = zeros(3)

    return wL, primL, hL, bcL, wR, primR, hR, bcR

end

function ib_rh(MaL::Real, gam::Real, u::AbstractArray{Float64,1}, K::Real)

    #--- calculate Rankine-Hugoniot relation ---#
    primL = [1.0, MaL * sqrt(gam / 2.0), 1.0]

    MaR = sqrt((MaL^2 * (gam - 1.0) + 2.0) / (2.0 * gam * MaL^2 - (gam - 1.0)))
	ratioT = (1.0 + (gam - 1.0) / 2.0 * MaL^2) * (2.0 * gam / (gam - 1.0) * MaL^2 - 1.0) / 
		     (MaL^2 * (2.0 * gam / (gam - 1.0) + (gam - 1.0) / 2.0))

    primR = [ primL[1] * (gam + 1.0) * MaL^2 / ((gam - 1.0) * MaL^2 + 2.0),
              MaR * sqrt(gam / 2.0) * sqrt(ratioT),
              primL[3] / ratioT ]

    wL = prim_conserve(primL, gam)
    wR = prim_conserve(primR, gam)

    hL = maxwellian(u, primL)
    hR = maxwellian(u, primR)

    bL = hL .* K ./ (2.0 .* primL[end])
    bR = hR .* K ./ (2.0 .* primR[end])

    bcL = zeros(3)
    bcR = zeros(3)

    return wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR

end


# ------------------------------------------------------------
# Initialize Brio-Wu MHD shock
# ------------------------------------------------------------
function ib_briowu(gam::Float64, uspace::AbstractArray{Float64,2}, mi::Float64, me::Float64)

    mr = mi / me

    # upstream
    primL = Array{Float64}(undef, 5, 2)
    primL[1,1] = 1.
    primL[2,1] = 0.
    primL[3,1] = 0.
    primL[4,1] = 0.
    primL[5,1] = 1.
    primL[1,2] = 1. / mr
    primL[2,2] = 0.
    primL[3,2] = 0.
    primL[4,2] = 0.
    primL[5,2] = 1. / mr

    wL = mixture_prim_conserve(primL, gam)
    h0L = mixture_maxwellian(uspace, primL)
    
    h1L = similar(h0L); h2L = similar(h0L), h3L = similar(h0L)    
    for j in axes(h0L, 2)
        h1L[:,j] .= primL[3,j] .* h0L[:,j]
        h2L[:,j] .= primL[4,j] .* h0L[:,j]
        h3L[:,j] .= (primL[3,j]^2 + primL[4,j]^2 + 2. / (2. * primL[end,j])) * h0L[:,j]
    end

    EL = zeros(3)
    BL = zeros(3); BL[1] = 0.75; BL[2] = 1.

    # downstream
    primR = Array{Float64}(undef, 5, 2)
    primR[1,1] = 0.125
    primR[2,1] = 0.
    primR[3,1] = 0.
    primR[4,1] = 0.
    primR[5,1] = 1.25
    primR[1,2] = 0.125 / mr
    primR[2,2] = 0.
    primR[3,2] = 0.
    primR[4,2] = 0.
    primR[5,2] = 1.25 / mr

    wR = mixture_prim_conserve(primR, gam)
    h0R = mixture_maxwellian(uspace, primR)
    
    h1R = similar(h0R); h2R = similar(h0R), h3R = similar(h0R)    
    for j in axes(h0L, 2)
        h1R[:,j] .= primR[3,j] .* h0R[:,j]
        h2R[:,j] .= primR[4,j] .* h0R[:,j]
        h3R[:,j] .= (primR[3,j]^2 + primR[4,j]^2 + 2. / (2. * primR[end,j])) * h0R[:,j]
    end

    ER = zeros(3)
    BR = zeros(3); BR[1] = 0.75; BR[2] = -1.

    lorenz = zeros(3, 2)
    bc = zeros(5, 2)
    
    return wL, primL, h0L, h1L, h2L, h3L, bc, EL, BL, lorenz, 
           wR, primR, h0R, h1R, h2R, h3R, bc, ER, BR, lorenz

end