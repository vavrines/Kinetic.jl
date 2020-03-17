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

    wL = prim_conserve(primL, gam)
    hL = mixture_maxwellian(uspace, primL)

    EL = zeros(3)
    BL = zeros(3); BL[1] = 0.75; BL[2] = 1.0

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

    wR = get_conserved(primR, gam)
    hR = get_maxwell(uspace, primR)
    
    ER = zeros(3)
    BR = zeros(3); BR[1] = 0.75; BR[2] = -1.0

    lorenz = zeros(3, 2)
    bc = zeros(5, 2)
    
    return wL, primL, hL, bc, EL, BL, lorenz, 
           wR, primR, hR, bc, ER, BR, lorenz

end