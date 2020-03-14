# ============================================================
# Initial and Boundary Conditions of Specific Problems
# ============================================================


export ib_rh


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