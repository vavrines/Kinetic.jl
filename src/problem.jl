# ============================================================
# Initial and Boundary Conditions of Specific Problems
# ============================================================


export ib_rh, ib_sod, ib_briowu, ib_cavity


```
Initialize Rankine-Hugoniot relation

```

# ------------------------------------------------------------
# 1D1F1V
# ------------------------------------------------------------
function ib_rh(MaL::Real, gam::Real, u::AbstractArray{<:AbstractFloat,1})

    #--- calculate Rankine-Hugoniot relation ---#
    primL = [1.0, MaL * sqrt(gam / 2.0), 1.0]

    MaR = sqrt((MaL^2 * (gam - 1.0) + 2.0) / (2.0 * gam * MaL^2 - (gam - 1.0)))
    ratioT =
        (1.0 + (gam - 1.0) / 2.0 * MaL^2) * (2.0 * gam / (gam - 1.0) * MaL^2 - 1.0) /
        (MaL^2 * (2.0 * gam / (gam - 1.0) + (gam - 1.0) / 2.0))

    primR = [
        primL[1] * (gam + 1.0) * MaL^2 / ((gam - 1.0) * MaL^2 + 2.0),
        MaR * sqrt(gam / 2.0) * sqrt(ratioT),
        primL[3] / ratioT,
    ]

    wL = prim_conserve(primL, gam)
    wR = prim_conserve(primR, gam)

    hL = maxwellian(u, primL)
    hR = maxwellian(u, primR)

    bcL = deepcopy(primL)
    bcR = deepcopy(primR)

    return wL, primL, hL, bcL, wR, primR, hR, bcR

end


# ------------------------------------------------------------
# 1D2F1V
# ------------------------------------------------------------
function ib_rh(MaL::Real, gam::Real, u::AbstractArray{<:AbstractFloat,1}, K::Real)

    #--- calculate Rankine-Hugoniot relation ---#
    primL = [1.0, MaL * sqrt(gam / 2.0), 1.0]

    MaR = sqrt((MaL^2 * (gam - 1.0) + 2.0) / (2.0 * gam * MaL^2 - (gam - 1.0)))
    ratioT =
        (1.0 + (gam - 1.0) / 2.0 * MaL^2) * (2.0 * gam / (gam - 1.0) * MaL^2 - 1.0) /
        (MaL^2 * (2.0 * gam / (gam - 1.0) + (gam - 1.0) / 2.0))

    primR = [
        primL[1] * (gam + 1.0) * MaL^2 / ((gam - 1.0) * MaL^2 + 2.0),
        MaR * sqrt(gam / 2.0) * sqrt(ratioT),
        primL[3] / ratioT,
    ]

    wL = prim_conserve(primL, gam)
    wR = prim_conserve(primR, gam)

    hL = maxwellian(u, primL)
    hR = maxwellian(u, primR)

    bL = hL .* K ./ (2.0 .* primL[end])
    bR = hR .* K ./ (2.0 .* primR[end])

    bcL = deepcopy(primL)
    bcR = deepcopy(primR)

    return wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR

end


# ------------------------------------------------------------
# 1D1F3V
# ------------------------------------------------------------
function ib_rh(
    MaL::Real,
    gam::Real,
    u::AbstractArray{Float64,3},
    v::AbstractArray{Float64,3},
    w::AbstractArray{Float64,3},
)

    #--- calculate Rankine-Hugoniot relation ---#
    primL = [1.0, MaL * sqrt(gam / 2.0), 0.0, 0.0, 1.0]

    MaR = sqrt((MaL^2 * (gam - 1.0) + 2.0) / (2.0 * gam * MaL^2 - (gam - 1.0)))
    ratioT =
        (1.0 + (gam - 1.0) / 2.0 * MaL^2) * (2.0 * gam / (gam - 1.0) * MaL^2 - 1.0) /
        (MaL^2 * (2.0 * gam / (gam - 1.0) + (gam - 1.0) / 2.0))

    primR = [
        primL[1] * (gam + 1.0) * MaL^2 / ((gam - 1.0) * MaL^2 + 2.0),
        MaR * sqrt(gam / 2.0) * sqrt(ratioT),
        0.0,
        0.0,
        primL[end] / ratioT,
    ]

    wL = prim_conserve(primL, gam)
    wR = prim_conserve(primR, gam)

    fL = maxwellian(u, v, w, primL)
    fR = maxwellian(u, v, w, primR)

    bcL = deepcopy(primL)
    bcR = deepcopy(primR)

    return wL, primL, fL, bcL, wR, primR, fR, bcR

end


```
Initialize Sod shock tube

```

# ------------------------------------------------------------
# 1D1F1V
# ------------------------------------------------------------
function ib_sod(γ::Real, u::AbstractArray{<:AbstractFloat,1})

    primL = [1.0, 0.0, 1.0]
    primR = [0.125, 0.0, 0.625]

    wL = prim_conserve(primL, γ)
    wR = prim_conserve(primR, γ)

    fL = maxwellian(u, primL)
    fR = maxwellian(u, primR)

    bcL = deepcopy(primL)
    bcR = deepcopy(primR)

    return wL, primL, fL, bcL, wR, primR, fR, bcR

end


# ------------------------------------------------------------
# 1D2F1V
# ------------------------------------------------------------
function ib_sod(γ::Real, u::AbstractArray{<:AbstractFloat,1}, K::Real)

    primL = [1.0, 0.0, 1.0]
    primR = [0.125, 0.0, 0.625]

    wL = prim_conserve(primL, γ)
    wR = prim_conserve(primR, γ)

    hL = maxwellian(u, primL)
    hR = maxwellian(u, primR)

    bL = hL .* K ./ (2.0 .* primL[end])
    bR = hR .* K ./ (2.0 .* primR[end])

    bcL = deepcopy(primL)
    bcR = deepcopy(primR)

    return wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR

end


# ------------------------------------------------------------
# 1D1F3V
# ------------------------------------------------------------
function ib_rh(
    γ::Real,
    u::AbstractArray{Float64,3},
    v::AbstractArray{Float64,3},
    w::AbstractArray{Float64,3},
)

    primL = [1.0, 0.0, 1.0]
    primR = [0.125, 0.0, 0.625]

    wL = prim_conserve(primL, γ)
    wR = prim_conserve(primR, γ)

    fL = maxwellian(u, v, w, primL)
    fR = maxwellian(u, v, w, primR)

    bcL = deepcopy(primL)
    bcR = deepcopy(primR)

    return wL, primL, fL, bcL, wR, primR, fR, bcR

end


```
Initialize lid-driven cavity

```

# ------------------------------------------------------------
# 2D1F2V
# ------------------------------------------------------------
function ib_cavity(
    gam::Real,
    Um::Real,
    Vm::Real,
    Tm::Real,
    u::AbstractArray{Float64,2},
    v::AbstractArray{Float64,2},
)

    primL = [1.0, 0.0, 0.0, 1.0]
    primR = deepcopy(primL)

    wL = prim_conserve(primL, gam)
    wR = prim_conserve(primR, gam)

    fL = maxwellian(u, v, primL)
    fR = maxwellian(u, v, primR)

    bcU = [1.0, Um, Vm, Tm]
    bcD = deepcopy(primR)
    bcL = deepcopy(primR)
    bcR = deepcopy(primR)

    return wL, primL, fL, bcL, wR, primR, fR, bcR, bcU, bcD

end

# ------------------------------------------------------------
# 2D2F2V
# ------------------------------------------------------------
function ib_cavity(
    gam::Real,
    Um::Real,
    Vm::Real,
    Tm::Real,
    u::AbstractArray{Float64,2},
    v::AbstractArray{Float64,2},
    K::Real,
)

    primL = [1.0, 0.0, 0.0, 1.0]
    primR = deepcopy(primL)

    wL = prim_conserve(primL, gam)
    wR = prim_conserve(primR, gam)

    hL = maxwellian(u, v, primL)
    hR = maxwellian(u, v, primR)

    bL = hL .* K ./ (2.0 * primL[end])
    bR = hR .* K ./ (2.0 * primR[end])

    bcU = [1.0, Um, Vm, Tm]
    bcD = deepcopy(primR)
    bcL = deepcopy(primR)
    bcR = deepcopy(primR)

    return wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR, bcU, bcD

end


```
Initialize Brio-Wu MHD shock

```

function ib_briowu(gam::Real, uspace::AbstractArray{Float64,2}, mi::Real, me::Real)

    # upstream
    primL = zeros(5, 2)
    primL[1, 1] = 1.0 * mi
    primL[2, 1] = 0.0
    primL[3, 1] = 0.0
    primL[4, 1] = 0.0
    primL[5, 1] = mi / 1.0
    primL[1, 2] = 1.0 * me
    primL[2, 2] = 0.0
    primL[3, 2] = 0.0
    primL[4, 2] = 0.0
    primL[5, 2] = me / 1.0

    wL = mixture_prim_conserve(primL, gam)
    h0L = mixture_maxwellian(uspace, primL)

    h1L = similar(h0L)
    h2L = similar(h0L)
    h3L = similar(h0L)
    for j in axes(h0L, 2)
        h1L[:, j] .= primL[3, j] .* h0L[:, j]
        h2L[:, j] .= primL[4, j] .* h0L[:, j]
        h3L[:, j] .=
            (primL[3, j]^2 + primL[4, j]^2 + 2.0 / (2.0 * primL[end, j])) .* h0L[:, j]
    end

    EL = zeros(3)
    BL = zeros(3)
    BL[1] = 0.75
    BL[2] = 1.0

    # downstream
    primR = zeros(5, 2)
    primR[1, 1] = 0.125 * mi
    primR[2, 1] = 0.0
    primR[3, 1] = 0.0
    primR[4, 1] = 0.0
    primR[5, 1] = mi * 1.25
    primR[1, 2] = 0.125 * me
    primR[2, 2] = 0.0
    primR[3, 2] = 0.0
    primR[4, 2] = 0.0
    primR[5, 2] = me * 1.25

    wR = mixture_prim_conserve(primR, gam)
    h0R = mixture_maxwellian(uspace, primR)

    h1R = similar(h0R)
    h2R = similar(h0R)
    h3R = similar(h0R)
    for j in axes(h0L, 2)
        h1R[:, j] .= primR[3, j] .* h0R[:, j]
        h2R[:, j] .= primR[4, j] .* h0R[:, j]
        h3R[:, j] .=
            (primR[3, j]^2 + primR[4, j]^2 + 2.0 / (2.0 * primR[end, j])) .* h0R[:, j]
    end

    ER = zeros(3)
    BR = zeros(3)
    BR[1] = 0.75
    BR[2] = -1.0

    lorenzL = zeros(3, 2)
    lorenzR = zeros(3, 2)
    bcL = zeros(5, 2)
    bcR = zeros(5, 2)

    return wL,
    primL,
    h0L,
    h1L,
    h2L,
    h3L,
    bcL,
    EL,
    BL,
    lorenzL,
    wR,
    primR,
    h0R,
    h1R,
    h2R,
    h3R,
    bcR,
    ER,
    BR,
    lorenzR

end
