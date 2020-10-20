"""
Kinetic flux vector splitting (KFVS) flux

* DOM: `flux_kfvs!(ff, fL, fR, u, dt, sfL, sfR)`
* 1D1F1V: `flux_kfvs!(fw, ff, fL, fR, u, ω, dt, sfL, sfR)`
* 1D1F3V: `flux_kfvs!(fw, ff, fL, fR, u, v, w, ω, dt, sfL, sfR)`
* 1D2F1V: `flux_kfvs!(fw, fh, fb, hL, bL, hR, bR, u, ω, dt, shL, sbL, shR, sbR)`
* 1D4F1V: `flux_kfvs!(fw, fh0, fh1, fh2, fh3, h0L, h1L, h2L, h3L, h0R, h1R, h2R, h3R, u, ω, dt, sh0L, sh1L, sh2L, sh3L, sh0R, sh1R, sh2R, sh3R)`
* 2D1F2V: `flux_kfvs!(fw, ff, fL, fR, u, v, ω, dt, len, sfL, sfR)`
* 2D2F2V: `flux_kfvs!(fw, fh, fb, hL, bL, hR, bR, u, v, ω, dt, len, shL, sbL, shR, sbR)`

* @arg: particle distribution functions and their left/right slopes
* @arg: particle velocity quadrature points and weights
* @arg: time step and cell size

"""
function flux_kfvs!(
    ff::X,
    fL::Y,
    fR::Y,
    u::Z,
    dt,
    sfL = zeros(eltype(fL), axes(fL))::Y,
    sfR = zeros(eltype(fR), axes(fR))::Y,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,1},
    Z<:AbstractArray{<:AbstractFloat,1},
} # 1F1V flux of pure DOM

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    f = @. fL * δ + fR * (1.0 - δ)
    sf = @. sfL * δ + sfR * (1.0 - δ)

    # --- calculate fluxes ---#
    @. ff = dt * u * f - 0.5 * dt^2 * u^2 * sf

    return nothing

end

# ------------------------------------------------------------
# 1F1V flux
# ------------------------------------------------------------
function flux_kfvs!(
    fw::X,
    ff::Y,
    fL::Z,
    fR::Z,
    u::A,
    ω::A,
    dt,
    sfL = zeros(eltype(fL), axes(fL))::Z,
    sfR = zeros(eltype(fR), axes(fR))::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,1},
    Z<:AbstractArray{<:AbstractFloat,1},
    A<:AbstractArray{<:AbstractFloat,1},
}

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    f = @. fL * δ + fR * (1.0 - δ)
    sf = @. sfL * δ + sfR * (1.0 - δ)

    # --- calculate fluxes ---#
    fw[1] = dt * sum(ω .* u .* f) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sf)
    fw[2] = dt * sum(ω .* u .^ 2 .* f) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sf)
    fw[3] = dt * 0.5 * sum(ω .* u .^ 3 .* f) - 0.5 * dt^2 * 0.5 * sum(ω .* u .^ 4 .* sf)

    @. ff = dt * u * f - 0.5 * dt^2 * u^2 * sf

    return nothing

end

# ------------------------------------------------------------
# 1F3V flux
# ------------------------------------------------------------
function flux_kfvs!(
    fw::X,
    ff::Y,
    fL::Z,
    fR::Z,
    u::A,
    v::A,
    w::A,
    ω::A,
    dt,
    sfL = zeros(eltype(fL), axes(fL))::Z,
    sfR = zeros(eltype(fR), axes(fR))::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,3},
    Z<:AbstractArray{<:AbstractFloat,3},
    A<:AbstractArray{<:AbstractFloat,3},
}

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    f = @. fL * δ + fR * (1.0 - δ)
    sf = @. sfL * δ + sfR * (1.0 - δ)

    # --- calculate fluxes ---#
    fw[1] = dt * sum(ω .* u .* f) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sf)
    fw[2] = dt * sum(ω .* u .^ 2 .* f) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sf)
    fw[3] = dt * sum(ω .* u .* v .* f) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* v .* sf)
    fw[4] = dt * sum(ω .* u .* w .* f) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* w .* sf)
    fw[5] =
        dt * 0.5 * sum(ω .* u .* (u .^ 2 .+ v .^ 2 .+ w .^ 2) .* f) -
        0.5 * dt^2 * 0.5 * sum(ω .* u .^ 2 .* (u .^ 2 .+ v .^ 2 .+ w .^ 2) .* sf)

    @. ff = dt * u * f - 0.5 * dt^2 * u^2 * sf

    return nothing

end

# ------------------------------------------------------------
# 2F1V flux
# ------------------------------------------------------------
function flux_kfvs!(
    fw::X,
    fh::Y,
    fb::Y,
    hL::Z,
    bL::Z,
    hR::Z,
    bR::Z,
    u::A,
    ω::A,
    dt,
    shL = zeros(eltype(hL), axes(hL))::Z,
    sbL = zeros(eltype(bL), axes(bL))::Z,
    shR = zeros(eltype(hR), axes(hR))::Z,
    sbR = zeros(eltype(bR), axes(bR))::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,1},
    Z<:AbstractArray{<:AbstractFloat,1},
    A<:AbstractArray{<:AbstractFloat,1},
}

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    h = @. hL * δ + hR * (1.0 - δ)
    b = @. bL * δ + bR * (1.0 - δ)

    sh = @. shL * δ + shR * (1.0 - δ)
    sb = @. sbL * δ + sbR * (1.0 - δ)

    # --- calculate fluxes ---#
    fw[1] = dt * sum(ω .* u .* h) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh)
    fw[2] = dt * sum(ω .* u .^ 2 .* h) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sh)
    fw[3] =
        dt * 0.5 * (sum(ω .* u .^ 3 .* h) + sum(ω .* u .* b)) -
        0.5 * dt^2 * 0.5 * (sum(ω .* u .^ 4 .* sh) + sum(ω .* u .^ 2 .* sb))

    @. fh = dt * u * h - 0.5 * dt^2 * u^2 * sh
    @. fb = dt * u * b - 0.5 * dt^2 * u^2 * sb

    return nothing

end

# ------------------------------------------------------------
# 4F1V flux (single component)
# ------------------------------------------------------------
function flux_kfvs!(
    fw::X,
    fh0::Y,
    fh1::Y,
    fh2::Y,
    fh3::Y,
    h0L::Z,
    h1L::Z,
    h2L::Z,
    h3L::Z,
    h0R::Z,
    h1R::Z,
    h2R::Z,
    h3R::Z,
    u::A,
    ω::A,
    dt,
    sh0L = zeros(eltype(h0L), axes(h0L))::Z,
    sh1L = zeros(eltype(h1L), axes(h1L))::Z,
    sh2L = zeros(eltype(h2L), axes(h2L))::Z,
    sh3L = zeros(eltype(h3L), axes(h3L))::Z,
    sh0R = zeros(eltype(h0R), axes(h0R))::Z,
    sh1R = zeros(eltype(h1R), axes(h1R))::Z,
    sh2R = zeros(eltype(h2R), axes(h2R))::Z,
    sh3R = zeros(eltype(h3R), axes(h3R))::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,1},
    Z<:AbstractArray{<:AbstractFloat,1},
    A<:AbstractArray{<:AbstractFloat,1},
}

    #--- upwind reconstruction ---#
    δ = heaviside.(u)

    h0 = @. h0L * δ + h0R * (1.0 - δ)
    h1 = @. h1L * δ + h1R * (1.0 - δ)
    h2 = @. h2L * δ + h2R * (1.0 - δ)
    h3 = @. h3L * δ + h3R * (1.0 - δ)

    sh0 = @. sh0L * δ + sh0R * (1.0 - δ)
    sh1 = @. sh1L * δ + sh1R * (1.0 - δ)
    sh2 = @. sh2L * δ + sh2R * (1.0 - δ)
    sh3 = @. sh3L * δ + sh3R * (1.0 - δ)

    # --- calculate fluxes ---#
    fw[1] = dt * sum(ω .* u .* h0) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh0)
    fw[2] = dt * sum(ω .* u .^ 2 .* h0) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sh0)
    fw[3] = dt * sum(ω .* u .* h1) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh1)
    fw[4] = dt * sum(ω .* u .* h2) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh2)
    fw[5] =
        dt * 0.5 * (sum(ω .* u .^ 3 .* h0) + sum(ω .* u .* h3)) -
        0.5 * dt^2 * 0.5 * (sum(ω .* u .^ 4 .* sh0) + sum(ω .* u .^ 2 .* sh3))

    @. fh0 = dt * u * h0 - 0.5 * dt^2 * u^2 * sh0
    @. fh1 = dt * u * h1 - 0.5 * dt^2 * u^2 * sh1
    @. fh2 = dt * u * h2 - 0.5 * dt^2 * u^2 * sh2
    @. fh3 = dt * u * h3 - 0.5 * dt^2 * u^2 * sh3

    return nothing

end

# ------------------------------------------------------------
# 4F1V flux (multiple component)
# ------------------------------------------------------------
function flux_kfvs!(
    fw::X,
    fh0::Y,
    fh1::Y,
    fh2::Y,
    fh3::Y,
    h0L::Z,
    h1L::Z,
    h2L::Z,
    h3L::Z,
    h0R::Z,
    h1R::Z,
    h2R::Z,
    h3R::Z,
    u::A,
    ω::A,
    dt,
    sh0L = zeros(eltype(h0L), axes(h0L))::Z,
    sh1L = zeros(eltype(h1L), axes(h1L))::Z,
    sh2L = zeros(eltype(h2L), axes(h2L))::Z,
    sh3L = zeros(eltype(h3L), axes(h3L))::Z,
    sh0R = zeros(eltype(h0R), axes(h0R))::Z,
    sh1R = zeros(eltype(h1R), axes(h1R))::Z,
    sh2R = zeros(eltype(h2R), axes(h2R))::Z,
    sh3R = zeros(eltype(h3R), axes(h3R))::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,2},
    Y<:AbstractArray{<:AbstractFloat,2},
    Z<:AbstractArray{<:AbstractFloat,2},
    A<:AbstractArray{<:AbstractFloat,2},
}

    for j in axes(fw, 2)
        _fw = @view fw[:, j]
        _fh0 = @view fh0[:, j]
        _fh1 = @view fh1[:, j]
        _fh2 = @view fh2[:, j]
        _fh3 = @view fh3[:, j]

        flux_kfvs!(
            _fw,
            _fh0,
            _fh1,
            _fh2,
            _fh3,
            h0L[:, j],
            h1L[:, j],
            h2L[:, j],
            h3L[:, j],
            h0R[:, j],
            h1R[:, j],
            h2R[:, j],
            h3R[:, j],
            u[:, j],
            ω[:, j],
            dt,
            sh0L[:, j],
            sh1L[:, j],
            sh2L[:, j],
            sh3L[:, j],
            sh0R[:, j],
            sh1R[:, j],
            sh2R[:, j],
            sh3R[:, j],
        )
    end

    return nothing

end

# ------------------------------------------------------------
# 1F2V flux
# ------------------------------------------------------------
function flux_kfvs!(
    fw::X,
    ff::Y,
    fL::Z,
    fR::Z,
    u::A,
    v::A,
    ω::A,
    dt,
    len,
    sfL = zeros(eltype(fL), axes(fL))::Z,
    sfR = zeros(eltype(fR), axes(fR))::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,2},
    Z<:AbstractArray{<:AbstractFloat,2},
    A<:AbstractArray{<:AbstractFloat,2},
}

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    f = @. fL * δ + fR * (1.0 - δ)
    sf = @. sfL * δ + sfR * (1.0 - δ)

    # --- calculate fluxes ---#
    fw[1] = dt * sum(ω .* u .* f) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sf)
    fw[2] = dt * sum(ω .* u .^ 2 .* f) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sf)
    fw[3] = dt * sum(ω .* v .* u .* f) - 0.5 * dt^2 * sum(ω .* v .* u .^ 2 .* sf)
    fw[4] =
        dt * 0.5 * sum(ω .* u .* (u .^ 2 .+ v .^ 2) .* f) -
        0.5 * dt^2 * 0.5 * sum(ω .* u .^ 2 .* (u .^ 2 .+ v .^ 2) .* sf)
    fw .*= len

    @. ff = (dt * u * f - 0.5 * dt^2 * u^2 * sf) * len

    return nothing

end

# ------------------------------------------------------------
# 2F2V flux
# ------------------------------------------------------------
function flux_kfvs!(
    fw::X,
    fh::Y,
    fb::Y,
    hL::Z,
    bL::Z,
    hR::Z,
    bR::Z,
    u::A,
    v::A,
    ω::A,
    dt,
    len,
    shL = zeros(eltype(hL), axes(hL))::Z,
    sbL = zeros(eltype(bL), axes(bL))::Z,
    shR = zeros(eltype(hR), axes(hR))::Z,
    sbR = zeros(eltype(bR), axes(bR))::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,2},
    Z<:AbstractArray{<:AbstractFloat,2},
    A<:AbstractArray{<:AbstractFloat,2},
}

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    h = @. hL * δ + hR * (1.0 - δ)
    b = @. bL * δ + bR * (1.0 - δ)
    sh = @. shL * δ + shR * (1.0 - δ)
    sb = @. sbL * δ + sbR * (1.0 - δ)

    # --- calculate fluxes ---#
    fw[1] = dt * sum(ω .* u .* h) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh)
    fw[2] = dt * sum(ω .* u .^ 2 .* h) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sh)
    fw[3] = dt * sum(ω .* v .* u .* h) - 0.5 * dt^2 * sum(ω .* v .* u .^ 2 .* sh)
    fw[4] =
        dt * 0.5 * (sum(ω .* u .* (u .^ 2 .+ v .^ 2) .* h) + sum(ω .* u .* b)) -
        0.5 *
        dt^2 *
        0.5 *
        (sum(ω .* u .^ 2 .* (u .^ 2 .+ v .^ 2) .* sh) + sum(ω .* u .^ 2 .* sb))
    fw .*= len

    @. fh = (dt * u * h - 0.5 * dt^2 * u^2 * sh) * len
    @. fb = (dt * u * b - 0.5 * dt^2 * u^2 * sb) * len

    return nothing

end

# ------------------------------------------------------------
# 3F2V flux (single component)
# ------------------------------------------------------------
function flux_kfvs!(
    fw::X,
    fh0::Y,
    fh1::Y,
    fh2::Y,
    h0L::Z,
    h1L::Z,
    h2L::Z,
    h0R::Z,
    h1R::Z,
    h2R::Z,
    u::A,
    v::A,
    ω::A,
    dt,
    len,
    sh0L = zeros(eltype(h0L), axes(h0L))::Z,
    sh1L = zeros(eltype(h1L), axes(h1L))::Z,
    sh2L = zeros(eltype(h2L), axes(h2L))::Z,
    sh0R = zeros(eltype(h0R), axes(h0R))::Z,
    sh1R = zeros(eltype(h1R), axes(h1R))::Z,
    sh2R = zeros(eltype(h2R), axes(h2R))::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,2},
    Z<:AbstractArray{<:AbstractFloat,2},
    A<:AbstractArray{<:AbstractFloat,2},
}

    #--- upwind reconstruction ---#
    δ = heaviside.(u)

    h0 = @. h0L * δ + h0R * (1.0 - δ)
    h1 = @. h1L * δ + h1R * (1.0 - δ)
    h2 = @. h2L * δ + h2R * (1.0 - δ)
    sh0 = @. sh0L * δ + sh0R * (1.0 - δ)
    sh1 = @. sh1L * δ + sh1R * (1.0 - δ)
    sh2 = @. sh2L * δ + sh2R * (1.0 - δ)

    #--- calculate fluxes ---#
    fw[1] = dt * sum(ω .* u .* h0) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh0)
    fw[2] = dt * sum(ω .* u .^ 2 .* h0) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sh0)
    fw[3] = dt * sum(ω .* u .* v .* h0) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* v .* sh0)
    fw[4] = dt * sum(ω .* u .* h1) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh1)
    fw[5] =
        dt * 0.5 * (sum(ω .* u .* (u .^ 2 .+ v .^ 2) .* h0) + sum(ω .* u .* h2)) -
        0.5 *
        dt^2 *
        0.5 *
        (sum(ω .* u .^ 2 .* (u .^ 2 .+ v .^ 2) .* sh0) + sum(ω .* u .^ 2 .* sh2))

    fw .*= len
    @. fh0 = (dt * u * h0 - 0.5 * dt^2 * u^2 * sh0) * len
    @. fh1 = (dt * u * h1 - 0.5 * dt^2 * u^2 * sh1) * len
    @. fh2 = (dt * u * h2 - 0.5 * dt^2 * u^2 * sh2) * len

    return nothing

end

# ------------------------------------------------------------
# 3F2V flux (multiple component)
# ------------------------------------------------------------
function flux_kfvs!(
    fw::X,
    fh0::Y,
    fh1::Y,
    fh2::Y,
    h0L::Z,
    h1L::Z,
    h2L::Z,
    h0R::Z,
    h1R::Z,
    h2R::Z,
    u::A,
    v::A,
    ω::A,
    dt,
    len,
    sh0L = zeros(eltype(h0L), axes(h0L))::Z,
    sh1L = zeros(eltype(h1L), axes(h1L))::Z,
    sh2L = zeros(eltype(h2L), axes(h2L))::Z,
    sh0R = zeros(eltype(h0R), axes(h0R))::Z,
    sh1R = zeros(eltype(h1R), axes(h1R))::Z,
    sh2R = zeros(eltype(h2R), axes(h2R))::Z,
) where {
    X<:AbstractArray{<:AbstractFloat,2},
    Y<:AbstractArray{<:AbstractFloat,3},
    Z<:AbstractArray{<:AbstractFloat,3},
    A<:AbstractArray{<:AbstractFloat,3},
}

    #--- reconstruct initial distribution ---#
    δ = heaviside.(u)

    h0 = @. h0L * δ + h0R * (1.0 - δ)
    h1 = @. h1L * δ + h1R * (1.0 - δ)
    h2 = @. h2L * δ + h2R * (1.0 - δ)

    sh0 = @. sh0L * δ + sh0R * (1.0 - δ)
    sh1 = @. sh1L * δ + sh1R * (1.0 - δ)
    sh2 = @. sh2L * δ + sh2R * (1.0 - δ)

    for j = 1:2
        fw[1, j] =
            dt * sum(ω[:, :, j] .* u[:, :, j] .* h0[:, :, j]) -
            0.5 * dt^2 * sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* sh0[:, :, j])
        fw[2, j] =
            dt * sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* h0[:, :, j]) -
            0.5 * dt^2 * sum(ω[:, :, j] .* u[:, :, j] .^ 3 .* sh0[:, :, j])
        fw[3, j] =
            dt * sum(ω[:, :, j] .* v[:, :, j] .* u[:, :, j] .* h0[:, :, j]) -
            0.5 * dt^2 * sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* v[:, :, j] .* sh0[:, :, j])
        fw[4, j] =
            dt * sum(ω[:, :, j] .* u[:, :, j] .* h1[:, :, j]) -
            0.5 * dt^2 * sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* sh1[:, :, j])
        fw[5, j] =
            dt *
            0.5 *
            (
                sum(
                    ω[:, :, j] .* u[:, :, j] .* (u[:, :, j] .^ 2 .+ v[:, :, j] .^ 2) .*
                    h0[:, :, j],
                ) + sum(ω[:, :, j] .* u[:, :, j] .* h2[:, :, j])
            ) -
            0.5 *
            dt^2 *
            0.5 *
            (
                sum(
                    ω[:, :, j] .* u[:, :, j] .^ 2 .* (u[:, :, j] .^ 2 .+ v[:, :, j] .^ 2) .*
                    sh0[:, :, j],
                ) + sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* sh2[:, :, j])
            )

        @. fh0[:, :, j] =
            dt * u[:, :, j] * h0[:, :, j] - 0.5 * dt^2 * u[:, :, j]^2 * sh0[:, :, j]
        @. fh1[:, :, j] =
            dt * u[:, :, j] * h1[:, :, j] - 0.5 * dt^2 * u[:, :, j]^2 * sh1[:, :, j]
        @. fh2[:, :, j] =
            dt * u[:, :, j] * h2[:, :, j] - 0.5 * dt^2 * u[:, :, j]^2 * sh2[:, :, j]
    end

    @. fw *= len
    @. fh0 *= len
    @. fh1 *= len
    @. fh2 *= len

    return nothing
    
end
