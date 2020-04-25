# ============================================================
# Mathematical Methods
# ============================================================


export heaviside,
       fortsign


# ------------------------------------------------------------
# Heaviside step function
# ------------------------------------------------------------
heaviside(x::Real) = ifelse(x >= 0, 1., 0.)


# ------------------------------------------------------------
# Fortran sign() function
# ------------------------------------------------------------
fortsign(x::Real, y::Real) = abs(x) * sign(y)


# ------------------------------------------------------------
# Gauss Legendre integral for fast spectral method
# ------------------------------------------------------------
function lgwt(N::Int, a::Real, b::Real)

    x = zeros(N)
    w = zeros(N)

    N1 = N
    N2 = N + 1

    y = zeros(N1); y0 = zeros(N1); Lp = zeros(N1)
    L = zeros(N1, N2)

    # initial guess
    for i=1:N1
        y[i] = cos( (2. * (i - 1.) + 1.) * 4. * atan(1.) / (2. * (N - 1.) + 2.)  ) +
               0.27 / N1 * sin( 4. * atan(1.) * (-1. + i * 2. / (N1 - 1.)) * (N - 1.) / N2 )
        y0[i] = 2.
    end

    # compute the zeros of the N+1 legendre Polynomial
    # using the recursion relation and the Newton method
    while maximum(abs.(y .- y0)) > 0.0000000000001
        L[:,1] .= 1.0
        L[:,2] .= y
        for k = 2:N1
            @. L[:,k+1] = ( (2. * k - 1.) * y * L[:,k] - (k - 1) * L[:,k-1] ) / k
        end
        @. Lp = N2 * ( L[:,N1] - y * L[:,N2] ) / (1. - y^2)
        @. y0 = y
        @. y = y0 - L[:,N2] / Lp
    end

    # linear map from [-1 1] to [a,b]
    @. x = ( a * (1. - y) + b * (1. + y) ) / 2.0
    @. w = N2^2 * (b - a) / ( (1.0 - y^2) * Lp^2 ) / N1^2

    return x, w

end
