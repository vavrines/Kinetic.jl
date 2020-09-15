# ============================================================
# General Structs
# ============================================================


"""
Structure of computational setup

"""
struct Setup{S,I,E,F} <: AbstractSetup

    case::S
    space::S
    nSpecies::I
    interpOrder::I
    limiter::S
    cfl::E
    maxTime::F

    function Setup(
        case::AbstractString,
        space::AbstractString,
        nSpecies::Int,
        interpOrder::Int,
        limiter::AbstractString,
        cfl::Real,
        maxTime::Real,
    )
        new{typeof(case),typeof(nSpecies),typeof(cfl),typeof(maxTime)}(
            case,
            space,
            nSpecies,
            interpOrder,
            limiter,
            cfl,
            maxTime,
        )
    end

end


"""
Structure of gas property

"""
struct Gas{A,B,C,D,E,F,G,H,I} <: AbstractProperty

    Kn::A
    Ma::B
    Pr::C
    K::D
    γ::E
    ω::F
    αᵣ::G
    ωᵣ::H
    μᵣ::I

    function Gas(
        Kn::Union{Real,AbstractArray}, # unified consideration of
        Ma::Union{Real,AbstractArray}, # 1. deterministic solution, and
        Pr::Union{Real,AbstractArray}, # 2. uncertainty quantification
        K::Union{Real,AbstractArray},
        γ::Union{Real,AbstractArray},
        ω::Union{Real,AbstractArray},
        αᵣ::Union{Real,AbstractArray},
        ωᵣ::Union{Real,AbstractArray},
        μᵣ::Union{Real,AbstractArray},
    )
        new{
            typeof(Kn),
            typeof(Ma),
            typeof(Pr),
            typeof(K),
            typeof(γ),
            typeof(ω),
            typeof(αᵣ),
            typeof(ωᵣ),
            typeof(μᵣ),
        }(
            Kn,
            Ma,
            Pr,
            K,
            γ,
            ω,
            αᵣ,
            ωᵣ,
            μᵣ,
        )
    end

end


"""
Structure of 1D plasma property

"""
struct Plasma1D{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P} <: AbstractProperty

    Kn::A
    Ma::B
    Pr::C
    K::D
    γ::E

    mi::F
    ni::G
    me::H
    ne::I
    lD::J
    rL::K

    sol::L
    χ::M
    ν::N
    Ap::O
    An::O
    D::P

    # unified consideration of deterministic and stochastic conditions
    # {mi, ni, me, ne, sol, χ, ν} keep the same as before
    function Plasma1D(
        Kn::Array,
        Ma::Union{Real,AbstractArray},
        Pr::Union{Real,AbstractArray},
        K::Union{Real,AbstractArray},
        γ::Union{Real,AbstractArray},
        mi::Real,
        ni::Real,
        me::Real,
        ne::Real,
        lD::Union{Real,AbstractArray},
        rL::Union{Real,AbstractArray},
        sol::Real,
        χ::Real,
        ν::Real,
    )
        # A^+
        Ap = Array{Float64}(undef, 8, 8)
        Ap[1, 1] = (sol * χ) / 2.0
        Ap[7, 1] = χ / 2.0
        Ap[2, 2] = sol / 2.0
        Ap[6, 2] = 0.5
        Ap[3, 3] = sol / 2.0
        Ap[5, 3] = -1.0 / 2.0
        Ap[4, 4] = (sol * ν) / 2.0
        Ap[8, 4] = (sol^2 * ν) / 2.0
        Ap[3, 5] = -sol^2 / 2.0
        Ap[5, 5] = sol / 2.0
        Ap[2, 6] = sol^2 / 2.0
        Ap[6, 6] = sol / 2.0
        Ap[1, 7] = (sol^2 * χ) / 2.0
        Ap[7, 7] = (sol * χ) / 2.0
        Ap[4, 8] = ν / 2.0
        Ap[8, 8] = (sol * ν) / 2.0

        # A^-
        An = Array{Float64}(undef, 8, 8)
        An[1, 1] = -(sol * χ) / 2.0
        An[7, 1] = χ / 2.0
        An[2, 2] = -sol / 2.0
        An[6, 2] = 1.0 / 2.0
        An[3, 3] = -sol / 2.0
        An[5, 3] = -1.0 / 2.0
        An[4, 4] = -(sol * ν) / 2.0
        An[8, 4] = (sol^2 * ν) / 2.0
        An[3, 5] = -sol^2 / 2.0
        An[5, 5] = -sol / 2.0
        An[2, 6] = sol^2 / 2.0
        An[6, 6] = -sol / 2.0
        An[1, 7] = (sol^2 * χ) / 2.0
        An[7, 7] = -(sol * χ) / 2.0
        An[4, 8] = ν / 2.0
        An[8, 8] = -(sol * ν) / 2.0

        # eigenvalues of A
        D = [sol, sol, sol * χ, sol * ν, -sol, -sol, -sol * χ, -sol * ν]

        # inner constructor method
        new{
            typeof(Kn),
            typeof(Ma),
            typeof(Pr),
            typeof(K),
            typeof(γ),
            typeof(mi),
            typeof(ni),
            typeof(me),
            typeof(ne),
            typeof(lD),
            typeof(rL),
            typeof(sol),
            typeof(χ),
            typeof(ν),
            typeof(Ap),
            typeof(D),
        }(
            Kn,
            Ma,
            Pr,
            K,
            γ,
            mi,
            ni,
            me,
            ne,
            lD,
            rL,
            sol,
            χ,
            ν,
            Ap,
            An,
            D,
        )
    end

end


"""
Structure of 2D plasma property

"""
struct Plasma2D{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P} <: AbstractProperty

    Kn::A
    Ma::B
    Pr::C
    K::D
    γ::E

    mi::F
    ni::G
    me::H
    ne::I
    lD::J
    rL::K

    sol::L
    χ::M
    ν::N
    A1p::O
    A1n::O
    A2p::O
    A2n::O
    D1::P
    D2::P

    # unified consideration of deterministic and stochastic conditions
    # {mi, ni, me, ne, sol, χ, ν} keep the same as before
    function Plasma2D(
        Kn::Array,
        Ma::Union{Real,AbstractArray},
        Pr::Union{Real,AbstractArray},
        K::Union{Real,AbstractArray},
        γ::Union{Real,AbstractArray},
        mi::Real,
        ni::Real,
        me::Real,
        ne::Real,
        lD::Union{Real,AbstractArray},
        rL::Union{Real,AbstractArray},
        sol::Real,
        χ::Real,
        ν::Real,
    )
        # A₁+
        A1p = zeros(8, 8)
        A1p[1, 1] = (sol * χ) / 2.0
        A1p[7, 1] = χ / 2.0
        A1p[2, 2] = sol / 2.0
        A1p[6, 2] = 0.5
        A1p[3, 3] = sol / 2.0
        A1p[5, 3] = -1.0 / 2.0
        A1p[4, 4] = (sol * ν) / 2.0
        A1p[8, 4] = (sol^2 * ν) / 2.0
        A1p[3, 5] = -sol^2 / 2.0
        A1p[5, 5] = sol / 2.0
        A1p[2, 6] = sol^2 / 2.0
        A1p[6, 6] = sol / 2.0
        A1p[1, 7] = (sol^2 * χ) / 2.0
        A1p[7, 7] = (sol * χ) / 2.0
        A1p[4, 8] = ν / 2.0
        A1p[8, 8] = (sol * ν) / 2.0

        # A₁-
        A1n = zeros(8, 8)
        A1n[1, 1] = -(sol * χ) / 2.0
        A1n[7, 1] = χ / 2.0
        A1n[2, 2] = -sol / 2.0
        A1n[6, 2] = 1.0 / 2.0
        A1n[3, 3] = -sol / 2.0
        A1n[5, 3] = -1.0 / 2.0
        A1n[4, 4] = -(sol * ν) / 2.0
        A1n[8, 4] = (sol^2 * ν) / 2.0
        A1n[3, 5] = -sol^2 / 2.0
        A1n[5, 5] = -sol / 2.0
        A1n[2, 6] = sol^2 / 2.0
        A1n[6, 6] = -sol / 2.0
        A1n[1, 7] = (sol^2 * χ) / 2.0
        A1n[7, 7] = -(sol * χ) / 2.0
        A1n[4, 8] = ν / 2.0
        A1n[8, 8] = -(sol * ν) / 2.0

        D1 = [sol, sol, sol * χ, sol * ν, -sol, -sol, -sol * χ, -sol * ν]

        # A₂+
        A2p = zeros(8, 8)
        A2p[1, 1] = sol / 2.0
        A2p[6, 1] = -1.0 / 2.0
        A2p[2, 2] = (sol * χ) / 2.0
        A2p[7, 2] = χ / 2.0
        A2p[3, 3] = sol / 2.0
        A2p[4, 3] = 1.0 / 2.0
        A2p[3, 4] = (sol^2) / 2.0
        A2p[4, 4] = sol / 2.0
        A2p[5, 5] = (sol * ν) / 2.0
        A2p[8, 5] = (sol^2 * ν) / 2.0
        A2p[1, 6] = -sol^2 / 2.0
        A2p[6, 6] = sol / 2.0
        A2p[2, 7] = (sol^2 * χ) / 2.0
        A2p[7, 7] = (sol * χ) / 2.0
        A2p[5, 8] = ν / 2.0
        A2p[8, 8] = (sol * ν) / 2.0

        # A₂-
        A2n = zeros(8, 8)
        A2n[1, 1] = -sol / 2.0
        A2n[6, 1] = -1.0 / 2.0
        A2n[2, 2] = -(sol * χ) / 2.0
        A2n[7, 2] = χ / 2.0
        A2n[3, 3] = -sol / 2.0
        A2n[4, 3] = 1.0 / 2.0
        A2n[3, 4] = (sol^2) / 2.0
        A2n[4, 4] = -sol / 2.0
        A2n[5, 5] = -(sol * ν) / 2.0
        A2n[8, 5] = (sol^2 * ν) / 2.0
        A2n[1, 6] = -sol^2 / 2.0
        A2n[6, 6] = -sol / 2.0
        A2n[2, 7] = (sol^2 * χ) / 2.0
        A2n[7, 7] = -(sol * χ) / 2.0
        A2n[5, 8] = ν / 2.0
        A2n[8, 8] = -(sol * ν) / 2.0

        D2 = [sol, sol, sol * χ, sol * ν, -sol, -sol, -sol * χ, -sol * ν]

        # inner constructor method
        new{
            typeof(Kn),
            typeof(Ma),
            typeof(Pr),
            typeof(K),
            typeof(γ),
            typeof(mi),
            typeof(ni),
            typeof(me),
            typeof(ne),
            typeof(lD),
            typeof(rL),
            typeof(sol),
            typeof(χ),
            typeof(ν),
            typeof(A1p),
            typeof(D1),
        }(
            Kn,
            Ma,
            Pr,
            K,
            γ,
            mi,
            ni,
            me,
            ne,
            lD,
            rL,
            sol,
            χ,
            ν,
            A1p,
            A1n,
            A2p,
            A2n,
            D1,
            D2,
        )
    end

end
