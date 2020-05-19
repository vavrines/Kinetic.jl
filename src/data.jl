# ============================================================
# Data Structures
# ============================================================


export Setup
export GasProperty, PlasmaProperty
export IB1D1F, IB1D2F, IB1D4F
export ControlVolume1D1F, ControlVolume1D2F, ControlVolume1D4F
export Interface1D1F, Interface1D2F, Interface1D4F
export Solution1D1F, Solution1D2F


# ------------------------------------------------------------
# Structure of computational setup
# ------------------------------------------------------------
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


# ------------------------------------------------------------
# Structure of property
# ------------------------------------------------------------
struct GasProperty{A,B,C,D,E,F,G,H,I} <: AbstractProperty

    Kn::A
    Ma::B
    Pr::C
    K::D
    γ::E
    ω::F
    αᵣ::G
    ωᵣ::H
    μᵣ::I

    function GasProperty(
        Kn::Real,
        Ma::Real,
        Pr::Real,
        K::Real,
        γ::Real,
        ω::Real,
        αᵣ::Real,
        ωᵣ::Real,
        μᵣ::Real,
    )

        # inner constructor method
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


struct PlasmaProperty{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P} <: AbstractProperty

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
    D1::P

    function PlasmaProperty(
        Kn::Array{<:Real,1},
        Ma::Real,
        Pr::Real,
        K::Real,
        γ::Real,
        mi::Real,
        ni::Real,
        me::Real,
        ne::Real,
        lD::Real,
        rL::Real,
        sol::Real,
        χ::Real,
        ν::Real,
    )

        # A^+
        A1p = Array{typeof(sol)}(undef, 8, 8)
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

        # A^-
        A1n = Array{typeof(sol)}(undef, 8, 8)
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
            D1,
        )

    end

end


# ------------------------------------------------------------
# Structure of initial and boundary conditions
# ------------------------------------------------------------
struct IB1D1F{A,B} <: AbstractCondition

    wL::A
    primL::A
    fL::B
    bcL::A

    wR::A
    primR::A
    fR::B
    bcR::A

    # works for both 1V/3V and single-/multi-component gases
    function IB1D1F(
        wL::Array,
        primL::Array,
        fL::AbstractArray,
        bcL::Array,
        wR::Array,
        primR::Array,
        fR::AbstractArray,
        bcR::Array,
    )

        new{typeof(wL),typeof(fL)}(wL, primL, fL, bcL, wR, primR, fR, bcR)

    end

end

struct IB1D2F{A,B} <: AbstractCondition

    # initial condition
    wL::A
    primL::A
    hL::B
    bL::B
    bcL::A

    wR::A
    primR::A
    hR::B
    bR::B
    bcR::A

    function IB1D2F(
        wL::Array,
        primL::Array,
        hL::AbstractArray,
        bL::AbstractArray,
        bcL::Array,
        wR::Array,
        primR::Array,
        hR::AbstractArray,
        bR::AbstractArray,
        bcR::Array,
    )

        new{typeof(wL),typeof(hL)}(
            wL,
            primL,
            hL,
            bL,
            bcL,
            wR,
            primR,
            hR,
            bR,
            bcR,
        )

    end

end


struct IB1D4F{A,B,C,D} <: AbstractCondition

    # initial/boundary condition
    wL::A
    primL::A
    h0L::B
    h1L::B
    h2L::B
    h3L::B
    bcL::A
    EL::C
    BL::C
    lorenzL::D

    wR::A
    primR::A
    h0R::B
    h1R::B
    h2R::B
    h3R::B
    bcR::A
    ER::C
    BR::C
    lorenzR::D

    function IB1D4F(
        wL::Array,
        primL::Array,
        h0L::AbstractArray,
        h1L::AbstractArray,
        h2L::AbstractArray,
        h3L::AbstractArray,
        bcL::Array,
        EL::Array,
        BL::Array,
        lorenzL::Array,
        wR::Array,
        primR::Array,
        h0R::AbstractArray,
        h1R::AbstractArray,
        h2R::AbstractArray,
        h3R::AbstractArray,
        bcR::Array,
        ER::Array,
        BR::Array,
        lorenzR::Array,
    )

        new{typeof(wL),typeof(h0L),typeof(EL),typeof(lorenzL)}(
            wL,
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
            lorenzR,
        )

    end

end


# ------------------------------------------------------------
# Structure of control volume
# ------------------------------------------------------------
mutable struct ControlVolume1D1F{F,A,B} <: AbstractControlVolume1D

    x::F
    dx::F

    w::A
    prim::A
    sw::A

    f::B
    sf::B

    function ControlVolume1D1F(
        X::Real,
        DX::Real,
        W::Array,
        PRIM::Array,
        F::AbstractArray,
    )

        x = deepcopy(X)
        dx = deepcopy(DX)

        w = deepcopy(W)
        prim = deepcopy(PRIM)
        sw = zeros(typeof(W[1]), axes(w))

        f = deepcopy(F)
        sf = zeros(typeof(F[1]), axes(f))

        new{typeof(x),typeof(w),typeof(f)}(x, dx, w, prim, sw, f, sf)

    end

end


mutable struct ControlVolume1D2F{F,A,B} <: AbstractControlVolume1D

    x::F
    dx::F

    w::A
    prim::A
    sw::A

    h::B
    b::B
    sh::B
    sb::B

    function ControlVolume1D2F(
        X::Real,
        DX::Real,
        W::Array,
        PRIM::Array,
        H::AbstractArray,
        B::AbstractArray,
    )

        x = deepcopy(X)
        dx = deepcopy(DX)

        w = deepcopy(W)
        prim = deepcopy(PRIM)
        sw = zeros(typeof(W[1]), axes(W))

        h = deepcopy(H)
        b = deepcopy(B)
        sh = zeros(typeof(H[1]), axes(H))
        sb = zeros(typeof(B[1]), axes(B))

        new{typeof(x),typeof(w),typeof(h)}(x, dx, w, prim, sw, h, b, sh, sb)

    end

end


mutable struct ControlVolume1D4F{F,A,B,C,D,E} <: AbstractControlVolume1D

    x::F
    dx::F

    w::A
    prim::A
    sw::A

    h0::B
    h1::B
    h2::B
    h3::B
    sh0::B
    sh1::B
    sh2::B
    sh3::B

    E::C
    B::C
    ϕ::D
    ψ::D
    lorenz::E

    # deterministic
    function ControlVolume1D4F(
        X::Real,
        DX::Real,
        W::Array{<:Real,2},
        PRIM::Array{<:Real,2},
        H0::AbstractArray{<:AbstractFloat,2},
        H1::AbstractArray{Float64,2},
        H2::AbstractArray{Float64,2},
        H3::AbstractArray{Float64,2},
        E0::Array{Float64,1},
        B0::Array{Float64,1},
        L::Array{Float64,2},
    )

        x = deepcopy(X)
        dx = deepcopy(DX)

        w = deepcopy(W)
        prim = deepcopy(PRIM)
        sw = zeros(typeof(W[1]), axes(W))

        h0 = deepcopy(H0)
        h1 = deepcopy(H1)
        h2 = deepcopy(H2)
        h3 = deepcopy(H3)
        sh0 = zeros(typeof(H0[1]), axes(H0))
        sh1 = zeros(typeof(H1[1]), axes(H1))
        sh2 = zeros(typeof(H2[1]), axes(H2))
        sh3 = zeros(typeof(H3[1]), axes(H3))

        E = deepcopy(E0)
        B = deepcopy(B0)
        ϕ = 0.0
        ψ = 0.0
        lorenz = deepcopy(L)

        new{typeof(x),typeof(w),typeof(h0),typeof(E),typeof(ϕ),typeof(lorenz)}(
            x,
            dx,
            w,
            prim,
            sw,
            h0,
            h1,
            h2,
            h3,
            sh0,
            sh1,
            sh2,
            sh3,
            E,
            B,
            ϕ,
            ψ,
            lorenz,
        )

    end

    # uncertainty quantification
    function ControlVolume1D4F(
        X::Real,
        DX::Real,
        W::Array{<:Real,3},
        PRIM::Array{<:Real,3},
        H0::AbstractArray{<:AbstractFloat,3},
        H1::AbstractArray{Float64,3},
        H2::AbstractArray{Float64,3},
        H3::AbstractArray{Float64,3},
        E0::Array{Float64,2},
        B0::Array{Float64,2},
        L::Array{Float64,3},
    )

        x = deepcopy(X)
        dx = deepcopy(DX)

        w = deepcopy(W)
        prim = deepcopy(PRIM)
        sw = zeros(typeof(W[1]), axes(W))

        h0 = deepcopy(H0)
        h1 = deepcopy(H1)
        h2 = deepcopy(H2)
        h3 = deepcopy(H3)
        sh0 = zeros(typeof(H0[1]), axes(H0))
        sh1 = zeros(typeof(H1[1]), axes(H1))
        sh2 = zeros(typeof(H2[1]), axes(H2))
        sh3 = zeros(typeof(H3[1]), axes(H3))

        E = deepcopy(E0)
        B = deepcopy(B0)
        ϕ = zeros(axes(E, 2))
        ψ = zeros(axes(B, 2))
        lorenz = deepcopy(L)

        new{typeof(x),typeof(w),typeof(h0),typeof(E),typeof(ϕ),typeof(lorenz)}(
            x,
            dx,
            w,
            prim,
            sw,
            h0,
            h1,
            h2,
            h3,
            sh0,
            sh1,
            sh2,
            sh3,
            E,
            B,
            ϕ,
            ψ,
            lorenz,
        )

    end

end


# ------------------------------------------------------------
# Structure of cell interface
# ------------------------------------------------------------
mutable struct Interface1D1F{A,B} <: AbstractInterface1D

    fw::A
    ff::B

    function Interface1D1F(w::Array, f::AbstractArray)

        fw = zeros(typeof(w[1]), axes(w))
        ff = zeros(typeof(f[1]), axes(f))

        new{typeof(fw),typeof(ff)}(fw, ff)

    end

end


mutable struct Interface1D2F{A,B} <: AbstractInterface1D

    fw::A
    fh::B
    fb::B

    function Interface1D2F(w::Array, f::AbstractArray)

        fw = zeros(typeof(w[1]), axes(w))
        fh = zeros(typeof(f[1]), axes(f))
        fb = zeros(typeof(f[1]), axes(f))

        new{typeof(fw),typeof(fh)}(fw, fh, fb)

    end

end


mutable struct Interface1D4F{A,B,C} <: AbstractInterface1D

    fw::A
    fh0::B
    fh1::B
    fh2::B
    fh3::B
    femL::C
    femR::C

    function Interface1D4F(w::Array, f::AbstractArray, E::Array{<:Real,1})

        fw = zeros(typeof(w[1]), axes(w))
        fh0 = zeros(typeof(f[1]), axes(f))
        fh1 = zeros(typeof(f[1]), axes(f))
        fh2 = zeros(typeof(f[1]), axes(f))
        fh3 = zeros(typeof(f[1]), axes(f))
        femL = zeros(typeof(E), 8)
        femR = zeros(typeof(E), 8)

        new{typeof(fw),typeof(fh0),typeof(femL)}(
            fw,
            fh0,
            fh1,
            fh2,
            fh3,
            femL,
            femR,
        )

    end

    function Interface1D4F(w::Array, f::AbstractArray, E::Array{<:Real,2})

        fw = zeros(typeof(w[1]), axes(w))
        fh0 = zeros(typeof(f[1]), axes(f))
        fh1 = zeros(typeof(f[1]), axes(f))
        fh2 = zeros(typeof(f[1]), axes(f))
        fh3 = zeros(typeof(f[1]), axes(f))
        femL = zeros(typeof(E), 8, axes(E, 2))
        femR = zeros(typeof(E), 8, axes(E, 2))

        new{typeof(fw),typeof(fh0),typeof(femL)}(
            fw,
            fh0,
            fh1,
            fh2,
            fh3,
            femL,
            femR,
        )

    end

end


mutable struct Solution1D1F{A,B} <: AbstractSolution

    w::A
    prim::A
    sw::A
    f::B
    sf::B

    function Solution1D1F(w::Array, prim::Array, f::AbstractArray)
        sw = zeros(typeof(w[1]), axes(w))
        sf = zeros(typeof(f[1]), axes(f))

        new{typeof(w),typeof(f)}(w, prim, sw, f, sf)
    end

    function Solution1D1F(
        w::Array,
        prim::Array,
        sw::Array,
        f::AbstractArray,
        sf::AbstractArray,
    )
        new{typeof(w),typeof(f)}(w, prim, sw, f, sf)
    end

end


mutable struct Solution1D2F{A,B} <: AbstractSolution

    w::A
    prim::A
    sw::A
    h::B
    b::B
    sh::B
    sb::B

    function Solution1D2F(
        w::Array,
        prim::Array,
        h::AbstractArray,
        b::AbstractArray,
    )
        sw = zeros(typeof(w[1]), axes(w))
        sh = zeros(typeof(h[1]), axes(h))
        sb = zeros(typeof(b[1]), axes(b))

        new{typeof(w),typeof(h)}(w, prim, sw, h, b, sh, sb)
    end

    function Solution1D2F(
        w::Array,
        prim::Array,
        sw::Array,
        h::AbstractArray,
        b::AbstractArray,
        sh::AbstractArray,
        sb::AbstractArray,
    )
        new{typeof(w),typeof(h)}(w, prim, sw, h, b, sh, sb)
    end

end
