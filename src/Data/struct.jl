# ============================================================
# Structs
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
Structure of property

"""
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

    # unified consideration of deterministic and stochastic conditions
    # {mi, ni, me, ne, sol, χ, ν} keep the same as before
    function PlasmaProperty(
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
        A1p = Array{Float64}(undef, 8, 8)
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
        A1n = Array{Float64}(undef, 8, 8)
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
struct IB{A} <: AbstractCondition

    wL::A
    primL::A
    bcL::A

    wR::A
    primR::A
    bcR::A

    bcU::A
    bcD::A

    # works for both 1V/3V and single-/multi-component gases
    function IB(
        wL::AbstractArray,
        primL::AbstractArray,
        bcL::AbstractArray,
        wR::AbstractArray,
        primR::AbstractArray,
        bcR::AbstractArray,
        bcU = deepcopy(bcR)::AbstractArray,
        bcD = deepcopy(bcR)::AbstractArray,
    )

        new{typeof(wL)}(wL, primL, bcL, wR, primR, bcR, bcU, bcD)

    end

end


struct IB1F{A,B} <: AbstractCondition

    wL::A
    primL::A
    fL::B
    bcL::A

    wR::A
    primR::A
    fR::B
    bcR::A

    bcU::A
    bcD::A

    # works for both 1V/3V and single-/multi-component gases
    function IB1F(
        wL::AbstractArray,
        primL::AbstractArray,
        fL::AbstractArray,
        bcL::AbstractArray,
        wR::AbstractArray,
        primR::AbstractArray,
        fR::AbstractArray,
        bcR::Array,
        bcU = deepcopy(bcR)::AbstractArray,
        bcD = deepcopy(bcR)::AbstractArray,
    )

        new{typeof(wL),typeof(fL)}(wL, primL, fL, bcL, wR, primR, fR, bcR, bcU, bcD)

    end

end


struct IB2F{A,B} <: AbstractCondition

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

    bcU::A
    bcD::A

    function IB2F(
        wL::AbstractArray,
        primL::AbstractArray,
        hL::AbstractArray,
        bL::AbstractArray,
        bcL::AbstractArray,
        wR::AbstractArray,
        primR::AbstractArray,
        hR::AbstractArray,
        bR::AbstractArray,
        bcR::AbstractArray,
        bcU = deepcopy(bcR)::AbstractArray,
        bcD = deepcopy(bcR)::AbstractArray,
    )

        new{typeof(wL),typeof(hL)}(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR, bcU, bcD)

    end

end


struct IB4F{A,B,C,D} <: AbstractCondition

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

    bcU::A
    bcD::A

    function IB4F(
        wL::AbstractArray,
        primL::AbstractArray,
        h0L::AbstractArray,
        h1L::AbstractArray,
        h2L::AbstractArray,
        h3L::AbstractArray,
        bcL::AbstractArray,
        EL::AbstractArray,
        BL::AbstractArray,
        lorenzL::AbstractArray,
        wR::AbstractArray,
        primR::AbstractArray,
        h0R::AbstractArray,
        h1R::AbstractArray,
        h2R::AbstractArray,
        h3R::AbstractArray,
        bcR::AbstractArray,
        ER::AbstractArray,
        BR::AbstractArray,
        lorenzR::AbstractArray,
        bcU = deepcopy(bcR)::AbstractArray,
        bcD = deepcopy(bcR)::AbstractArray,
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
            bcU,
            bcD,
        )

    end

end


"""
Structure of control volume -> array of struct in flow simulation

"""
mutable struct ControlVolume1D{F,A} <: AbstractControlVolume1D

    x::F
    dx::F

    w::A
    prim::A
    sw::A

    function ControlVolume1D(X::Real, DX::Real, W::AbstractArray, PRIM::AbstractArray)

        x = deepcopy(X)
        dx = deepcopy(DX)

        w = deepcopy(W)
        prim = deepcopy(PRIM)
        sw = zeros(typeof(W[1]), axes(w))

        new{typeof(x),typeof(w)}(x, dx, w, prim, sw)

    end

end


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
        W::AbstractArray,
        PRIM::AbstractArray,
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
        W::AbstractArray,
        PRIM::AbstractArray,
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
        W::AbstractArray{<:Real,2},
        PRIM::AbstractArray{<:Real,2},
        H0::AbstractArray{<:AbstractFloat,2},
        H1::AbstractArray{Float64,2},
        H2::AbstractArray{Float64,2},
        H3::AbstractArray{Float64,2},
        E0::AbstractArray{Float64,1},
        B0::AbstractArray{Float64,1},
        L::AbstractArray{Float64,2},
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
        W::AbstractArray{<:Real,3},
        PRIM::AbstractArray{<:Real,3},
        H0::AbstractArray{<:AbstractFloat,3},
        H1::AbstractArray{Float64,3},
        H2::AbstractArray{Float64,3},
        H3::AbstractArray{Float64,3},
        E0::AbstractArray{Float64,2},
        B0::AbstractArray{Float64,2},
        L::AbstractArray{Float64,3},
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


mutable struct ControlVolume2D{F,A,B} <: AbstractControlVolume2D

    x::F
    y::F
    dx::F
    dy::F

    w::A
    prim::A
    sw::B

    function ControlVolume2D(
        X::Real,
        DX::Real,
        Y::Real,
        DY::Real,
        W::AbstractArray,
        PRIM::AbstractArray,
    )

        x = deepcopy(X)
        dx = deepcopy(DX)
        y = deepcopy(Y)
        dy = deepcopy(DY)

        w = deepcopy(W)
        prim = deepcopy(PRIM)
        sw = zeros(eltype(W), (axes(W)..., Base.OneTo(2)))

        new{typeof(x),typeof(w),typeof(sw)}(x, dx, y, dy, w, prim, sw)

    end

end


mutable struct ControlVolume2D1F{F,A,B,C,D} <: AbstractControlVolume2D

    x::F
    y::F
    dx::F
    dy::F

    w::A
    prim::A
    sw::B

    f::C
    sf::D

    function ControlVolume2D1F(
        X::Real,
        DX::Real,
        Y::Real,
        DY::Real,
        W::AbstractArray,
        PRIM::AbstractArray,
        F::AbstractArray,
    )

        x = deepcopy(X)
        dx = deepcopy(DX)
        y = deepcopy(Y)
        dy = deepcopy(DY)

        w = deepcopy(W)
        prim = deepcopy(PRIM)
        sw = zeros(eltype(W), (axes(W)..., Base.OneTo(2)))

        f = deepcopy(F)
        sf = zeros(eltype(F), (axes(F)..., Base.OneTo(2)))

        new{typeof(x),typeof(w),typeof(sw),typeof(f),typeof(sf)}(
            x,
            dx,
            y,
            dy,
            w,
            prim,
            sw,
            f,
            sf,
        )

    end

end


mutable struct ControlVolume2D2F{F,A,B,C,D} <: AbstractControlVolume2D

    x::F
    y::F
    dx::F
    dy::F

    w::A
    prim::A
    sw::B

    h::C
    b::C
    sh::D
    sb::D

    function ControlVolume2D2F(
        X::Real,
        DX::Real,
        Y::Real,
        DY::Real,
        W::AbstractArray,
        PRIM::AbstractArray,
        H::AbstractArray,
        B::AbstractArray,
    )

        x = deepcopy(X)
        dx = deepcopy(DX)
        y = deepcopy(Y)
        dy = deepcopy(DY)

        w = deepcopy(W)
        prim = deepcopy(PRIM)
        sw = zeros(eltype(W), (axes(W)..., Base.OneTo(2)))

        h = deepcopy(H)
        b = deepcopy(B)
        sh = zeros(eltype(H), (axes(H)..., Base.OneTo(2)))
        sb = zeros(eltype(B), (axes(B)..., Base.OneTo(2)))

        new{typeof(x),typeof(w),typeof(sw),typeof(h),typeof(sh)}(
            x,
            dx,
            y,
            dy,
            w,
            prim,
            sw,
            h,
            b,
            sh,
            sb,
        )

    end

end


mutable struct ControlVolume2D3F{F,A,B,C,D,E,F,G} <: AbstractControlVolume2D

    x::F
    y::F
    dx::F
    dy::F

    w::A
    prim::A
    sw::B

    h0::C
    h1::C
    h2::C
    sh0::D
    sh1::D
    sh2::D

    E::E
    B::E
    ϕ::F
    ψ::F
    lorenz::G

    # deterministic & stochastic
    function ControlVolume2D3F(
        X::Real,
        DX::Real,
        Y::Real,
        DY::Real,
        W::AbstractArray,
        PRIM::AbstractArray,
        H0::AbstractArray,
        H1::AbstractArray,
        H2::AbstractArray,
        E0::AbstractArray,
        B0::AbstractArray,
        L::AbstractArray,
    )
        x = deepcopy(X)
        dx = deepcopy(DX)
        y = deepcopy(Y)
        dy = deepcopy(DY)

        w = deepcopy(W)
        prim = deepcopy(PRIM)
        sw = zeros(eltype(W), (axes(W)..., Base.OneTo(2))) # 2D

        h0 = deepcopy(H0)
        h1 = deepcopy(H1)
        h2 = deepcopy(H2)
        sh0 = zeros(eltype(H0), (axes(H0)..., Base.OneTo(2)))
        sh1 = zeros(eltype(H1), (axes(H1)..., Base.OneTo(2)))
        sh2 = zeros(eltype(H2), (axes(H2)..., Base.OneTo(2)))

        E = deepcopy(E0)
        B = deepcopy(B0)
        ϕ = 0.0
        ψ = 0.0
        lorenz = deepcopy(L)

        new{typeof(x),typeof(w),typeof(sw),typeof(h0),typeof(sh2),typeof(E),typeof(ϕ),typeof(lorenz)}(
            x,
            dx,
            y,
            dy,
            w,
            prim,
            sw,
            h0,
            h1,
            h2,
            sh0,
            sh1,
            sh2,
            E,
            B,
            ϕ,
            ψ,
            lorenz,
        )
    end

end


mutable struct Interface1D{A} <: AbstractInterface1D

    fw::A

    function Interface1D(w::AbstractArray)

        fw = zeros(eltype(w), axes(w))

        new{typeof(fw)}(fw)

    end

end


mutable struct Interface1D1F{A,B} <: AbstractInterface1D

    fw::A
    ff::B

    function Interface1D1F(w::AbstractArray, f::AbstractArray)

        fw = zeros(eltype(w), axes(w))
        ff = zeros(eltype(f), axes(f))

        new{typeof(fw),typeof(ff)}(fw, ff)

    end

end


mutable struct Interface1D2F{A,B} <: AbstractInterface1D

    fw::A
    fh::B
    fb::B

    function Interface1D2F(w::AbstractArray, f::AbstractArray)

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

    function Interface1D4F(w::AbstractArray, f::AbstractArray, E::AbstractArray{<:Real,1})
        fw = zeros(eltype(w), axes(w))
        fh0 = zeros(eltype(f), axes(f))
        fh1 = zeros(eltype(f), axes(f))
        fh2 = zeros(eltype(f), axes(f))
        fh3 = zeros(eltype(f), axes(f))
        femL = zeros(eltype(E), 8)
        femR = zeros(eltype(E), 8)

        new{typeof(fw),typeof(fh0),typeof(femL)}(fw, fh0, fh1, fh2, fh3, femL, femR)
    end

    function Interface1D4F(w::AbstractArray, f::AbstractArray, E::AbstractArray{<:Real,2})
        fw = zeros(eltype(w), axes(w))
        fh0 = zeros(eltype(f), axes(f))
        fh1 = zeros(eltype(f), axes(f))
        fh2 = zeros(eltype(f), axes(f))
        fh3 = zeros(eltype(f), axes(f))
        femL = zeros(eltype(E), 8, axes(E, 2))
        femR = zeros(eltype(E), 8, axes(E, 2))

        new{typeof(fw),typeof(fh0),typeof(femL)}(fw, fh0, fh1, fh2, fh3, femL, femR)
    end

end


mutable struct Interface2D{A,B,C} <: AbstractInterface2D

    len::A
    n::B
    fw::C

    function Interface2D(L::Real, C::Real, S::Real, w::AbstractArray)

        len = L
        n = [C, S]

        fw = zeros(eltype(w), axes(w))

        new{typeof(len),typeof(n),typeof(fw)}(len, n, fw)

    end

end


mutable struct Interface2D1F{A,B,C,D} <: AbstractInterface2D

    len::A
    n::B

    fw::C
    ff::D

    function Interface2D1F(L::Real, C::Real, S::Real, w::AbstractArray, f::AbstractArray)

        len = L
        n = [C, S]

        fw = zeros(eltype(w), axes(w))
        ff = zeros(eltype(f), axes(f))

        new{typeof(len),typeof(n),typeof(fw),typeof(ff)}(len, n, fw, ff)

    end

end


mutable struct Interface2D2F{A,B,C,D} <: AbstractInterface2D

    len::A
    n::B

    fw::C
    fh::D
    fb::D

    function Interface2D2F(L::Real, C::Real, S::Real, w::AbstractArray, f::AbstractArray)

        len = L
        n = [C, S]

        fw = zeros(eltype(w), axes(w))
        fh = zeros(eltype(f), axes(f))
        fb = zeros(eltype(f), axes(f))

        new{typeof(len),typeof(n),typeof(fw),typeof(fh)}(len, n, fw, fh, fb)

    end

end


"""
Data structure 2
solution & flux : array of arrays

"""
mutable struct Solution1D{A} <: AbstractSolution1D

    w::A
    prim::A
    sw::A

    function Solution1D(
        w::AbstractArray,
        prim::AbstractArray,
        sw = [zeros(axes(w[1])) for i in axes(w, 1)]::AbstractArray,
    )
        new{typeof(w)}(w, prim, sw)
    end

end


mutable struct Solution1D1F{A,B} <: AbstractSolution1D

    w::A
    prim::A
    sw::A
    f::B
    sf::B

    function Solution1D1F(w::AbstractArray, prim::AbstractArray, f::AbstractArray)
        sw = [zeros(axes(w[1])) for i in axes(w, 1)]
        sf = [zeros(axes(f[1])) for i in axes(f, 1)]

        new{typeof(w),typeof(f)}(w, prim, sw, f, sf)
    end

    function Solution1D1F(
        w::AbstractArray,
        prim::AbstractArray,
        sw::AbstractArray,
        f::AbstractArray,
        sf::AbstractArray,
    )
        new{typeof(w),typeof(f)}(w, prim, sw, f, sf)
    end

end


mutable struct Solution1D2F{A,B} <: AbstractSolution1D

    w::A
    prim::A
    sw::A
    h::B
    b::B
    sh::B
    sb::B

    function Solution1D2F(
        w::AbstractArray,
        prim::AbstractArray,
        h::AbstractArray,
        b::AbstractArray,
    )
        sw = [zeros(axes(w[1])) for i in axes(w, 1)]
        sh = [zeros(axes(h[1])) for i in axes(h, 1)]
        sb = [zeros(axes(b[1])) for i in axes(b, 1)]

        new{typeof(w),typeof(h)}(w, prim, sw, h, b, sh, sb)
    end

    function Solution1D2F(
        w::AbstractArray,
        prim::AbstractArray,
        sw::AbstractArray,
        h::AbstractArray,
        b::AbstractArray,
        sh::AbstractArray,
        sb::AbstractArray,
    )
        new{typeof(w),typeof(h)}(w, prim, sw, h, b, sh, sb)
    end

end


mutable struct Solution2D{A,B} <: AbstractSolution2D

    w::A
    prim::A
    sw::B

    function Solution2D(
        w::AbstractArray,
        prim::AbstractArray,
        sw = [
            zeros((axes(w[1])..., Base.OneTo(2))) for i in axes(w, 1), j in axes(w, 2)
        ]::AbstractArray,
    )

        new{typeof(w),typeof(sw)}(w, prim, sw)
    end

end


mutable struct Solution2D1F{A,B,C,D} <: AbstractSolution2D

    w::A
    prim::A
    sw::B
    f::C
    sf::D

    function Solution2D1F(w::AbstractArray, prim::AbstractArray, f::AbstractArray)
        sw = [zeros((axes(w[1])..., Base.OneTo(2))) for i in axes(w, 1), j in axes(w, 2)]
        sf = [zeros((axes(f[1])..., Base.OneTo(2))) for i in axes(f, 1), j in axes(f, 2)]

        new{typeof(w),typeof(sw),typeof(f),typeof(sf)}(w, prim, sw, f, sf)
    end

    function Solution2D1F(
        w::AbstractArray,
        prim::AbstractArray,
        sw::AbstractArray,
        f::AbstractArray,
        sf::AbstractArray,
    )
        new{typeof(w),typeof(sw),typeof(f),typeof(sf)}(w, prim, sw, f, sf)
    end

end


mutable struct Solution2D2F{A,B,C,D} <: AbstractSolution2D

    w::A
    prim::A
    sw::B
    h::C
    b::C
    sh::D
    sb::D

    function Solution2D2F(
        w::AbstractArray,
        prim::AbstractArray,
        h::AbstractArray,
        b::AbstractArray,
    )
        sw = [zeros((axes(w[1])..., Base.OneTo(2))) for i in axes(w, 1), j in axes(w, 2)]
        sh = [zeros((axes(h[1])..., Base.OneTo(2))) for i in axes(h, 1), j in axes(h, 2)]
        sb = [zeros((axes(b[1])..., Base.OneTo(2))) for i in axes(b, 1), j in axes(b, 2)]

        new{typeof(w),typeof(sw),typeof(h),typeof(sh)}(w, prim, sw, h, b, sh, sb)
    end

    function Solution2D2F(
        w::AbstractArray,
        prim::AbstractArray,
        sw::AbstractArray,
        h::AbstractArray,
        b::AbstractArray,
        sh::AbstractArray,
        sb::AbstractArray,
    )
        new{typeof(w),typeof(sw),typeof(h),typeof(sh)}(w, prim, sw, h, b, sh, sb)
    end

end


mutable struct Flux1D{A,B} <: AbstractFlux1D

    w::A
    fw::B

    function Flux1D(w::AbstractArray, fw::AbstractArray)
        new{typeof(w),typeof(fw)}(w, fw)
    end

end


mutable struct Flux1D1F{A,B,C} <: AbstractFlux1D

    w::A
    fw::B
    ff::C

    function Flux1D1F(w::AbstractArray, fw::AbstractArray, ff::AbstractArray)
        new{typeof(w),typeof(fw),typeof(ff)}(w, fw, ff)
    end

end


mutable struct Flux1D2F{A,B,C} <: AbstractFlux1D

    w::A
    fw::B
    fh::C
    fb::C

    function Flux1D2F(
        w::AbstractArray,
        fw::AbstractArray,
        fh::AbstractArray,
        fb::AbstractArray,
    )
        new{typeof(w),typeof(fw),typeof(fh)}(w, fw, fh, fb)
    end

end


mutable struct Flux2D{A,B,C} <: AbstractFlux2D

    n1::A
    w1::B
    fw1::C

    n2::A
    w2::B
    fw2::C

    function Flux2D(
        n1::AbstractArray,
        w1::AbstractArray,
        fw1::AbstractArray,
        n2::AbstractArray,
        w2::AbstractArray,
        fw2::AbstractArray,
    )
        new{typeof(n1),typeof(w1),typeof(fw1)}(n1, w1, fw1, n2, w2, fw2)
    end

end


mutable struct Flux2D1F{A,B,C,D} <: AbstractFlux2D

    n1::A
    w1::B
    fw1::C
    ff1::D

    n2::A
    w2::B
    fw2::C
    ff2::D

    function Flux2D1F(
        n1::AbstractArray,
        w1::AbstractArray,
        fw1::AbstractArray,
        ff1::AbstractArray,
        n2::AbstractArray,
        w2::AbstractArray,
        fw2::AbstractArray,
        ff2::AbstractArray,
    )
        new{typeof(n1),typeof(w1),typeof(fw1),typeof(ff1)}(
            n1,
            w1,
            fw1,
            ff1,
            n2,
            w2,
            fw2,
            ff2,
        )
    end

end


mutable struct Flux2D2F{A,B,C,D} <: AbstractFlux2D

    n1::A
    w1::B
    fw1::C
    fh1::D
    fb1::D

    n2::A
    w2::B
    fw2::C
    fh2::D
    fb2::D

    function Flux2D2F(
        n1::AbstractArray,
        w1::AbstractArray,
        fw1::AbstractArray,
        fh1::AbstractArray,
        fb1::AbstractArray,
        n2::AbstractArray,
        w2::AbstractArray,
        fw2::AbstractArray,
        fh2::AbstractArray,
        fb2::AbstractArray,
    )
        new{typeof(n1),typeof(w1),typeof(fw1),typeof(fh1)}(
            n1,
            w1,
            fw1,
            fh1,
            fb1,
            n2,
            w2,
            fw2,
            fh2,
            fb2,
        )
    end

end
