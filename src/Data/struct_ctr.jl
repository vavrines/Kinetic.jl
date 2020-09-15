# ============================================================
# Structs of Control Volume
# constituting array of struct in numerical simulation
# ============================================================


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


mutable struct ControlVolume1D3F{F,A,B,C,D,E} <: AbstractControlVolume1D

    x::F
    dx::F

    w::A
    prim::A
    sw::A

    h0::B
    h1::B
    h2::B
    sh0::B
    sh1::B
    sh2::B

    E::C
    B::C
    ϕ::D
    ψ::D
    lorenz::E

    # deterministic
    function ControlVolume1D3F(
        X::Real,
        DX::Real,
        W::AbstractArray{<:Real,2},
        PRIM::AbstractArray{<:Real,2},
        H0::AbstractArray{<:AbstractFloat,3},
        H1::AbstractArray{<:AbstractFloat,3},
        H2::AbstractArray{<:AbstractFloat,3},
        E0::AbstractArray{<:AbstractFloat,1},
        B0::AbstractArray{<:AbstractFloat,1},
        L::AbstractArray{<:AbstractFloat,2},
    )
        x = deepcopy(X)
        dx = deepcopy(DX)

        w = deepcopy(W)
        prim = deepcopy(PRIM)
        sw = zeros(eltype(H0), axes(W))

        h0 = deepcopy(H0)
        h1 = deepcopy(H1)
        h2 = deepcopy(H2)
        sh0 = zeros(eltype(H0), axes(H0))
        sh1 = zeros(eltype(H1), axes(H1))
        sh2 = zeros(eltype(H2), axes(H2))

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

    # stochastic
    function ControlVolume1D3F(
        X::Real,
        DX::Real,
        W::AbstractArray{<:Real,3},
        PRIM::AbstractArray{<:Real,3},
        H0::AbstractArray{<:AbstractFloat,4},
        H1::AbstractArray{<:AbstractFloat,4},
        H2::AbstractArray{<:AbstractFloat,4},
        E0::AbstractArray{<:AbstractFloat,2},
        B0::AbstractArray{<:AbstractFloat,2},
        L::AbstractArray{<:AbstractFloat,3},
    )
        x = deepcopy(X)
        dx = deepcopy(DX)

        w = deepcopy(W)
        prim = deepcopy(PRIM)
        sw = zeros(eltype(H0), axes(W))

        h0 = deepcopy(H0)
        h1 = deepcopy(H1)
        h2 = deepcopy(H2)
        sh0 = zeros(eltype(H0), axes(H0))
        sh1 = zeros(eltype(H1), axes(H1))
        sh2 = zeros(eltype(H2), axes(H2))

        E = deepcopy(E0)
        B = deepcopy(B0)
        ϕ = zeros(axes(E, 2)) # here is difference
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

        new{
            typeof(x),
            typeof(w),
            typeof(sw),
            typeof(h0),
            typeof(sh2),
            typeof(E),
            typeof(ϕ),
            typeof(lorenz),
        }(
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