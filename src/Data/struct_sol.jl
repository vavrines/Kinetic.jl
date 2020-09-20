# ============================================================
# Structs of Solution
# -> array of arrays in numerical simulation
# ============================================================


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