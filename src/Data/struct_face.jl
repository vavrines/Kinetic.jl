# ============================================================
# Structs of Cell Interface
# compatible with control volume structs
# ============================================================

# ------------------------------------------------------------
# 1D
# ------------------------------------------------------------

"""
1D cell interface with no distribution function

    @vars: fw

`Interface1D(w::AbstractArray)`

"""
mutable struct Interface1D{A} <: AbstractInterface1D

    fw::A

    function Interface1D(w::AbstractArray)
        fw = zeros(eltype(w), axes(w))
        new{typeof(fw)}(fw)
    end

end


"""
1D cell interface with 1 distribution function

    @vars: fw, ff

`Interface1D1F(w::AbstractArray, f::AbstractArray)`

"""
mutable struct Interface1D1F{A,B} <: AbstractInterface1D

    fw::A
    ff::B

    function Interface1D1F(w::AbstractArray, f::AbstractArray)
        fw = zeros(eltype(w), axes(w))
        ff = zeros(eltype(f), axes(f))

        new{typeof(fw),typeof(ff)}(fw, ff)
    end

end


"""
1D cell interface with 2 distribution functions

    @vars: fw, fh, fb

`Interface1D2F(w::AbstractArray, f::AbstractArray)`

"""
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


"""
1D cell interface with 3 distribution functions

    @vars: fw, fh0, fh1, fh2, femL, femR,

- deterministic: `Interface1D3F(w::AbstractArray, f::AbstractArray, E::AbstractArray{<:Real,1})`
- stochastic: `Interface1D3F(w::AbstractArray, f::AbstractArray, E::AbstractArray{<:Real,2})`

"""
mutable struct Interface1D3F{A,B,C} <: AbstractInterface1D

    fw::A
    fh0::B
    fh1::B
    fh2::B
    femL::C
    femR::C

    function Interface1D3F(
        w::AbstractArray,
        f::AbstractArray,
        E::AbstractArray{<:Real,1},
    )
        fw = zeros(eltype(w), axes(w))
        fh0 = zeros(eltype(f), axes(f))
        fh1 = zeros(eltype(f), axes(f))
        fh2 = zeros(eltype(f), axes(f))
        femL = zeros(eltype(E), 8)
        femR = zeros(eltype(E), 8)

        new{typeof(fw),typeof(fh0),typeof(femL)}(fw, fh0, fh1, fh2, femL, femR)
    end

    function Interface1D3F(
        w::AbstractArray,
        f::AbstractArray,
        E::AbstractArray{<:Real,2},
    )
        fw = zeros(eltype(w), axes(w))
        fh0 = zeros(eltype(f), axes(f))
        fh1 = zeros(eltype(f), axes(f))
        fh2 = zeros(eltype(f), axes(f))
        femL = zeros(eltype(E), 8, axes(E, 2))
        femR = zeros(eltype(E), 8, axes(E, 2))

        new{typeof(fw),typeof(fh0),typeof(femL)}(fw, fh0, fh1, fh2, femL, femR)
    end

end


"""
1D cell interface with 4 distribution functions

    @vars: fw, fh0, fh1, fh2, fh3, femL, femR,

- deterministic: `Interface1D4F(w::AbstractArray, f::AbstractArray, E::AbstractArray{<:Real,1})`
- stochastic: `Interface1D4F(w::AbstractArray, f::AbstractArray, E::AbstractArray{<:Real,2})`

"""
mutable struct Interface1D4F{A,B,C} <: AbstractInterface1D

    fw::A
    fh0::B
    fh1::B
    fh2::B
    fh3::B
    femL::C
    femR::C

    function Interface1D4F(
        w::AbstractArray,
        f::AbstractArray,
        E::AbstractArray{<:Real,1},
    )
        fw = zeros(eltype(w), axes(w))
        fh0 = zeros(eltype(f), axes(f))
        fh1 = zeros(eltype(f), axes(f))
        fh2 = zeros(eltype(f), axes(f))
        fh3 = zeros(eltype(f), axes(f))
        femL = zeros(eltype(E), 8)
        femR = zeros(eltype(E), 8)

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

    function Interface1D4F(
        w::AbstractArray,
        f::AbstractArray,
        E::AbstractArray{<:Real,2},
    )
        fw = zeros(eltype(w), axes(w))
        fh0 = zeros(eltype(f), axes(f))
        fh1 = zeros(eltype(f), axes(f))
        fh2 = zeros(eltype(f), axes(f))
        fh3 = zeros(eltype(f), axes(f))
        femL = zeros(eltype(E), 8, axes(E, 2))
        femR = zeros(eltype(E), 8, axes(E, 2))

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

# ------------------------------------------------------------
# 2D
# ------------------------------------------------------------

"""
2D cell interface with no distribution function

    @vars: len, n, fw

`Interface2D(L::Real, C::Real, S::Real, w::AbstractArray)`

"""
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


"""
2D cell interface with 1 distribution function

    @vars: len, n, fw, ff

`Interface2D1F(L::Real, C::Real, S::Real, w::AbstractArray, f::AbstractArray)`

"""
mutable struct Interface2D1F{A,B,C,D} <: AbstractInterface2D

    len::A
    n::B

    fw::C
    ff::D

    function Interface2D1F(
        L::Real,
        C::Real,
        S::Real,
        w::AbstractArray,
        f::AbstractArray,
    )
        len = L
        n = [C, S]

        fw = zeros(eltype(w), axes(w))
        ff = zeros(eltype(f), axes(f))

        new{typeof(len),typeof(n),typeof(fw),typeof(ff)}(len, n, fw, ff)
    end

end


"""
2D cell interface with 2 distribution functions

    @vars: len, n, fw, fh, fb

`Interface2D2F(L::Real, C::Real, S::Real, w::AbstractArray, f::AbstractArray)`

"""
mutable struct Interface2D2F{A,B,C,D} <: AbstractInterface2D

    len::A
    n::B

    fw::C
    fh::D
    fb::D

    function Interface2D2F(
        L::Real,
        C::Real,
        S::Real,
        w::AbstractArray,
        f::AbstractArray,
    )
        len = L
        n = [C, S]

        fw = zeros(eltype(w), axes(w))
        fh = zeros(eltype(f), axes(f))
        fb = zeros(eltype(f), axes(f))

        new{typeof(len),typeof(n),typeof(fw),typeof(fh)}(len, n, fw, fh, fb)
    end

end
