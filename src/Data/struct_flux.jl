# ============================================================
# Structs of Interface Fluxes
# compatible with solution structs
# ============================================================


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