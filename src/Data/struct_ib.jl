# ============================================================
# Structs of Initial and Boundary Conditions
# ============================================================

"""
Initial & boundary condition with no distribution function

    @consts: wL, primL, bcL, wR, primR, bcR, bcU, bcD

"""
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


"""
Initial & boundary condition with 1 distribution function

    @consts: wL, primL, fL, bcL, wR, primR, fR, bcR, bcU, bcD

"""
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


"""
Initial & boundary condition with 2 distribution functions

    @consts: wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR, bcU, bcD

"""
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


"""
Initial & boundary condition with 3 distribution functions

    @consts: wL, primL, h0L, h1L, h2L, bcL, EL, BL, lorenzL, wR, primR, h0R, h1R, h2R, bcR, ER, BR, lorenzR, bcU, bcD

"""
struct IB3F{A,B,C,D} <: AbstractCondition

    # initial/boundary condition
    wL::A
    primL::A
    h0L::B
    h1L::B
    h2L::B
    bcL::A
    EL::C
    BL::C
    lorenzL::D

    wR::A
    primR::A
    h0R::B
    h1R::B
    h2R::B
    bcR::A
    ER::C
    BR::C
    lorenzR::D

    bcU::A
    bcD::A

    function IB3F(
        wL::AbstractArray,
        primL::AbstractArray,
        h0L::AbstractArray,
        h1L::AbstractArray,
        h2L::AbstractArray,
        bcL::AbstractArray,
        EL::AbstractArray,
        BL::AbstractArray,
        lorenzL::AbstractArray,
        wR::AbstractArray,
        primR::AbstractArray,
        h0R::AbstractArray,
        h1R::AbstractArray,
        h2R::AbstractArray,
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
            bcL,
            EL,
            BL,
            lorenzL,
            wR,
            primR,
            h0R,
            h1R,
            h2R,
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
Initial & boundary condition with 4 distribution functions

    @consts: wL, primL, h0L, h1L, h2L, h3L, bcL, EL, BL, lorenzL, wR, primR, h0R, h1R, h2R, h3R, bcR, ER, BR, lorenzR, bcU, bcD

"""
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