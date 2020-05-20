# ============================================================
# Abstract Types
# ============================================================


export AbstractVelocitySpace,
    AbstractPhysicalSpace,
    AbstractSetup,
    AbstractProperty,
    AbstractCondition,
    AbstractSolverSet,
    AbstractControlVolume,
    AbstractControlVolume1D,
    AbstractInterface,
    AbstractInterface1D,
    AbstractSolution


abstract type AbstractPhysicalSpace end
abstract type AbstractVelocitySpace end
abstract type AbstractSetup end
abstract type AbstractProperty end
abstract type AbstractCondition end
abstract type AbstractSolverSet end
abstract type AbstractControlVolume end
abstract type AbstractInterface end

abstract type AbstractControlVolume1D <: AbstractControlVolume end
abstract type AbstractInterface1D <: AbstractInterface end

abstract type AbstractSolution end
