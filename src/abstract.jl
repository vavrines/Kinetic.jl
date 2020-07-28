# ============================================================
# Abstract Types
# ============================================================


export AbstractSolverSet
export AbstractSetup, AbstractProperty, AbstractCondition
export AbstractVelocitySpace, AbstractPhysicalSpace
export AbstractControlVolume, AbstractControlVolume1D, AbstractControlVolume2D
export AbstractInterface, AbstractInterface1D, AbstractInterface2D
export AbstractSolution, AbstractSolution1D, AbstractSolution2D
export AbstractFlux, AbstractFlux1D, AbstractFlux2D


abstract type AbstractPhysicalSpace end
abstract type AbstractVelocitySpace end
abstract type AbstractSetup end
abstract type AbstractProperty end
abstract type AbstractCondition end
abstract type AbstractSolverSet end

abstract type AbstractControlVolume end
abstract type AbstractInterface end
abstract type AbstractControlVolume1D <: AbstractControlVolume end
abstract type AbstractControlVolume2D <: AbstractControlVolume end
abstract type AbstractInterface1D <: AbstractInterface end
abstract type AbstractInterface2D <: AbstractInterface end

abstract type AbstractSolution end
abstract type AbstractFlux end
abstract type AbstractSolution1D <: AbstractSolution end
abstract type AbstractSolution2D <: AbstractSolution end
abstract type AbstractFlux1D <: AbstractFlux end
abstract type AbstractFlux2D <: AbstractFlux end
