# ============================================================
# Abstract Types
# ============================================================


export AbstractVelocityMesh,
       AbstractPhysicalMesh,
       AbstractSetup,
       AbstractProperty,
       AbstractCondition,
       AbstractSolverSet,
       AbstractControlVolume,
       AbstractControlVolume1D,
       AbstractInterface,
       AbstractInterface1D


abstract type AbstractPhysicalMesh end      
abstract type AbstractVelocityMesh end
abstract type AbstractSetup end
abstract type AbstractProperty end
abstract type AbstractCondition end
abstract type AbstractSolverSet end
abstract type AbstractControlVolume end
abstract type AbstractInterface end

abstract type AbstractControlVolume1D <: AbstractControlVolume end
abstract type AbstractInterface1D <: AbstractInterface end