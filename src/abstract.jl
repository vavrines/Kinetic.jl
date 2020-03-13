# ============================================================
# Abstract Types
# ============================================================


export AbstractVelocityMesh,
       AbstractPhysicalMesh,
       AbstractSetup,
       AbstractProperty,
       AbstractCondition,
       AbstractSolverSet


abstract type AbstractPhysicalMesh end      
abstract type AbstractVelocityMesh end
abstract type AbstractSetup end
abstract type AbstractProperty end
abstract type AbstractCondition end
abstract type AbstractSolverSet end