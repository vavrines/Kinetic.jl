# ============================================================
# Data Structure
# ============================================================

export AbstractSolverSet
export AbstractSetup, AbstractProperty, AbstractCondition
export AbstractVelocitySpace, AbstractPhysicalSpace
export AbstractControlVolume, AbstractControlVolume1D, AbstractControlVolume2D
export AbstractInterface, AbstractInterface1D, AbstractInterface2D
export AbstractSolution, AbstractSolution1D, AbstractSolution2D
export AbstractFlux, AbstractFlux1D, AbstractFlux2D

export Setup
export GasProperty, PlasmaProperty
export IB, IB1F, IB2F, IB4F
export ControlVolume1D, ControlVolume1D1F, ControlVolume1D2F, ControlVolume1D4F
export ControlVolume2D, ControlVolume2D1F, ControlVolume2D2F
export Interface1D, Interface1D1F, Interface1D2F, Interface1D4F
export Interface2D, Interface2D1F, Interface2D2F
export Solution1D, Solution1D1F, Solution1D2F
export Solution2D, Solution2D1F, Solution2D2F
export Flux1D, Flux1D1F, Flux1D2F
export Flux2D, Flux2D1F, Flux2D2F

include("abstract.jl")
include("struct.jl")
