# ============================================================
# Kinetic.jl : Kinetic Theory, Modeling and Simulation Toolbox
# ============================================================

module Kinetic

using Dates
using OffsetArrays
using SpecialFunctions
using FFTW
using OrdinaryDiffEq
using FileIO
using JLD2
using Plots
using PyCall

include("Data/data.jl")
include("IO/io.jl")
include("Math/math.jl")
include("Geometry/geometry.jl")
include("Theory/theory.jl")
include("Phase/phase.jl")
include("Reconstruction/reconstruction.jl")
include("Flux/flux.jl")
include("Config/config.jl")
include("Solver/solver.jl")

end