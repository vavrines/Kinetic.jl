# ============================================================
# Kinetic.jl : Theory, Experiments and Simulations
# ============================================================


module Kinetic


using Dates
using OffsetArrays
using SpecialFunctions
using FFTW
using FileIO
using JLD2
using Plots


include("abstract.jl")
include("io.jl")
include("math.jl")
include("geo.jl")
include("theory.jl")
include("velocity.jl")
include("slope.jl")
include("flux.jl")
include("data.jl")
include("problem.jl")
include("solver.jl")
include("initialize.jl")


end
