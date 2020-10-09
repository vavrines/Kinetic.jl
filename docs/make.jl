#push!(LOAD_PATH,"../src/")
import Pkg
Pkg.add("Documenter")

using Documenter, Kinetic

makedocs(sitename="Kinetic.jl")