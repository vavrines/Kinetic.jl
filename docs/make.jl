#push!(LOAD_PATH,"../src/")
import Pkg
Pkg.add("Documenter")

using Documenter, Kinetic

makedocs(
    sitename="Kinetic.jl",
    modules = [Kinetic],
    pages = Any[
        "index.md",
        "type.md",
        "algorithm.md",
        "API" => [
            "api_io.md",
            "api_math.md",
            "api_geo.md",
            "api_theory.md",
            "api_phase.md",
            "api_reconstruct.md",
            "api_flux.md",
            "api_config.md",
            "api_solver.md",
            ],
    ]
)

deploydocs(
    repo = "github.com/vavrines/Kinetic.jl.git",
)
