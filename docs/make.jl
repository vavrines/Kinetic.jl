#push!(LOAD_PATH,"../src/")
import Pkg
Pkg.add("Documenter")

using Documenter, Kinetic
using Kinetic: KitBase, KitML

makedocs(
    sitename= "Kinetic.jl",
    modules = [Kinetic, KitBase, KitML],
    pages = Any[
        "Home" => "index.md",
        "Installation" => "install.md",
        "Tutorial" => "tutorial.md",
        "Physics" => "physics.md",
        "Type" => "type.md",
        "Solver" => [
            "solver.md",
            "solver_pre.md",
            "solver_timestep.md",
            "solver_reconstruction.md",
            "solver_flux.md",
            "solver_update.md",
            "solver_post.md",
            ],
        "SciML" => "sciml.md",
        "Auxiliary" => [
            "api_io.md",
            "api_math.md",
            "api_theory.md",
            "api_geo.md",
            "api_phase.md",
            "api_config.md",
            "api_step.md",
            ],
    ]
)

deploydocs(
    repo = "github.com/vavrines/Kinetic.jl.git",
)
