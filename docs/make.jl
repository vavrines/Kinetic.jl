import Pkg
Pkg.add("Documenter")
Pkg.add("KitFort")

using Documenter, Kinetic, KitFort
using Kinetic: KitBase, KitML

tutorial_page = [
    "Shock tube problem" => "eg_shock.md",
    "Lid-driven cavity" => "eg_cavity.md",
]

type_page = [
    "General" => "type.md",
    "One-dimensional data" => "type_1d.md",
    "Two-dimensional data" => "type_2d.md",
]

solver_page = [
    "General" => "solver.md",
    "Preprocess" => "solver_pre.md",
    "Timestep" => "solver_timestep.md",
    "Reconstruction" => "solver_reconstruction.md",
    "Flux" => "solver_flux.md",
    "Update" => "solver_update.md",
    "Postprocess" => "solver_post.md",
]

format = Documenter.HTML(
    collapselevel = 1,
)

makedocs(
    sitename= "Kinetic.jl",
    modules = [Kinetic, KitBase, KitML, KitFort],
    pages = Any[
        "Home" => "index.md",
        "Installation" => "install.md",
        "Tutorial" => tutorial_page,
        "Physics" => "physics.md",
        "Type" => type_page,
        "Solver" => solver_page,
        "Utility" => [
            "api_io.md",
            "api_math.md",
            "api_theory.md",
            "api_geo.md",
            "api_phase.md",
            "api_config.md",
            "api_step.md",
            ],
        "KitML" => "kitml.md",
        "KitFort" => "kitfort.md",
        "Index" => "function_index.md",
        "Contribution" => "contribution.md",
        "Reference" => "reference.md",
    ],
    format = format,
)

deploydocs(
    repo = "github.com/vavrines/Kinetic.jl.git",
)
