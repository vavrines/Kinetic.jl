import Pkg
Pkg.add("Documenter")
Pkg.add("KitFort")

using Documenter, Kinetic, KitFort
using Kinetic: KitBase, KitML

tutorial_page = [
    "Examples" => "tutorial.md",
    "Advection-diffusion" => "eg_advection.md",
    "Shock tube" => "eg_shock.md",
    "Lid-driven cavity" => "eg_cavity.md",
]

type_page = [
    "Configuration" => "type_config.md",
    "Setup" => "type_setup.md",
    "Domain" => "type_domain.md",
    "Velocity" => "type_velocity.md",
    "Property" => "type_property.md",
    "Condition" => "type_ib.md",
    "One-dimensional FVM" => "type_1d.md",
    "Two-dimensional FVM" => "type_2d.md",
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

utility_page = [
    "I / O" => "api_io.md",
    "Math" => "api_math.md",
    "Theory" => "api_theory.md",
    "Physical space" => "api_geo.md",
    "Phase space" => "api_phase.md",
    "Configuration" => "api_config.md",
    "Stepper" => "api_step.md",
]

fortran_page = [
    "KitFort.jl" => "fortran1.md",
    "Benchmark" => "fortran2.md",
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
        "Physics" => "physics.md",
        "Type" => type_page,
        "Solver" => solver_page,
        "Tutorial" => tutorial_page,
        "Parallelization" => "parallel.md",
        "Utility" => utility_page,
        "SciML" => "kitml.md",
        "Fortran" => fortran_page,
        "Index" => "function_index.md",
        "Python" => "python.md",
        "Contribution" => "contribution.md",
        "Reference" => "reference.md",
    ],
    format = format,
)

deploydocs(
    repo = "github.com/vavrines/Kinetic.jl.git",
)
