push!(LOAD_PATH, "..")

import Pkg
Pkg.add("Documenter")

using Documenter, Kinetic
using Kinetic: KitBase, KitML
using Kinetic.KitBase: ib_rh, ib_sod, ib_briowu, ib_cavity
using Kinetic.KitBase: newton_cotes, triangle_weights
using Kinetic.KitBase: em_coefficients, advection_flux, burgers_flux, euler_flux, euler_jacobi
using Kinetic.KitBase: gauss_moments, mixture_gauss_moments, pdf_slope, mixture_pdf_slope, moments_conserve_slope, mixture_moments_conserve_slope
using Kinetic.KitBase: aap_hs_diffeq!
using Kinetic.KitBase: plot_line, plot_contour, write_jld
using Kinetic.KitML.Solaris
using Kinetic.KitML.Solaris: load_data, save_model
using Kinetic.KitBase.FiniteMesh
using Kinetic.KitBase.FiniteMesh: read_mesh, mesh_connectivity_2D, mesh_cell_center, mesh_cell_area_2D

DocMeta.setdocmeta!(KitBase, :DocTestSetup, :(using KitBase); recursive = true)
DocMeta.setdocmeta!(KitML, :DocTestSetup, :(using KitML); recursive = true)

tutorial_page = [
    "Examples" => "tutorial.md",
    "Advection diffusion" => "eg_advection.md",
    "Burgers" => "eg_burgers.md",
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
    "FVM" => "type_fvm.md",
]

solver_page = [
    "Framework" => "solver.md",
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

parallel_page = [
    "General" => "parallel.md",
    "Threading" => "para_thread.md",
    "Distributed" => "para_dist.md",
    "CUDA" => "para_cuda.md",
]

ml_page = ["KitML" => "kitml1.md", "UBE" => "kitml2.md"]

fortran_page = ["KitFort" => "fortran1.md", "Benchmark" => "fortran2.md"]

format = Documenter.HTML(assets = ["assets/favicon.ico"], collapselevel = 1)

makedocs(
    sitename = "Kinetic.jl",
    modules = [Kinetic, KitBase, KitML, FiniteMesh, Solaris],
    pages = Any[
        "Home"=>"index.md",
        "Installation"=>"install.md",
        "Physics"=>"physics.md",
        "Type"=>type_page,
        "Solver"=>solver_page,
        "Tutorial"=>tutorial_page,
        "Parallelization"=>parallel_page,
        "Utility"=>utility_page,
        "SciML"=>ml_page,
        "Fortran"=>fortran_page,
        "Index"=>"function_index.md",
        "Python"=>"python.md",
        "Contribution"=>"contribution.md",
        "Reference"=>"reference.md",
    ],
    format = format,
    checkdocs = :none,
)

deploydocs(repo = "github.com/vavrines/Kinetic.jl.git")
