using Pkg; Pkg.add("Documenter")
using Documenter, Kinetic

makedocs(
    sitename = "Kinetic.jl",
    format = Documenter.HTML(assets = ["assets/myfont.css"]),
    modules = [Kinetic],
    authors = "tianbaixiao@gmail.com",
    doctest = true,
    pages = Any[
        "index.md",
        "type_hierarchy.md",
        "Usage" => [
            "numerical_integration.md",
            "quadrature_rules.md",
            "orthogonal_polynomials_canonical.md",
            "gaussian_mixture_model.md",
            "multiple_discretization.md",
            "scalar_products.md",
            "Polynomial Chaos" => [ "Basic Usage" => "pce_tutorial.md",
                                    "Chi Squared, One DOF" => "chi_squared_k1.md",
                                    "Chi Squared, Several DOFs" => "chi_squared_k_greater1.md",
                                    "Random ODE" => "random_ode.md"
                                  ],
            "Optimal Power Flow" => "DCsOPF.md"
            ],
    "math.md",
    "functions.md"
    ]
)


deploydocs(
    repo = "github.com/vavrines/Kinetic.jl.git",
    target = "build",
)
