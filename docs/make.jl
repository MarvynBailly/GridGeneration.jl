# docs/make.jl
using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

# Optional hot-reload during local dev
try
    using Revise
catch
end

using Documenter
using GridGeneration


# Make doctest examples run with `using GridGeneration`
DocMeta.setdocmeta!(GridGeneration, :DocTestSetup, :(using GridGeneration); recursive=true)

makedocs(
    modules  = [GridGeneration],
    sitename = "GridGeneration.jl",
    authors  = "Marvyn Bailly",
    format   = Documenter.HTML(
        prettyurls = !isempty(get(ENV, "CI", "")),  # false locally (Windows), true on CI
    ),
    
    
    
    pages = [
        "Home" => "index.md",
        "Ordinary Differential Equations" => Any[
            "ODE Formulation" => "pages/ODEFormulation.md",
            "ODE Numerical Methods"  => "pages/ODENumericalMethods.md",
            ],
        
        "1D Projection onto 2D" => Any[
            "Distribution of Points" => "pages/DistributionOfPoints.md",
        ],

        "Appendix" => Any[
            "Mathematical Work" => "pages/MathematicalWork.md",
        ]
    ],



    # Optional quality gates once you’re ready:
    # doctest  = true,
    checkdocs = :exports,
)

deploydocs(
    repo      = "github.com/MarvynBailly/GridGeneration.jl",
    devbranch = "main",
    # versions = ["stable" => "v^", "dev" => "main"],  # enable when you start tagging
)