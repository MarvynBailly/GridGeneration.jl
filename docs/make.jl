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
            "ODE Formulation" => "pages/ODE/ODEFormulation.md",
            "Mathematical Work" => "pages/ODE/MathematicalWork.md",
            ],
        
        "Numerical Methods" => Any[
            "First Order System" => "pages/NumericalMethods/FirstOrderSystem.md",
            "Second Order BVP ODE" => "pages/NumericalMethods/SecondOrderBVP.md",
            "Semi-Analytical Method" => "pages/NumericalMethods/SemiAnalyticalMethod.md",
        ],

        "2D to 1D Reformulation " => Any[
            "Mapping 2D to 1D" => "pages/2Dto1D/Mapping2Dto1D.md",
            "Metric Reformulation" => "pages/2Dto1D/MetricReformulation.md",
            "Projecting Points" => "pages/2Dto1D/PointProjection.md",
        ],
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