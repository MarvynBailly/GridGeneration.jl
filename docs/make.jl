using Documenter
using GridGeneration


# Make doctest examples run with `using GridGeneration`
DocMeta.setdocmeta!(GridGeneration, :DocTestSetup, :(using GridGeneration); recursive=true)

makedocs(
    modules  = [GridGeneration],
    sitename = "GridGeneration.jl",
    authors  = "Marvyn Bailly",
    # format = Documenter.HTML(
    #     prettyurls = !isempty(get(ENV, "CI", "")),
    # ),
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://MarvynBailly.github.io/GridGeneration",
        assets=String[],
    ),
    
    
    
    pages = [
        "Home" => "index.md",
        "Ordinary Differential Equations" => Any[
            "ODE Formulation" => "pages/ODE/ODEFormulation.md",
            "Mathematical Work" => "pages/ODE/MathematicalWork.md",
            ],

        "Examples" => Any[
            "Airfoil" => "pages/Examples/airfoil.md",
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

        "Grid Format" => Any["Grid Format" => "pages/GridFormat.md"],

        "Single Block Grid Input" => Any[
            "pages/SingleBlock/nosplitting.md",
            "pages/SingleBlock/splitting.md"
        ],

        "Multi-Block Grid Input" => Any["Multi-Block Input" => "pages/MultiBlock/multiblock.md"],
    ],

    # Optional quality gates once you’re ready:
    checkdocs = :exports,
)

deploydocs(
    repo      = "github.com/MarvynBailly/GridGeneration",
    devbranch = "main",
)