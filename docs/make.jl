using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

try
    using Revise
catch
end

using Documenter
using GridGeneration

DocMeta.setdocmeta!(GridGeneration, :DocTestSetup, :(using GridGeneration); recursive=true)

makedocs(
    modules = [GridGeneration],
    sitename = "GridGeneration.jl",
    format = Documenter.HTML(prettyurls = !isempty(get(ENV, "CI", ""))),
    authors = "Marvyn Bailly",
    pages = ["Home" => "index.md", 
             "ODE Formulation" => "pages/ODEFormulation.md",
             "ODE Numerical Methods" => "pages/ODENumericalMethods.md",
             "Distribution of Points" => "pages/DistributionOfPoints.md",],
)

deploydocs(
    repo      = "github.com/MarvynBailly/GridGeneration.jl",
    devbranch = "main",
    # Optional: versioned docs once you start tagging releases
    # versions = ["stable" => "v^", "dev" => "main"],
)
