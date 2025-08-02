using Documenter, GridGeneration

makedocs(
    sitename = "GridGeneration.jl",
    format = Documenter.HTML(),
    modules  = [GridGeneration],
)

deploydocs(
    repo      = "github.com/MarvynBailly/GridGeneration.git",
)