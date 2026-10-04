using Documenter
using GridGeneration


DocMeta.setdocmeta!(GridGeneration, :DocTestSetup, :(using GridGeneration); recursive=true)

# --- Theory pages: generated from math_notes/ (the single source) -----------------------
# Converts $$...$$ to ```math blocks and $...$ to ``...`` (Documenter's math syntax),
# leaving fenced code blocks untouched.
const THEORY_NOTES = [
    "steger_sorenson_theory.md",
    "elliptic_smoothing.md",
    "comparison_steger_sorenson.md",
    "steger_sorenson_implementation_spec.md",
]
const THEORY_PREAMBLE = Dict(
    "steger_sorenson_implementation_spec.md" => """
        !!! note "Planned design"
            This is a design specification for a line Gauss-Seidel Steger & Sorenson smoother.
            It is not implemented yet; the functions and `SSParams` type it describes do not exist
            in the package. The current smoother is described in the elliptic smoothing page.

        """,
)

function convert_math(text)
    out = IOBuffer()
    infence = false
    chunk = IOBuffer()   # prose accumulated between code fences
    flush_prose() = begin
        s = String(take!(chunk))
        s = replace(s, r"\$\$(.+?)\$\$"s => m -> "\n```math\n" * strip(chop(m; head=2, tail=2)) * "\n```\n")
        s = replace(s, r"(?<![\$\\])\$([^\$\n]+?)\$" => m -> "``" * chop(m; head=1, tail=1) * "``")
        print(out, s)
    end
    for line in eachline(IOBuffer(replace(text, "\r\n" => "\n")); keep=true)
        if startswith(lstrip(line), "```")
            infence || flush_prose()
            infence = !infence
            print(out, line)
        elseif infence
            print(out, line)
        else
            print(chunk, line)
        end
    end
    flush_prose()
    return String(take!(out))
end

let src = joinpath(@__DIR__, "..", "math_notes"), dst = joinpath(@__DIR__, "src", "pages", "Theory")
    mkpath(dst)
    for f in THEORY_NOTES
        text = convert_math(read(joinpath(src, f), String))
        lines = split(text, '\n'; limit=2)   # keep the H1 title first
        body = get(THEORY_PREAMBLE, f, "")
        write(joinpath(dst, f), lines[1] * "\n\n" * body * (length(lines) > 1 ? lines[2] : ""))
    end
end

makedocs(
    modules  = [GridGeneration],
    sitename = "GridGeneration.jl",
    authors  = "Marvyn Bailly",
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://MarvynBailly.github.io/GridGeneration.jl/stable/",
        assets=String[],
    ),
    
    
    
    pages = [
        "Home" => "index.md",
        "Getting Started" => "pages/GettingStarted.md",
        "Examples" => Any[
            "Gallery" => "pages/Examples/gallery.md",
            "Airfoil" => "pages/Examples/airfoil.md",
        ],

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

        "Grid Format" => Any["Grid Format" => "pages/GridFormat.md"],

        "Single Block Grid Input" => Any[
            "pages/SingleBlock/nosplitting.md",
            "pages/SingleBlock/splitting.md"
        ],

        "Multi-Block Grid Input" => Any["Multi-Block Input" => "pages/MultiBlock/multiblock.md"],

        "Theory: Elliptic Smoothing" => Any[
            "Steger & Sorenson Theory" => "pages/Theory/steger_sorenson_theory.md",
            "Elliptic Smoothing (current)" => "pages/Theory/elliptic_smoothing.md",
            "Comparison with Steger & Sorenson" => "pages/Theory/comparison_steger_sorenson.md",
            "Line Gauss-Seidel Spec (planned)" => "pages/Theory/steger_sorenson_implementation_spec.md",
        ],

        "API Reference" => "pages/api.md",
    ],

    # Optional quality gates once you’re ready:
    checkdocs = :exports,
)

deploydocs(
    repo      = "github.com/MarvynBailly/GridGeneration.jl",
    devbranch = "main",
    versions = ["stable" => "v^", "v#.#.#"] 
)