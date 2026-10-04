#=
Regenerate every figure of the example gallery (docs/src/pages/Examples/gallery.md).

    julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'   # once
    julia --project=examples examples/gallery/run_gallery.jl

Figures are written to docs/src/assets/images/gallery/ and a summary table is printed.
=#
using Logging
include(joinpath(@__DIR__, "gallery_tools.jl"))
for f in ("bump_channel", "annulus", "wavy_channel", "backward_step", "square_metrics", "smoothing")
    include(joinpath(@__DIR__, "cases", f * ".jl"))
end

const OUT = mkpath(joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "images", "gallery"))
global_logger(ConsoleLogger(stderr, Logging.Error))   # keep the output to the summary

rows = []
for (name, case) in (("bump_channel", bump_channel_case), ("annulus", annulus_case),
                     ("wavy_channel", wavy_channel_case), ("backward_step", backward_step_case))
    c = case()
    kw = haskey(c, :splitRequests) ? (; splitRequests = c.splitRequests) : (;)
    result, stats = run_case(c.input, c.bndInfo, c.interInfo, c.M; params = c.params, kw...)
    initial = c.input isa AbstractVector ? c.input : [c.input]
    gallery_figure(joinpath(OUT, name * ".png"), initial, result, c.M; zoom = get(c, :zoom, nothing))
    push!(rows, name => stats)
end
for (n, s) in square_metrics_figure(joinpath(OUT, "square_metrics.png")); push!(rows, "square: " * n => s); end
for (n, s) in smoothing_figure(joinpath(OUT, "smoothing.png")); push!(rows, "smoothing: " * n => s); end

println("| Case | Blocks | Cells | Mean / max angle deviation | Time |")
println("|---|---|---|---|---|")
for (n, s) in rows
    println("| $n | $(s.blocks) | $(s.cells) | $(round(s.meandev, digits=1))° / $(round(s.maxdev, digits=1))° | $(round(s.seconds, digits=2)) s |")
end
