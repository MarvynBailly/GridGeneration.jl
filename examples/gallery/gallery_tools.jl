#=
Shared helpers for the example gallery: building blocks from boundary curves, a few
reusable metric "ingredients", and the figure layout used for every case.
=#
using GridGeneration
using Plots

# ----------------------------------------------------------------------------------------
# Domains
# ----------------------------------------------------------------------------------------

"""
    block_from_curves(bottom, top, left, right; Ni, Nj)

Build a `[2, Ni, Nj]` block by transfinite interpolation from four parametric curves
`c(t) -> (x, y)`, `t ∈ [0, 1]`. `bottom`/`top` run left→right and `left`/`right` run
bottom→top, with matching corners.
"""
function block_from_curves(bottom, top, left, right; Ni, Nj)
    edge(c, n) = reduce(vcat, [collect(c(t))' for t in range(0, 1, length=n)])
    return TFI([edge(top, Ni), edge(right, Nj), edge(bottom, Ni), edge(left, Nj)])
end

line(p, q) = t -> (p[1] + t * (q[1] - p[1]), p[2] + t * (q[2] - p[2]))

# ----------------------------------------------------------------------------------------
# Metric ingredients: each returns M(x, y) -> (M11, M22); combine them with `metric_sum`
# ----------------------------------------------------------------------------------------

metric_sum(fs...) = (x, y) -> reduce((a, b) -> a .+ b, (f(x, y) for f in fs))

"Constant (isotropic) background metric: target spacing 1/sqrt(m)."
background(m) = (x, y) -> (m, m)

"Gaussian hotspot of amplitude `A` and radius `r` centred at `c`."
hotspot(c, A, r) = (x, y) -> (w = A * exp(-((x - c[1])^2 + (y - c[2])^2) / r^2); (w, w))

"Refinement along the line through `p` with direction `d` (e.g. an oblique shock), width `w`."
function feature_line(p, d, A, w)
    n = (-d[2], d[1]) ./ hypot(d...)
    return (x, y) -> (v = A * exp(-(((x - p[1]) * n[1] + (y - p[2]) * n[2]) / w)^2); (v, v))
end

"Refinement along the circle of radius `R` about `c` (a front), width `w`."
front(c, R, A, w) = (x, y) -> (v = A * exp(-((hypot(x - c[1], y - c[2]) - R) / w)^2); (v, v))

"""
    wall_layer(dist, A, ℓ; normal=:y)

Boundary-layer clustering `A exp(-d/ℓ)` with `d = dist(x, y)` the distance to a wall.
`normal = :both` refines isotropically; `normal = :x` or `:y` refines only that metric
component, so spacing shrinks across the wall but not along it (an anisotropic metric).
"""
function wall_layer(dist, A, ℓ; normal=:both)
    return function (x, y)
        v = A * exp(-dist(x, y) / ℓ)
        normal === :x ? (v, 0.0) : normal === :y ? (0.0, v) : (v, v)
    end
end

# ----------------------------------------------------------------------------------------
# Running a case and drawing it
# ----------------------------------------------------------------------------------------

"""
    run_case(input, bndInfo, interInfo, M; params, splitRequests=nothing)

Run `GenerateGrid` (single block or multi-block) and collect summary statistics.
"""
function run_case(input, bndInfo, interInfo, M; params, splitRequests=nothing)
    t = @elapsed result = splitRequests === nothing ?
        GenerateGrid(input, bndInfo, interInfo, M; params) :
        GenerateGrid(input, bndInfo, interInfo, M; params, splitRequests)
    devs, _, maxdev, _, _ = ComputeAngleDeviation(result.smoothBlocks)
    meandev = sum(sum(filter(isfinite, d[1:end-1, 1:end-1])) for d in devs) /
              sum(length(d[1:end-1, 1:end-1]) for d in devs)
    stats = (blocks = length(result.smoothBlocks),
             cells = sum((size(b, 2) - 1) * (size(b, 3) - 1) for b in result.smoothBlocks),
             maxdev = maxdev, meandev = meandev, seconds = t)
    return result, stats
end

# grid lines of a set of blocks, block outlines drawn thicker
function draw_blocks!(p, blocks; lw = 0.35, color = :black, outline = :crimson, outline_lw = 1.4, skip = 1)
    for b in blocks
        x, y = b[1, :, :], b[2, :, :]
        for j in 1:skip:size(x, 2); plot!(p, x[:, j], y[:, j]; lw, color, label = ""); end
        for i in 1:skip:size(x, 1); plot!(p, x[i, :], y[i, :]; lw, color, label = ""); end
    end
    if outline !== nothing
        for b in blocks
            x, y = b[1, :, :], b[2, :, :]
            ox = [x[:, 1]; x[end, :]; reverse(x[:, end]); reverse(x[1, :])]
            oy = [y[:, 1]; y[end, :]; reverse(y[:, end]); reverse(y[1, :])]
            plot!(p, ox, oy; lw = outline_lw, color = outline, label = "")
        end
    end
    return p
end

block_outline(b) = (x = b[1, :, :]; y = b[2, :, :];
    ([x[:, 1]; x[end, :]; reverse(x[:, end]); reverse(x[1, :])],
     [y[:, 1]; y[end, :]; reverse(y[:, end]); reverse(y[1, :])]))

function inside_polygon(px, py, xs, ys)
    inside = false
    j = length(xs)
    @inbounds for i in eachindex(xs)
        if ((ys[i] > py) != (ys[j] > py)) && (px < (xs[j] - xs[i]) * (py - ys[i]) / (ys[j] - ys[i]) + xs[i])
            inside = !inside
        end
        j = i
    end
    return inside
end

# metric magnitude sqrt(M11 + M22) as a heatmap, masked to the union of the blocks
function draw_metric!(p, blocks, M; colormap = :inferno, npix = 700)
    outlines = map(block_outline, blocks)
    xmin = minimum(minimum(o[1]) for o in outlines); xmax = maximum(maximum(o[1]) for o in outlines)
    ymin = minimum(minimum(o[2]) for o in outlines); ymax = maximum(maximum(o[2]) for o in outlines)
    nx = npix; ny = max(2, round(Int, npix * (ymax - ymin) / (xmax - xmin)))
    xs = range(xmin, xmax, length = nx); ys = range(ymin, ymax, length = ny)
    z = fill(NaN, ny, nx)
    for (jx, x) in enumerate(xs), (jy, y) in enumerate(ys)
        if any(o -> inside_polygon(x, y, o[1], o[2]), outlines)
            z[jy, jx] = log10(sqrt(sum(M(x, y))))
        end
    end
    heatmap!(p, xs, ys, z; c = colormap, colorbar = true, colorbar_title = "log₁₀ √(M₁₁+M₂₂)")
    return p
end

"""
    gallery_figure(path, title, initialBlocks, result, M; zoom=nothing, kwargs...)

Save a figure with three panels: the metric field over the domain with the initial block
layout, the final grid, and (optionally) a zoomed view of the final grid
(`zoom = (xlims, ylims)`).
"""
function gallery_figure(path, initialBlocks, result, M; zoom = nothing, skip = 1, width = 1500)
    outlines = map(block_outline, result.smoothBlocks)
    xmin = minimum(minimum(o[1]) for o in outlines); xmax = maximum(maximum(o[1]) for o in outlines)
    ymin = minimum(minimum(o[2]) for o in outlines); ymax = maximum(maximum(o[2]) for o in outlines)
    aspect = (xmax - xmin) / (ymax - ymin)
    pad = 0.02 * max(xmax - xmin, ymax - ymin)
    lims = (xlims = (xmin - pad, xmax + pad), ylims = (ymin - pad, ymax + pad))
    kw = (aspect_ratio = :equal, framestyle = :box, grid = false, legend = false,
          titlefontsize = 12, tickfontsize = 8, guidefontsize = 9)

    p1 = plot(; title = "Metric (log scale) and block layout", lims..., kw...)
    draw_metric!(p1, result.smoothBlocks, M)
    for b in result.smoothBlocks
        plot!(p1, block_outline(b)...; color = :white, lw = 1.0, label = "")
    end

    p2 = plot(; title = "Generated grid ($(length(result.smoothBlocks)) blocks)", lims..., kw...)
    draw_blocks!(p2, result.smoothBlocks; skip)

    panels = Any[p1, p2]
    if zoom !== nothing
        p3 = plot(; title = "Detail", xlims = zoom[1], ylims = zoom[2], kw...)
        draw_blocks!(p3, result.smoothBlocks; lw = 0.6, outline_lw = 1.6)
        push!(panels, p3)
    end

    n = length(panels)
    if aspect > 1.8      # wide domain: stack the panels
        h = round(Int, width / aspect * 1.05) + 70
        fig = plot(panels...; layout = (n, 1), size = (width, n * h), left_margin = 4Plots.mm)
    else
        w = width ÷ n
        fig = plot(panels...; layout = (1, n), size = (width, round(Int, w / aspect * 0.95) + 90),
                   left_margin = 4Plots.mm, bottom_margin = 4Plots.mm)
    end
    savefig(fig, path)
    return fig
end

print_stats(name, s) = println(rpad(name, 34), " blocks=$(s.blocks)  cells=$(s.cells)  ",
    "max angle dev=$(round(s.maxdev, digits=1))°  mean=$(round(s.meandev, digits=2))°  time=$(round(s.seconds, digits=1)) s")
