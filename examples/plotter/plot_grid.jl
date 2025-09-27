using Plots

function plot_grid(x, y, title_str; skip = 2, plt = nothing, c = nothing)
    p = plt == nothing ? plot() : plt
    plot!(p, title=title_str, aspect_ratio=:equal, legend=false, framestyle=:box, grid = false)
    # Plot η-lines (lines of constant j)
    for j in 1:size(x, 2)
        plot!(p, x[:, j], y[:, j], lw=0.1, color=c)
    end
    # Plot ξ-lines (lines of constant i)
    for i in 1:skip:size(x, 1)
        plot!(p, x[i, :], y[i, :], lw=0.1, color=c)
    end
    return p
end