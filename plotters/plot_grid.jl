using Plots

function plot_grid(x, y, title_str; skip = 2, plt = nothing, c = :black)
    p = plt == nothing ? plot() : plt
    plot!(p, title=title_str, aspect_ratio=:equal, legend=false, framestyle=:box, grid = false)
    # Plot η-lines (lines of constant j)
    for j in 1:skip:size(x, 2)
        plot!(p, x[:, j], y[:, j], lw=0.1, color=c)
    end
    # Plot ξ-lines (lines of constant i)
    for i in 1:skip:size(x, 1)
        plot!(p, x[i, :], y[i, :], lw=0.1, color=c)
    end
    return p
end

function plot_blocks(blocks, title; plot_boundary=false, boundary_width = 2, boundary_color = :black, skip = 1)
    p = plot(title=title, aspect_ratio=:equal, legend=false, framestyle=:box, grid = false)
    for block in blocks
        plot_grid(block[1, :, :], block[2, :, :], ""; skip=skip, plt=p, c = RGB(0.0, 0.0, 0.0))
        # plot the boundary of the block
        if plot_boundary
            plot!(p, [block[1,1,1], block[1,end,1], block[1,end,end], block[1,1,end], block[1,1,1]],
                  [block[2,1,1], block[2,end,1], block[2,end,end], block[2,1,end], block[2,1,1]], lw=boundary_width, color=boundary_color)
        end
    end
    return p
end