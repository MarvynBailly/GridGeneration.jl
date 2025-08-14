using Plots
using Plots: Shape

"""
    plot_scalar_field(f, xs, ys; boundary=nothing,
                      boundary_style=:outline,      # :outline | :fill | :both
                      boundary_linecolor=:black,
                      boundary_linewidth=1.5,
                      boundary_linealpha=1.0,
                      boundary_fillcolor=:white,
                      boundary_fillalpha=1.0,
                      mask_inside=false,            # set W=NaN inside polygon
                      title="", xlabel="x", ylabel="y",
                      cb_label="w(x,y)", clim=nothing,
                      colormap=:viridis, equal_aspect=true)

Evaluate and plot a scalar field f(x,y) on grid xs×ys.
If `boundary` (2×N) is provided, it can be outlined and/or filled.
Returns (plot_handle, W::Matrix).
"""
function plot_scalar_field(f, xs::AbstractVector, ys::AbstractVector;
    boundary=nothing,
    boundary_style::Symbol=:outline,
    boundary_linecolor=:black,
    boundary_linewidth::Real=1.5,
    boundary_linealpha::Real=1.0,
    boundary_fillcolor=:white,
    boundary_fillalpha::Real=1.0,
    mask_inside::Bool=false,
    title="", xlabel="x", ylabel="y",
    cb_label="w(x,y)", clim=nothing,
    colormap=:viridis, equal_aspect=true)

    # evaluate f on (ys × xs) grid
    W = [let v=f(x,y); v isa Tuple ? first(v) : v end for y in ys, x in xs]

    # optional masking (NaN-out values inside polygon)
    if mask_inside && boundary !== nothing
        xb = boundary[1, :]; yb = boundary[2, :]
        function point_in_polygon(x, y)
            # ray-casting
            inside = false
            j = length(xb)
            for i in 1:length(xb)
                yi, yj = yb[i], yb[j]
                xi, xj = xb[i], xb[j]
                if ((yi > y) != (yj > y)) &&
                   (x < (xj - xi) * (y - yi) / (yj - yi + eps()) + xi)
                    inside = !inside
                end
                j = i
            end
            return inside
        end
        for (j,y) in enumerate(ys), (i,x) in enumerate(xs)
            if point_in_polygon(x,y)
                W[j,i] = NaN
            end
        end
    end

    # heatmap
    p = heatmap(xs, ys, W; xlabel, ylabel, title,
                c=colormap, colorbar_title=cb_label,
                clim, aspect_ratio=(equal_aspect ? 1 : :auto))

    # boundary overlay
    if boundary !== nothing
        xb = boundary[1, :]; yb = boundary[2, :]
        if boundary_style == :fill || boundary_style == :both
            plot!(p, Shape(xb, yb);
                  color=boundary_fillcolor, alpha=boundary_fillalpha,
                  linecolor=boundary_linecolor, linewidth=boundary_linewidth,
                  label="")
        end
        if boundary_style == :outline || boundary_style == :both
            plot!(p, xb, yb; lw=boundary_linewidth, c=boundary_linecolor,
                  alpha=boundary_linealpha, label="boundary")
        end
    end
    return p, W
end

# Meshgrid-form overload
function plot_scalar_field(f, XX::AbstractMatrix, YY::AbstractMatrix; kwargs...)
    @assert size(XX) == size(YY) "XX and YY must have the same size"
    W = similar(XX, Float64)
    @inbounds for j in axes(XX,1), i in axes(XX,2)
        v = f(XX[j,i], YY[j,i]); W[j,i] = v isa Tuple ? first(v) : v
    end
    xs = range(extrema(XX)[1], extrema(XX)[2], length=size(XX, 2))
    ys = range(extrema(YY)[1], extrema(YY)[2], length=size(XX, 1))
    plot_scalar_field((x,y)->W[searchsortedfirst(ys,y), searchsortedfirst(xs,x)],
                      xs, ys; kwargs...)
end



"""
    plot_scalar_field(f, xs, ys; boundary=nothing, title="",
                      xlabel="x", ylabel="y", cb_label="w(x,y)",
                      clim=nothing, colormap=:viridis, equal_aspect=true)

Evaluate and plot a scalar field f(x,y) on the grid defined by vectors xs, ys.
Returns (plot_handle, W::Matrix).
"""
function plot_scalar_field(f, xs::AbstractVector, ys::AbstractVector;
                           boundary=nothing, title="", xlabel="x", ylabel="y",
                           cb_label="w(x,y)", clim=nothing,
                           colormap=:viridis, equal_aspect=true)

    # Evaluate f on the (ys × xs) grid (rows = y, cols = x)
    W = [let val = f(x, y); val isa Tuple ? first(val) : val end
         for y in ys, x in xs]

    p = heatmap(xs, ys, W; xlabel, ylabel, title,
                c=colormap, colorbar_title=cb_label,
                clim, aspect_ratio=(equal_aspect ? 1 : :auto))

    if boundary !== nothing
        plot!(p, boundary[1, :], boundary[2, :], lw=1.5, c=:black, label="boundary")
    end
    return p, W
end
