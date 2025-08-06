# Build M(x) as linear interpolation, Mx(x) as piecewise-constant slope.
function build_interps_linear(xv, Mv)
    @assert length(xv) == length(Mv) ≥ 2
    @assert issorted(xv)

    # precompute slopes on each interval
    Δx = diff(xv)

    # linear interpolation with clamping to [xv[1], xv[end]]
    M_of = function (x::Real)
        if x ≤ xv[1]; return Mv[1] end
        if x ≥ xv[end]; return Mv[end] end
        i = searchsortedlast(xv, x)       # i such that xv[i] ≤ x < xv[i+1]
        t = (x - xv[i]) / Δx[i]
        return (1 - t)*Mv[i] + t*Mv[i+1]
    end

    return M_of
end

function InterpolateBoundaryManually(x_spacing, boundary)
    # boundary: 2×N (x row, y row)
    x_bnd = collect(boundary[1, :])
    y_bnd = collect(boundary[2, :])
    n = length(x_bnd)
    @assert size(boundary, 1) == 2 && n ≥ 2 "boundary must be 2×N with N≥2"

    # Step 1: cumulative arc length along boundary order
    arc_lengths = zeros(eltype(x_bnd), n)
    @inbounds for i in 2:n
        dx = x_bnd[i] - x_bnd[i-1]
        dy = y_bnd[i] - y_bnd[i-1]
        arc_lengths[i] = arc_lengths[i-1] + hypot(dx, dy)
    end
    L = arc_lengths[end]
    @assert L > 0 "boundary has zero total length"

    # Step 2: map x_spacing to [0, L]
    xs = (x_spacing .- first(x_spacing)) ./ (last(x_spacing) - first(x_spacing)) .* L

    # Clamp to [0,L] to avoid tiny numerical overshoots
    @inbounds for j in eachindex(xs)
        xs[j] = clamp(xs[j], 0, L)
    end

    # Step 3: piecewise linear interpolation along arc length
    projected_x = similar(xs, eltype(x_bnd))
    projected_y = similar(xs, eltype(y_bnd))

    @inbounds for j in eachindex(xs)
        s = xs[j]
        i = searchsortedlast(arc_lengths, s)
        if i ≥ n
            projected_x[j] = x_bnd[end]
            projected_y[j] = y_bnd[end]
            continue
        elseif i == 0
            projected_x[j] = x_bnd[1]
            projected_y[j] = y_bnd[1]
            continue
        end

        s0 = arc_lengths[i]
        s1 = arc_lengths[i+1]
        denom = s1 - s0
        t = denom == 0 ? 0.0 : (s - s0) / denom 

        projected_x[j] = (1 - t)*x_bnd[i] + t*x_bnd[i+1]
        projected_y[j] = (1 - t)*y_bnd[i] + t*y_bnd[i+1]
    end

    return projected_x, projected_y
end
