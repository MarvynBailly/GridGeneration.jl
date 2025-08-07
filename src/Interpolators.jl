# Build M(x) as linear interpolation, Mx(x) as piecewise-constant slope.
function build_interps_linear(xv, Mv)
    @assert length(xv) == length(Mv) ≥ 2
    @assert issorted(xv)

    # precompute slopes on each interval
    Δx = diff(xv)

    # linear interpolation with clamping to [xv[1], xv[end]]
    M_of = function (x)
        if x ≤ xv[1]; return Mv[1] end
        if x ≥ xv[end]; return Mv[end] end
        i = searchsortedlast(xv, x)       # i such that xv[i] ≤ x < xv[i+1]
        t = (x - xv[i]) / Δx[i]
        return (1 - t)*Mv[i] + t*Mv[i+1]
    end

    return M_of
end

