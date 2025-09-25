""" 
Function to convert a 2D boundary to a 1D representation
"""
function ProjectBoundary2Dto1D(boundary)
    # boundarySection is 2×N (rows: x,y; columns: points in order)
    x = boundary[1, :]
    y = boundary[2, :]

    # segment lengths (N-1)
    Δx = diff(x)
    Δy = diff(y)
    Δs = sqrt.(Δx.^2 .+ Δy.^2)

    xs = [0.0; cumsum(Δs)]   # length N, xs[1]=0, xs[end]=arclength

    return xs
end


""" 
Function to convert a 1D boundary to a 2D representation
- project `points` onto the `boundary` defined by `xs`
"""
function ProjectBoundary1Dto2D(boundary, sol)
    N = size(boundary, 2)
    # step 1: arc length cumulative xs
    seglen = sqrt.((diff(boundary[1, :])).^2 .+ (diff(boundary[2, :])).^2)
    xs = [0.0; cumsum(seglen)]
    totalL = xs[end]

    # sanity check: sol should be in [0, totalL]
    @assert all(sol .>= 0) && all(sol .<= totalL + 1e-12)

    projected = zeros(2, length(sol))

    j = 1
    for (k, s) in enumerate(sol)
        # step 3: find interval
        while j < length(xs) && xs[j+1] < s
            j += 1
        end
        if j == length(xs)
            projected[:, k] .= boundary[:, end]
        else
            # step 4: relative position within interval
            θ = (s - xs[j]) / (xs[j+1] - xs[j])
            # step 5: interpolate boundary point
            projected[:, k] .= (1-θ) * boundary[:, j] + θ * boundary[:, j+1]
        end
    end

    return projected
end
