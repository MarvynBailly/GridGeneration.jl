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
Project 1D solution points onto 2D boundary curve.
Converts arc-length parametrization back to physical coordinates.
"""
function ProjectBoundary1Dto2D(boundary, sol)
    N = size(boundary, 2)
    
    # Compute cumulative arc length
    seglen = sqrt.((diff(boundary[1, :])).^2 .+ (diff(boundary[2, :])).^2)
    xs = [0.0; cumsum(seglen)]
    totalL = xs[end]

    @assert all(sol .>= 0) && all(sol .<= totalL + 1e-12)

    projected = zeros(2, length(sol))

    j = 1
    for (k, s) in enumerate(sol)
        # Find interval containing s
        while j < length(xs) && xs[j+1] < s
            j += 1
        end
        if j == length(xs)
            projected[:, k] .= boundary[:, end]
        else
            # Linear interpolation within interval
            θ = (s - xs[j]) / (xs[j+1] - xs[j])
            projected[:, k] .= (1-θ) * boundary[:, j] + θ * boundary[:, j+1]
        end
    end

    return projected
end
