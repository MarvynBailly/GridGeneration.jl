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

    # normalize 
    xs = xs ./ xs[end]  # now xs is in [0, 1]
    return xs
end

""" 
Function to convert a 1D boundary to a 2D representation
- project `points` onto the `boundary` defined by `xs`
"""
function ProjectBoundary1Dto2D(boundary, points, xs)
    projectedPoints = zeros(2, length(points))

    for (i, pnt) in enumerate(points)
        # clamp
        if i == 1
            intervalIndex = 1
            normalDist = 0
        elseif i == length(points)
            intervalIndex = length(xs) - 1
            normalDist = 1
        else
            intervalIndex = FindContainingIntervalIndex(pnt, xs)
            normalDist = ComputeNormalDistance(pnt, xs, intervalIndex)
        end


        projectPoint = ProjectPointOntoBoundary(normalDist, intervalIndex, boundary)

        projectedPoints[:, i] = projectPoint

        # println("Projecting point $i: $pnt onto boundary at index $intervalIndex with distance $dist, normalDist = $normalDist")

    end
    return projectedPoints
end

function ComputeNormalDistance(pnt, dist, int)
    # assert that dist is increasing
    @assert int >= 1 && int < length(dist) "int must be in range [1, $(length(dist) - 1)]"
    return (pnt - dist[int]) / (dist[int + 1] - dist[int])
end


function FindContainingIntervalIndex(pnt, dist)
    # assert that dist is increasing
    
    # clamp
    if pnt < dist[1]
        return  1
    end

    if pnt > dist[end]
        return length(dist) - 1
    end

    for (i,d) in enumerate(dist)
        if d > pnt
            return i - 1
        end
    end

    return -1
end

function ProjectPointOntoBoundary(s, ind, boundary)
    # interpolate boundary[ind] and boundary[ind + 1] with s
    projectedPoint = s * boundary[:,ind + 1] + (1 - s) * boundary[:,ind]
    return projectedPoint
end
