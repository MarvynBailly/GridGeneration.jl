""" 
Function to convert a 2D boundary to a 1D representation
"""
function Boundary2Dto1D(boundary)
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