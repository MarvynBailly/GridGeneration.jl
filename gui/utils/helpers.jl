"""
Helper utility functions for the Grid Generation GUI.
"""

"""
    to_linesegments(X, Y)

Convert 2D grid arrays to line segments for plotting.
Returns a vector of Point2f pairs representing grid lines.
"""
function to_linesegments(X, Y)
    nx, ny = size(X)
    points = Point2f[]
    
    # Horizontal lines
    for j in 1:ny
        for i in 1:nx-1
            push!(points, Point2f(X[i,j], Y[i,j]), Point2f(X[i+1,j], Y[i+1,j]))
        end
    end
    
    # Vertical lines
    for i in 1:nx, j in 1:ny-1
        push!(points, Point2f(X[i,j], Y[i,j]), Point2f(X[i,j+1], Y[i,j+1]))
    end
    
    return points
end

"""
    parse_indices(s::Union{String, Nothing})

Parse a comma-separated string of integers into a vector.
Returns an empty vector if the string is nothing or empty.
Skips any parts that cannot be parsed as integers.

# Examples
```julia
parse_indices("10,25,50")  # Returns [10, 25, 50]
parse_indices("8")          # Returns [8]
parse_indices("")           # Returns []
parse_indices(nothing)      # Returns []
```
"""
function parse_indices(s::Union{String, Nothing})
    indices = Int[]
    isnothing(s) && return indices
    isempty(s) && return indices
    
    parts = split(s, ',', keepempty=false)
    for part in parts
        try
            push!(indices, parse(Int, strip(part)))
        catch e
            println("Warning: Could not parse '$(strip(part))' as an integer. Skipping.")
        end
    end
    
    return indices
end

"""
    log_to_console(console_obs::Observable, message::String)

Append a message to the console observable with timestamp.
"""
function log_to_console(console_obs::Observable, message::String)
    timestamp = Dates.format(now(), "HH:MM:SS")
    current_text = console_obs[]
    new_text = current_text * "[$timestamp] " * message * "\n"
    console_obs[] = new_text
end

"""
    clear_console(console_obs::Observable)

Clear all text from the console.
"""
function clear_console(console_obs::Observable)
    console_obs[] = ""
end
