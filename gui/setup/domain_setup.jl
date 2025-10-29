"""
Domain setup and initialization for the Grid Generation GUI.
"""

"""
    setup_rectangle_domain(width, height, num_points_width, num_points_height)

Create a simple rectangular domain with TFI initial grid.

# Arguments
- `width`: Width of the rectangle
- `height`: Height of the rectangle
- `num_points_width`: Number of grid points in the width direction
- `num_points_height`: Number of grid points in the height direction

# Returns
- `initialGrid`: Vector containing the initial block
- `initialBndInfo`: Boundary information dictionary
- `initialInterfaceInfo`: Interface information dictionary (empty for single block)
"""
function setup_rectangle_domain(width, height, num_points_width, num_points_height)
    # Define boundary curves
    top = hcat(range(0, stop=width, length=num_points_width), fill(height, num_points_width))
    bottom = hcat(range(0, stop=width, length=num_points_width), fill(0.0, num_points_width))
    right = hcat(fill(width, num_points_height), range(0, stop=height, length=num_points_height))
    left = hcat(fill(0.0, num_points_height), range(0, stop=height, length=num_points_height))
    
    # Create initial block using TFI
    initial_block = GridGeneration.TFI([top, right, bottom, left])
    
    # Setup boundary information
    initialBndInfo = []
    push!(initialBndInfo, Dict(
        "name" => "bottom", 
        "faces" => [Dict("block" => 1, "start" => [1,1,1], "end" => [num_points_width,1,1])]
    ))
    push!(initialBndInfo, Dict(
        "name" => "right", 
        "faces" => [Dict("block" => 1, "start" => [num_points_width,1,1], "end" => [num_points_width,num_points_height,1])]
    ))
    push!(initialBndInfo, Dict(
        "name" => "top", 
        "faces" => [Dict("block" => 1, "start" => [1,num_points_height,1], "end" => [num_points_width,num_points_height,1])]
    ))
    push!(initialBndInfo, Dict(
        "name" => "left", 
        "faces" => [Dict("block" => 1, "start" => [1,1,1], "end" => [1,num_points_height,1])]
    ))
    
    # No interfaces for single block
    initialInterfaceInfo = []
    
    return [initial_block], initialBndInfo, initialInterfaceInfo
end

"""
    setup_turtle_grid_domain(metricFieldFile, gridFolder)

Load domain from Turtle grid files.

# Arguments
- `metricFieldFile`: Path to the metric field file
- `gridFolder`: Path to the folder containing grid files

# Returns
- `initialGrid`: Vector of initial blocks
- `initialBndInfo`: Boundary information dictionary
- `initialInterfaceInfo`: Interface information dictionary
- `M`: Metric function for the domain
"""
function setup_turtle_grid_domain(metricFieldFile, gridFolder)
    metricData, datatype, gridfile = GridGeneration.readTurtleFields(metricFieldFile)
    gridfile = joinpath(gridFolder, gridfile)
    blocks, centers, Xfa, interfaceInfo, bndInfo = GridGeneration.ImportTurtleGrid(gridfile)
    
    tree, refs = GridGeneration.setup_metric_tree(centers)
    M = (x, y) -> GridGeneration.find_nearest_kd(metricData, tree, refs, x, y)
    
    return blocks, bndInfo, interfaceInfo, M
end

"""
    create_default_metric(s=1000)

Create a simple default metric function.

# Arguments
- `s`: Scaling factor for the metric

# Returns
- `M`: Metric function M(x,y) that returns a 2-element vector
"""
function create_default_metric(s=1000)
    return (x, y) -> [s * x^2, s]
end
