include("C:\\Users\\admin\\Documents\\GitHub\\GridGeneration\\src\\GridGeneration.jl")

using .GridGeneration


using GLMakie
using Makie: Polygon
using DelimitedFiles



function populate_control_panel!(parent_layout::GridLayout)
    widgets = Dict{Symbol, Any}()

    Label(parent_layout[1, 1:2], "Parameters", font = :bold, halign = :center, padding=(10,0,0,10), valign=:top)
    Label(parent_layout[2, 1:2], "Splits Locations", font = :bold, halign = :left, padding=(10,0,10,20), valign=:top)
    Label(parent_layout[3, 1], "I-splits (e.g., 10,25)", halign=:left, valign=:top, padding=(10,0,0,0))
    Label(parent_layout[3, 2], "J-splits (e.g., 8)")
    widgets[:i_splits] = Textbox(parent_layout[4, 1], placeholder="none")
    widgets[:j_splits] = Textbox(parent_layout[4, 2], placeholder="none")

    Label(parent_layout[5, 1:2], "Boundary Solver Parameters", font = :bold, halign = :left, padding=(10,0,10,20), valign=:top)
    Label(parent_layout[6, 1], "Edge Solver Type:", halign=:left, padding=(10,0,0,0))
    widgets[:edge_solver] = Menu(parent_layout[6, 2], options = ["analytic", "numerical"], default = "analytic")

    Label(parent_layout[7, 1:2], "Smoothing Parameters", font = :bold, halign = :left, padding=(10,0,10,20), valign=:top)
    Label(parent_layout[8, 1], "Smoothing Type:", halign=:left, padding =(10,0,0,0))
    widgets[:smoothing_type] = Menu(parent_layout[8, 2], options = ["Elliptic-SS"], default = "Elliptic-SS")

    Label(parent_layout[9, 1], "Max Iterations:", padding=(10,0,0,0))
    widgets[:max_iter] = Textbox(parent_layout[9, 2], placeholder="5000", validator=Int)
    Label(parent_layout[10, 1], "Tolerance:", padding=(10,0,0,0))
    widgets[:tolerance] = Textbox(parent_layout[10, 2], placeholder="1e-5", validator=Float64)

    Label(parent_layout[11, 1], ("Omega: "), halign=:left, padding=(10,0,0,0))
    widgets[:omega] = Textbox(parent_layout[11, 2], placeholder="0.2", validator=Float64)

    Label(parent_layout[12, 1:2], "Wall Forcing Parameters", font = :bold, halign = :left, padding=(10,0,10,20))
    forcing_grid = parent_layout[13, 1:2] = GridLayout()
    Label(forcing_grid[1, 2], "a", halign=:center)
    Label(forcing_grid[1, 3], "b", halign=:center)
    Label(forcing_grid[2, 1], "Left:", halign=:right)
    widgets[:forcing_left_a] = Textbox(forcing_grid[2, 2], placeholder="0.4", validator=Float64)
    widgets[:forcing_left_b] = Textbox(forcing_grid[2, 3], placeholder="0.4", validator=Float64)
    Label(forcing_grid[3, 1], "Right:", halign=:right)
    widgets[:forcing_right_a] = Textbox(forcing_grid[3, 2], placeholder="0.4", validator=Float64)
    widgets[:forcing_right_b] = Textbox(forcing_grid[3, 3], placeholder="0.4", validator=Float64)
    Label(forcing_grid[4, 1], "Bottom:", halign=:right)
    widgets[:forcing_bottom_a] = Textbox(forcing_grid[4, 2], placeholder="0.4", validator=Float64)
    widgets[:forcing_bottom_b] = Textbox(forcing_grid[4, 3], placeholder="0.4", validator=Float64)
    Label(forcing_grid[5, 1], "Top:", halign=:right)
    widgets[:forcing_top_a] = Textbox(forcing_grid[5, 2], placeholder="0.4", validator=Float64)
    widgets[:forcing_top_b] = Textbox(forcing_grid[5, 3], placeholder="0.4", validator=Float64)
    Label(forcing_grid[6, 1], " ", halign=:right)

    return widgets
end

function populate_button_panel!(parent_layout)
    btns = Dict{Symbol, Any}()
    
    btns[:split_domain] = Button(parent_layout[1, 1], label = "Split Domain")
    btns[:edge_solve] = Button(parent_layout[1, 2], label = "Solve Edges")
    btns[:smooth_grid] = Button(parent_layout[1, 3], label = "Smooth Grid")
    btns[:save_grid] = Button(parent_layout[1, 4], label = "Save Grid")
    btns[:reset_view] = Button(parent_layout[1, 5], label = "Reset View")
    Label(parent_layout[1, 6], " ", halign=:right, padding =(20, 20, 20, 20))
    return btns
end


function to_linesegments(X,Y)
    nx,ny=size(X)
    points=Point2f[]
    for j in 1:ny
        for i in 1:nx-1
            push!(points,Point2f(X[i,j],Y[i,j]),Point2f(X[i+1,j],Y[i+1,j]))
        end
    end
    for i in 1:nx, j in 1:ny-1
        push!(points,Point2f(X[i,j],Y[i,j]),Point2f(X[i,j+1],Y[i,j+1]))
    end
    return points
end


function split_domain(i_splits_str, j_splits_str)    
    i_splits_str = isnothing(i_splits_str) ? "" : i_splits_str
    j_splits_str = isnothing(j_splits_str) ? "" : j_splits_str

    # Parse the input strings into arrays of integers
    i_splits = isempty(i_splits_str) ? Int[] : parse.(Int, split(i_splits_str, ','))
    j_splits = isempty(j_splits_str) ? Int[] : parse.(Int, split(j_splits_str, ','))

    if i_splits == Int[] && j_splits == Int[]
        println("No splits provided.")
        return
    end
    println("Splitting domain with I-splits: [$(i_splits), $(j_splits)]")
    return
end







########################################
########################################
########################################
##########
# parameters
# metricFieldFile = "step/BFstepTest_entropy.metric"
# gridFolder = "step/coarseGrids"

# metricData,datatype,gridfile = GridGeneration.readTurtleFields(metricFieldFile)
# gridfile = joinpath(gridFolder, gridfile)
# blocks, centers, Xfa, interfaceInfo, bndInfo = GridGeneration.ImportTurtleGrid(gridfile)

# numBlocks = length(blocks)


# initialGrid = [blocks[1]]
# initialBndInfo = bndInfo
# initialInterfaceInfo = interfaceInfo

# tree, refs = GridGeneration.setup_metric_tree(centers)
# M = (x,y) -> GridGeneration.find_nearest_kd(metricData, tree, refs, x, y)

# simple example for now
s = 1000
M = (x,y) -> [s, s]

height = 2.0
width  = 4.0

# Define number of points
num_points_height = 50
num_points_width  = 100

top = hcat(range(0, stop=width, length=num_points_width), fill(height, num_points_width))
bottom = hcat(range(0, stop=width, length=num_points_width), fill(0.0, num_points_width))
right = hcat(fill(width, num_points_height), range(0, stop=height, length=num_points_height))
left = hcat(fill(0.0, num_points_height), range(0, stop=height, length=num_points_height))

initialGrid = GridGeneration.TFI([top, right, bottom, left])
initialBndInfo = []
push!(initialBndInfo, Dict("name"=>"bottom", "faceInfo" => [Dict("block" => 1, "start"=> [1,1,1], "end"=>[num_points_width,1,1])]))
push!(initialBndInfo, Dict("name"=>"right", "faceInfo" => [Dict("block" => 1, "start"=> [num_points_width,1,1], "end"=>[num_points_width,num_points_height,1])]))
push!(initialBndInfo, Dict("name"=>"top", "faceInfo" => [Dict("block" => 1, "start"=> [1,num_points_height,1], "end"=>[num_points_width,num_points_height,1])]))
push!(initialBndInfo, Dict("name"=>"left", "faceInfo" => [Dict("block" => 1, "start"=> [1,1,1], "end"=>[1,num_points_height,1])]))

initialInterfaceInfo = []

initialGrid = [initialGrid]
########################################
########################################
########################################



set_theme!(theme_light())


fig = Figure(size = (1600, 1200))

# --- Main Layout ---
right_layout = fig[1, 2] = GridLayout()
left_layout = fig[1, 1] = GridLayout(tellwidth = false, halign = :left, valign = :top)
colsize!(fig.layout, 1, Relative(1/5))
colsize!(fig.layout, 2, Relative(4/5))

control_panel = left_layout[1,1] = GridLayout()
plot_panel = right_layout[1,1] = GridLayout()
tool_bar = right_layout[2, 1] = GridLayout()
rowsize!(right_layout, 1, Relative(9/10))
colsize!(right_layout, 1, Relative(4/5))

Box(fig[1, 1], linestyle = :solid, strokecolor = :black, cornerradius = 20)
Box(right_layout[1, 1], linestyle = :solid, strokecolor = :black, cornerradius = 20)
Box(right_layout[2, 1], linestyle = :solid, strokecolor = :black, cornerradius = 20)



# --- Visualization Panels (Axes) ---
ax_initial = Axis(plot_panel[1, 1], aspect = DataAspect(), title = "Initial Grid (TFI)")
ax_generated = Axis(plot_panel[2, 1], aspect = DataAspect(), title = "Generated Grid")
ax_smoothed = Axis(plot_panel[3, 1], aspect = DataAspect(), title = "Smoothed Grid")

generated_grid = Observable{Any}(nothing)
generated_bndInfo = Observable{Any}(nothing)
generated_interfaceInfo = Observable{Any}(nothing)
smoothed_grid = Observable{Any}(nothing)
i_split_lines = Observable(Point2f[])
j_split_lines = Observable(Point2f[])

# set generated and smooth grid to initial for now
generated_grid[] = initialGrid[]
smoothed_grid[] = initialGrid[1]
generated_bndInfo[] = initialBndInfo
generated_interfaceInfo[] = initialInterfaceInfo


segments_initial = to_linesegments(initialGrid[1][1,:,:], initialGrid[1][2,:,:])
segments_generated = @lift(to_linesegments($generated_grid[1,:,:], $generated_grid[2,:,:]))

linesegments!(ax_initial, segments_initial, color = :black, linewidth=1.0)
linesegments!(ax_initial, i_split_lines, color = :red, linewidth = 1.5, overdraw=true)
linesegments!(ax_initial, i_split_lines, color = :red, linewidth = 1.5, overdraw=true)
linesegments!(ax_initial, j_split_lines, color = :green, linewidth = 1.5, overdraw=true)

linesegments!(ax_generated, segments_generated, color = :black, linewidth=1.0)

autolimits!(ax_initial)
autolimits!(ax_generated)







# --- Control Panel ---
controls = populate_control_panel!(control_panel)

# --- Tool Panel ---
tool = populate_button_panel!(tool_bar)



# --- Event Handling ---
on(tool[:split_domain].clicks) do _
    # Get the current values from the textboxes
    i_splits = controls[:i_splits].stored_string[]
    j_splits = controls[:j_splits].stored_string[]
    i_indices = i_splits === nothing ? Int[] : parse.(Int, split(i_splits, ','))
    j_indices = j_splits === nothing ? Int[] : parse.(Int, split(j_splits, ','))

    if i_indices == Int[] && j_indices == Int[]
        println("No splits provided.")
        return
    end

    # Call the placeholder function with the textbox values
    # split_domain(i_splits, j_splits)
    split_locations = [i_indices, j_indices]
    println("Splitting domain with I-splits: [$(i_indices), $(j_indices)]")
    println("Type of initialGrid: ", typeof(initialGrid))
    println("Type of initialBndInfo: ", typeof(initialBndInfo))
    println("Type of initialInterfaceInfo: ", typeof(initialInterfaceInfo))
    println("initialBndInfo: ", initialBndInfo)
    println("initialInterfaceInfo: ", initialInterfaceInfo)
    generatedBlocks, generatedBndInfo, generatedInterInfo = GridGeneration.SplitBlock(deepcopy(initialGrid[1]), split_locations, deepcopy(initialBndInfo), deepcopy(initialInterfaceInfo))
    # update the observables
    generated_grid[] = generatedBlocks[1]
    generated_bndInfo[] = generatedBndInfo
    generated_interfaceInfo[] = generatedInterInfo
end

on(tool[:edge_solve].clicks) do _
    edge_type = controls[:edge_solver].selection[]
    if edge_type == "analytic"
        edge_type = :analytic
    elseif edge_type == "numerical"
        edge_type = :numerical
    end

    # get the generated info
    generatedBlocks = generated_grid[]
    generatedBndInfo = generated_bndInfo[]
    generatedInterInfo = generated_interfaceInfo[]

    println(typeof(generatedBlocks))
    println(typeof(generatedBndInfo))
    println(typeof(generatedInterInfo))

    println("length of generatedBlocks: ", length(generatedBlocks))
    # println(generatedInterInfo)

    # SolvedBlocks, SolvedBndInfo, solvedInterInfo = 
    GridGeneration.SolveAllBlocks(M, deepcopy(generatedBlocks), deepcopy(generatedBndInfo), deepcopy(generatedInterInfo); solver=edge_type)
    # # update the current generated Info
    # generated_grid[] = SolvedBlocks[1]
    # generatedBndInfo[] = SolvedBndInfo  
    # generatedInterInfo[] = solvedInterInfo
end



on(tool[:reset_view].clicks) do _
    autolimits!(ax_initial)
    autolimits!(ax_generated)
    autolimits!(ax_smoothed)
end

# Event handlers for reactive split line highlighting
on(controls[:i_splits].stored_string) do s
    indices = s
    if indices == "none"
        i_split_lines[] = Point2f[]
        return
    end
    indices_list = parse.(Int, split(indices, ',')) 
    points = Point2f[]
    X, Y = initialGrid[1][1, :, :], initialGrid[1][2, :, :]
    _, ny = size(X)
    for i in indices_list
        # Ensure index is within bounds before trying to access
        if 1 <= i <= size(X, 1)
            for j in 1:(ny-1)
                p1 = Point2f(X[i, j], Y[i, j])
                p2 = Point2f(X[i, j+1], Y[i, j+1])
                push!(points, p1, p2)
            end
        end
    end
    i_split_lines[] = points
end

on(controls[:j_splits].stored_string) do s
    indices = s
    if indices == "none"
        j_split_lines[] = Point2f[]
        return
    end
    indices_list = parse.(Int, split(indices, ',')) 
    points = Point2f[]
    X, Y = initialGrid[1][1, :, :], initialGrid[1][2, :, :]
    nx, _ = size(X)
    for j in indices_list
        # Ensure index is within bounds before trying to access
        if 1 <= j <= size(X, 2)
            for i in 1:(nx-1)
                p1 = Point2f(X[i, j], Y[i, j])
                p2 = Point2f(X[i+1, j], Y[i+1, j])
                push!(points, p1, p2)
            end
        end
    end
    j_split_lines[] = points
end


fig