# It's assumed that your own GridGeneration library is in the path.
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
    
    # highlight_layout = parent_layout[2, 1] = GridLayout()
    btns[:highlight_boundaries] = Toggle(parent_layout[1, 6], active = false)
    Label(parent_layout[1, 7], "Highlight Boundaries", halign=:left, padding=(0,0,0,0))

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
    plot_blocks_with_highlighting!(ax, blocks_obs, highlight_obs)

Helper function to create a reactive plot of grid blocks.
Toggles boundary highlighting based on an observable.
"""
function plot_blocks_with_highlighting!(ax::Axis, blocks_obs::Observable, highlight_obs::Observable{Any})
    segments_obs = Observable(Point2f[])
    colors_obs = Observable(RGBAf[])

    onany(blocks_obs, highlight_obs) do blocks, highlight_on
        all_points = Point2f[]
        all_colors = RGBAf[]
        default_color = RGBAf(0,0,0,1)
        highlight_color = RGBAf(1,0,1,1) # Magenta

        !isnothing(blocks) && isa(blocks, Vector) || return

        for block in blocks
            X, Y = block[1,:,:], block[2,:,:]
            nx, ny = size(X)
            
            # Horizontal lines
            for j in 1:ny
                is_boundary = (j == 1 || j == ny)
                color = (highlight_on && is_boundary) ? highlight_color : default_color
                for i in 1:(nx-1)
                    push!(all_points, Point2f(X[i,j], Y[i,j]), Point2f(X[i+1,j], Y[i+1,j]))
                    push!(all_colors, color)
                end
            end
            
            # Vertical lines
            for i in 1:nx
                is_boundary = (i == 1 || i == nx)
                color = (highlight_on && is_boundary) ? highlight_color : default_color
                for j in 1:(ny-1)
                    push!(all_points, Point2f(X[i,j], Y[i,j]), Point2f(X[i,j+1], Y[i,j+1]))
                    push!(all_colors, color)
                end
            end
        end
        
        segments_obs[] = all_points
        colors_obs[] = all_colors
        autolimits!(ax)
    end
    
    linesegments!(ax, segments_obs, color = colors_obs, linewidth=1.0)
    notify(blocks_obs) # Trigger initial plot
end

########################################
# --- Data Loading ---
########################################
s = 1000
M = (x,y) -> [s * x, s]
height = 2.0; width  = 4.0
num_points_height = 50; num_points_width  = 100
top = hcat(range(0, stop=width, length=num_points_width), fill(height, num_points_width))
bottom = hcat(range(0, stop=width, length=num_points_width), fill(0.0, num_points_width))
right = hcat(fill(width, num_points_height), range(0, stop=height, length=num_points_height))
left = hcat(fill(0.0, num_points_height), range(0, stop=height, length=num_points_height))
initial_block = GridGeneration.TFI([top, right, bottom, left])
initialBndInfo = []
push!(initialBndInfo, Dict("name"=>"bottom", "faceInfo" => [Dict("block" => 1, "start"=> [1,1,1], "end"=>[num_points_width,1,1])]))
push!(initialBndInfo, Dict("name"=>"right", "faceInfo" => [Dict("block" => 1, "start"=> [num_points_width,1,1], "end"=>[num_points_width,num_points_height,1])]))
push!(initialBndInfo, Dict("name"=>"top", "faceInfo" => [Dict("block" => 1, "start"=> [1,num_points_height,1], "end"=>[num_points_width,num_points_height,1])]))
push!(initialBndInfo, Dict("name"=>"left", "faceInfo" => [Dict("block" => 1, "start"=> [1,1,1], "end"=>[1,num_points_height,1])]))
initialInterfaceInfo = []
initialGrid = [initial_block]


# metricFieldFile = "step/BFstepTest_entropy.metric"
# gridFolder = "step/coarseGrids"

# metricData,datatype,gridfile = GridGeneration.readTurtleFields(metricFieldFile)
# gridfile = joinpath(gridFolder, gridfile)
# blocks, centers, Xfa, interfaceInfo, bndInfo = GridGeneration.ImportTurtleGrid(gridfile)

# numBlocks = length(blocks)


# initialGrid = blocks
# initialBndInfo = bndInfo
# initialInterfaceInfo = interfaceInfo

# tree, refs = GridGeneration.setup_metric_tree(centers)
# M = (x,y) -> GridGeneration.find_nearest_kd(metricData, tree, refs, x, y)

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
rowsize!(right_layout, 2, Relative(1/10))
Box(fig[1, 1], linestyle = :solid, strokecolor = :black, cornerradius = 20)
Box(right_layout[1, 1], linestyle = :solid, strokecolor = :black, cornerradius = 20)
Box(right_layout[2, 1], linestyle = :solid, strokecolor = :black, cornerradius = 20)

# --- Visualization Panels (Axes) ---
ax_initial = Axis(plot_panel[1, 1], aspect = DataAspect(), title = "Initial Grid (TFI)")
ax_generated = Axis(plot_panel[2, 1], aspect = DataAspect(), title = "Generated Grid")
ax_smoothed = Axis(plot_panel[3, 1], aspect = DataAspect(), title = "Smoothed Grid")

# --- Observables ---
initial_blocks_obs = Observable(deepcopy(initialGrid))
generated_blocks = Observable(deepcopy(initialGrid))
generated_bndInfo = Observable(deepcopy(initialBndInfo))
generated_interfaceInfo = Observable(deepcopy(initialInterfaceInfo))
smoothed_blocks = Observable(deepcopy(initialGrid))
i_split_lines = Observable(Point2f[])
j_split_lines = Observable(Point2f[])

# --- Control Panel & Tool Panel Setup ---
controls = populate_control_panel!(control_panel)
tool = populate_button_panel!(tool_bar)
highlight_toggle = tool[:highlight_boundaries]

# --- Plotting ---
plot_blocks_with_highlighting!(ax_initial, initial_blocks_obs, highlight_toggle.active)
plot_blocks_with_highlighting!(ax_generated, generated_blocks, highlight_toggle.active)
plot_blocks_with_highlighting!(ax_smoothed, smoothed_blocks, highlight_toggle.active)

# Add split lines to the initial plot
linesegments!(ax_initial, i_split_lines, color = :red, linewidth = 1.5, overdraw=true)
linesegments!(ax_initial, j_split_lines, color = :green, linewidth = 1.5, overdraw=true)

# --- Event Handling ---
on(tool[:split_domain].clicks) do _
    i_splits_str = controls[:i_splits].stored_string[]
    j_splits_str = controls[:j_splits].stored_string[]
    i_indices = parse_indices(i_splits_str)
    j_indices = parse_indices(j_splits_str)

    if isempty(i_indices) && isempty(j_indices); println("No splits provided."); return; end

    split_locations = [i_indices, j_indices]
    println("Splitting domain...")
    new_blocks, new_bnd_info, new_inter_info = GridGeneration.SplitBlock(deepcopy(initialGrid[1]), split_locations, deepcopy(initialBndInfo), deepcopy(initialInterfaceInfo))
    generated_blocks[] = new_blocks
    generated_bndInfo[] = new_bnd_info
    generated_interfaceInfo[] = new_inter_info
end

on(tool[:edge_solve].clicks) do _
    edge_type = controls[:edge_solver].selection[]
    solver_sym = edge_type == "analytic" ? :analytic : :numerical
    current_blocks = generated_blocks[]
    current_bnd_info = generated_bndInfo[]
    current_inter_info = generated_interfaceInfo[]
    println("Solving edges for $(length(current_blocks)) block(s)...")
    SolvedBlocks, SolvedBndInfo, solvedInterInfo = GridGeneration.SolveAllBlocks(M, deepcopy(current_blocks), deepcopy(current_bnd_info), deepcopy(current_inter_info); solver=solver_sym)
    # SolvedBlocks, SolvedBndInfo, solvedInterInfo = GridGeneration.SolveAllBlocks(M, deepcopy(current_blocks), [], []; solver=solver_sym)
    generated_blocks[] = SolvedBlocks
    generated_bndInfo[] = SolvedBndInfo
    generated_interfaceInfo[] = solvedInterInfo
    println("Edge solving complete.")
end

on(tool[:reset_view].clicks) do _
    autolimits!(ax_initial)
    autolimits!(ax_generated)
    autolimits!(ax_smoothed)
end

on(controls[:i_splits].stored_string) do s
    indices = parse_indices(s)
    points = Point2f[]
    if !isempty(indices)
        X, Y = initialGrid[1][1, :, :], initialGrid[1][2, :, :]
        _, ny = size(X)
        for i in indices; if 1 <= i <= size(X, 1); for j in 1:(ny-1)
            push!(points, Point2f(X[i, j], Y[i, j]), Point2f(X[i, j+1], Y[i, j+1]))
        end; end; end
    end
    i_split_lines[] = points
end

on(controls[:j_splits].stored_string) do s
    indices = parse_indices(s)
    points = Point2f[]
    if !isempty(indices)
        X, Y = initialGrid[1][1, :, :], initialGrid[1][2, :, :]
        nx, _ = size(X)
        for j in indices; if 1 <= j <= size(X, 2); for i in 1:(nx-1)
            push!(points, Point2f(X[i, j], Y[i, j]), Point2f(X[i+1, j], Y[i+1, j]))
        end; end; end
    end
    j_split_lines[] = points
end

fig

