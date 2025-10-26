"""
Interactive GUI for GridGeneration.jl

This GUI provides:
- Visualization of initial grid
- Interactive split location selection (click on grid)
- Dropdown menus for solver and smoothing options
- Parameter controls via sliders and checkboxes
- Real-time grid generation and visualization

Requirements:
    using Pkg
    Pkg.add(["GLMakie", "DelimitedFiles", "MAT"])
"""

using GLMakie
using DelimitedFiles
using MAT: matread
using Printf

# Include GridGeneration module
include("../../src/GridGeneration.jl")
using .GridGeneration

# Include example setup functions
include("../airfoil/airfoil.jl")
include("../rectangle/rectangle.jl")

# ============================================================================
# Grid Visualization Functions
# ============================================================================

"""
Plot a 2D grid block on a given axis
"""
function plot_grid!(ax, block; color=:black, linewidth=1.0, show_points=false)
    X = block[1, :, :]
    Y = block[2, :, :]
    
    # Plot η-lines (lines of constant j)
    for j in 1:size(Y, 2)
        lines!(ax, X[:, j], Y[:, j], color=color, linewidth=linewidth)
    end
    
    # Plot ξ-lines (lines of constant i)
    for i in 1:size(X, 1)
        lines!(ax, X[i, :], Y[i, :], color=color, linewidth=linewidth)
    end
    
    if show_points
        scatter!(ax, vec(X), vec(Y), color=color, markersize=3)
    end
end

"""
Plot multiple blocks
"""
function plot_blocks!(ax, blocks; kwargs...)
    for block in blocks
        plot_grid!(ax, block; kwargs...)
    end
end

"""
Highlight grid lines at specified i and j indices
"""
function highlight_splits!(ax, block, i_splits, j_splits; color=:red, linewidth=3.0)
    X = block[1, :, :]
    Y = block[2, :, :]
    
    # Highlight i-splits (vertical lines)
    for i in i_splits
        if 1 <= i <= size(X, 1)
            lines!(ax, X[i, :], Y[i, :], color=color, linewidth=linewidth, linestyle=:dash)
        end
    end
    
    # Highlight j-splits (horizontal lines)
    for j in j_splits
        if 1 <= j <= size(Y, 2)
            lines!(ax, X[:, j], Y[:, j], color=color, linewidth=linewidth, linestyle=:dash)
        end
    end
end

# ============================================================================
# GUI State Management
# ============================================================================

mutable struct GUIState
    # Domain data
    case::Symbol
    initialGrid::Array{Float64, 3}
    bndInfo::Vector{Any}
    interInfo::Vector{Any}
    metric::Function
    
    # Split locations (indices on initial grid)
    i_splits::Vector{Int}
    j_splits::Vector{Int}
    
    # Parameters
    useSplitting::Bool
    useEdgeSolver::Bool
    boundarySolver::Symbol
    useSmoothing::Bool
    smoothMethod::Symbol
    
    # Elliptic parameters
    max_iter::Int
    tol::Float64
    ω::Float64
    useTopWall::Bool
    useBottomWall::Bool
    useLeftWall::Bool
    useRightWall::Bool
    a_decay::Float64
    b_decay::Float64
    verbose::Bool
    
    # Results
    resultBlocks::Union{Nothing, Vector{Array{Float64, 3}}}
    smoothBlocks::Union{Nothing, Vector{Array{Float64, 3}}}
    
    function GUIState()
        new(:airfoil, zeros(2,2,2), [], [], (x,y) -> [1.0, 1.0],
            Int[], Int[], 
            true, true, :analytic, true, :ellipticSS,
            5000, 1e-6, 0.2, true, true, true, true, 0.4, 0.4, false,
            nothing, nothing)
    end
end

# ============================================================================
# Main GUI Application
# ============================================================================

function launch_gui()
    # Create GUI state
    state = GUIState()
    
    @info "Initializing GUI components..."
    
    # Create figure and layout
    fig = Figure(resolution=(1800, 1000), backgroundcolor=:white)
    
    @info "Building control panel..."
    
    # ========================================================================
    # Left Panel - Controls
    # ========================================================================
    controls = fig[1, 1] = GridLayout()
    
    # Domain Selection
    Label(controls[1, 1], "Domain Type:", fontsize=16, halign=:left)
    domain_menu = Menu(controls[1, 2], 
        options=["Airfoil", "Rectangle"],
        default="Airfoil")
    
    # Load Domain Button
    load_btn = Button(controls[2, 1:2], label="Load Domain", fontsize=14)
    
    # Metric Selection
    Label(controls[3, 1], "Metric Problem:", fontsize=16, halign=:left)
    metric_slider = Slider(controls[3, 2], range=1:10, startvalue=6)
    metric_label = Label(controls[3, 3], "6", fontsize=14)
    
    Label(controls[4, 1], "Metric Scale:", fontsize=16, halign=:left)
    scale_slider = Slider(controls[4, 2], range=0.001:0.001:1.0, startvalue=0.05)
    scale_label = Label(controls[4, 3], "0.050", fontsize=14)
    
    # Splitting Controls
    Label(controls[5, 1:2], "═══ Splitting ═══", fontsize=18, halign=:center)
    
    use_splitting = Toggle(controls[6, 1], active=true)
    Label(controls[6, 2], "Use Block Splitting", fontsize=14, halign=:left)
    
    Label(controls[7, 1], "I-splits (vertical):", fontsize=14, halign=:left)
    i_splits_box = Textbox(controls[7, 2], placeholder="e.g., 300,400")
    
    Label(controls[8, 1], "J-splits (horizontal):", fontsize=14, halign=:left)
    j_splits_box = Textbox(controls[8, 2], placeholder="e.g., 30")
    
    # Solver Controls
    Label(controls[9, 1:2], "═══ Solvers ═══", fontsize=18, halign=:center)
    
    use_edge_solver = Toggle(controls[10, 1], active=true)
    Label(controls[10, 2], "Use Edge Solver", fontsize=14, halign=:left)
    
    Label(controls[11, 1], "Boundary Solver:", fontsize=14, halign=:left)
    boundary_solver_menu = Menu(controls[11, 2],
        options=["analytic", "numeric"],
        default="analytic")
    
    # Smoothing Controls
    Label(controls[12, 1:2], "═══ Smoothing ═══", fontsize=18, halign=:center)
    
    use_smoothing = Toggle(controls[13, 1], active=true)
    Label(controls[13, 2], "Use Smoothing", fontsize=14, halign=:left)
    
    Label(controls[14, 1], "Smooth Method:", fontsize=14, halign=:left)
    smooth_method_menu = Menu(controls[14, 2],
        options=["ellipticSS"],
        default="ellipticSS")
    
    # Elliptic Parameters
    Label(controls[15, 1:2], "═══ Elliptic Parameters ═══", fontsize=18, halign=:center)
    
    Label(controls[16, 1], "Max Iterations:", fontsize=14, halign=:left)
    max_iter_slider = Slider(controls[16, 2], range=100:100:10000, startvalue=5000)
    max_iter_label = Label(controls[16, 3], "5000", fontsize=14)
    
    Label(controls[17, 1], "Tolerance:", fontsize=14, halign=:left)
    tol_slider = Slider(controls[17, 2], range=1e-8:1e-8:1e-4, startvalue=1e-6)
    tol_label = Label(controls[17, 3], "1e-6", fontsize=14)
    
    Label(controls[18, 1], "Relaxation ω:", fontsize=14, halign=:left)
    omega_slider = Slider(controls[18, 2], range=0.1:0.01:2.0, startvalue=0.2)
    omega_label = Label(controls[18, 3], "0.20", fontsize=14)
    
    Label(controls[19, 1], "Decay Parameter:", fontsize=14, halign=:left)
    decay_slider = Slider(controls[19, 2], range=0.1:0.05:1.0, startvalue=0.4)
    decay_label = Label(controls[19, 3], "0.40", fontsize=14)
    
    # Wall forcing toggles
    use_top_wall = Toggle(controls[20, 1], active=true)
    Label(controls[20, 2], "Top Wall Forcing", fontsize=14, halign=:left)
    
    use_bottom_wall = Toggle(controls[21, 1], active=true)
    Label(controls[21, 2], "Bottom Wall Forcing", fontsize=14, halign=:left)
    
    use_left_wall = Toggle(controls[22, 1], active=true)
    Label(controls[22, 2], "Left Wall Forcing", fontsize=14, halign=:left)
    
    use_right_wall = Toggle(controls[23, 1], active=true)
    Label(controls[23, 2], "Right Wall Forcing", fontsize=14, halign=:left)
    
    verbose_toggle = Toggle(controls[24, 1], active=false)
    Label(controls[24, 2], "Verbose Output", fontsize=14, halign=:left)
    
    # Generate Button
    generate_btn = Button(controls[25, 1:2], label="Generate Grid", 
        fontsize=16, buttoncolor=:green)
    
    # Status Label
    status_label = Label(controls[26, 1:2], "Ready", fontsize=14, 
        halign=:center, color=:blue)
    
    # ========================================================================
    # Right Panel - Visualization
    # ========================================================================
    viz_layout = fig[1, 2] = GridLayout()
    
    # Initial Grid View
    Label(viz_layout[1, 1], "Initial Grid", fontsize=18, halign=:center)
    ax_initial = Axis(viz_layout[2, 1], aspect=DataAspect())
    
    # Result Grid View
    Label(viz_layout[1, 2], "Generated Grid", fontsize=18, halign=:center)
    ax_result = Axis(viz_layout[2, 2], aspect=DataAspect())
    
    # Smoothed Grid View (if applicable)
    Label(viz_layout[3, 1:2], "Smoothed Grid", fontsize=18, halign=:center)
    ax_smooth = Axis(viz_layout[4, 1:2], aspect=DataAspect())
    
    # ========================================================================
    # Event Handlers
    # ========================================================================
    
    # Update metric labels
    on(metric_slider.value) do val
        metric_label.text[] = string(val)
    end
    
    on(scale_slider.value) do val
        scale_label.text[] = @sprintf("%.3f", val)
    end
    
    on(max_iter_slider.value) do val
        max_iter_label.text[] = string(val)
    end
    
    on(tol_slider.value) do val
        tol_label.text[] = @sprintf("%.1e", val)
    end
    
    on(omega_slider.value) do val
        omega_label.text[] = @sprintf("%.2f", val)
    end
    
    on(decay_slider.value) do val
        decay_label.text[] = @sprintf("%.2f", val)
    end
    
    # Load Domain Handler
    on(load_btn.clicks) do n
        try
            status_label.text[] = "Loading domain..."
            status_label.color[] = :orange
            
            # Determine case
            state.case = domain_menu.selection[] == "Airfoil" ? :airfoil : :rectangle
            
            # Load domain
            if state.case == :airfoil
                state.initialGrid, state.bndInfo, state.interInfo = 
                    GetAirfoilSetup(
                        airfoilPath="examples/airfoil/data/A-airfoil.txt", 
                        radius=3, 
                        type=:cgrid
                    )
                # Get metric
                problem = metric_slider.value[]
                state.metric = GetAirfoilMetric(problem; scale=scale_slider.value[])
            else
                initDomain = GetRectangleDomain()
                state.initialGrid = TFI(initDomain)
                state.bndInfo = Any[]
                state.interInfo = Any[]
                # Get metric
                problem = metric_slider.value[]
                state.metric = GetRectangleMetric(problem; scale=10000)
            end
            
            # Clear axes and plot initial grid
            empty!(ax_initial)
            plot_grid!(ax_initial, state.initialGrid, color=:black, linewidth=0.5)
            autolimits!(ax_initial)
            
            status_label.text[] = "Domain loaded successfully"
            status_label.color[] = :green
        catch e
            status_label.text[] = "Error loading domain: $(e)"
            status_label.color[] = :red
            @error "Domain loading failed" exception=(e, catch_backtrace())
        end
    end
    
    # Parse split locations from text boxes
    function parse_splits(text)
        if isempty(strip(text))
            return Int[]
        end
        try
            return parse.(Int, split(strip(text), ','))
        catch
            return Int[]
        end
    end
    
    # Update split visualization when text changes
    on(i_splits_box.stored_string) do text
        state.i_splits = parse_splits(text)
        if !isnothing(state.initialGrid) && size(state.initialGrid) != (2,2,2)
            empty!(ax_initial)
            plot_grid!(ax_initial, state.initialGrid, color=:black, linewidth=0.5)
            highlight_splits!(ax_initial, state.initialGrid, state.i_splits, state.j_splits)
            autolimits!(ax_initial)
        end
    end
    
    on(j_splits_box.stored_string) do text
        state.j_splits = parse_splits(text)
        if !isnothing(state.initialGrid) && size(state.initialGrid) != (2,2,2)
            empty!(ax_initial)
            plot_grid!(ax_initial, state.initialGrid, color=:black, linewidth=0.5)
            highlight_splits!(ax_initial, state.initialGrid, state.i_splits, state.j_splits)
            autolimits!(ax_initial)
        end
    end
    
    # Generate Grid Handler
    on(generate_btn.clicks) do n
        try
            status_label.text[] = "Generating grid..."
            status_label.color[] = :orange
            
            # Build parameters
            elliptic_params = EllipticParams(
                max_iter=max_iter_slider.value[],
                tol=tol_slider.value[],
                ω=omega_slider.value[],
                useTopWall=use_top_wall.active[],
                useBottomWall=use_bottom_wall.active[],
                useLeftWall=use_left_wall.active[],
                useRightWall=use_right_wall.active[],
                a_decay_left=decay_slider.value[],
                b_decay_left=decay_slider.value[],
                a_decay_right=decay_slider.value[],
                b_decay_right=decay_slider.value[],
                a_decay_top=decay_slider.value[],
                b_decay_top=decay_slider.value[],
                a_decay_bottom=decay_slider.value[],
                b_decay_bottom=decay_slider.value[],
                verbose=verbose_toggle.active[]
            )
            
            splitLocs = [state.i_splits, state.j_splits]
            
            params = SimParams(
                useSplitting=use_splitting.active[],
                splitLocations=splitLocs,
                useEdgeSolver=use_edge_solver.active[],
                boundarySolver=Symbol(boundary_solver_menu.selection[]),
                useSmoothing=use_smoothing.active[],
                smoothMethod=Symbol(smooth_method_menu.selection[]),
                elliptic=elliptic_params
            )
            
            # Generate grid
            result = GenerateGrid(
                state.initialGrid, 
                state.bndInfo, 
                state.interInfo, 
                state.metric, 
                params=params
            )
            
            # Parse results based on what was returned
            if length(result) == 6
                state.smoothBlocks, state.resultBlocks, _, _, _, _ = result
            elseif length(result) == 3
                state.resultBlocks, _, _ = result
                state.smoothBlocks = nothing
            end
            
            # Visualize results
            empty!(ax_result)
            if !isnothing(state.resultBlocks)
                plot_blocks!(ax_result, state.resultBlocks, color=:blue, linewidth=0.5)
                autolimits!(ax_result)
            end
            
            empty!(ax_smooth)
            if !isnothing(state.smoothBlocks)
                plot_blocks!(ax_smooth, state.smoothBlocks, color=:darkgreen, linewidth=0.5)
                autolimits!(ax_smooth)
            end
            
            status_label.text[] = "Grid generated successfully!"
            status_label.color[] = :green
            
        catch e
            status_label.text[] = "Error generating grid: $(e)"
            status_label.color[] = :red
            @error "Grid generation failed" exception=(e, catch_backtrace())
        end
    end
    
    # ========================================================================
    # Display
    # ========================================================================
    
    @info "Opening GUI window..."
    
    # Display the figure and get the screen handle
    screen = display(fig)
    
    @info "GUI ready! Window should now be visible."
    @info "To close: Press Ctrl+C in terminal or close the window."
    
    return fig, state, screen
end

# ============================================================================
# Launch Application
# ============================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    @info "Launching GridGeneration GUI..."
    fig, state, screen = launch_gui()
    
    # Keep the window open by waiting for it to close
    @info "GUI window opened. Close the window or press Ctrl+C to exit."
    
    try
        # Wait for the window to be closed
        # This checks if the screen/window is still open
        while isopen(screen)
            sleep(0.1)  # Small sleep to prevent busy waiting
        end
        @info "GUI window closed."
    catch e
        if isa(e, InterruptException)
            @info "GUI closed by user (Ctrl+C)."
        else
            rethrow(e)
        end
    end
end
