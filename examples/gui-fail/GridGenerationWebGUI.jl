"""
Web-based Interactive GUI for GridGeneration.jl using Interact.jl and Blink.jl

This provides a browser-based interface with:
- Domain selection and loading
- Interactive parameter controls
- Real-time visualization using Plots.jl
- Simpler setup than GLMakie version

Requirements:
    using Pkg
    Pkg.add(["Interact", "Blink", "Plots", "WebIO"])
"""

using Interact
using Blink
using Plots
using DelimitedFiles
using MAT: matread

# Include GridGeneration module
include("../../src/GridGeneration.jl")
using .GridGeneration

# Include example setup functions
include("../airfoil/airfoil.jl")
include("../rectangle/rectangle.jl")

# ============================================================================
# Plotting Functions
# ============================================================================

"""
Plot a single grid block
"""
function plot_grid(block; kwargs...)
    p = plot(aspect_ratio=:equal, legend=false, grid=false, framestyle=:box)
    X = block[1, :, :]
    Y = block[2, :, :]
    
    # Plot η-lines (constant j)
    for j in 1:size(Y, 2)
        plot!(p, X[:, j], Y[:, j]; kwargs...)
    end
    
    # Plot ξ-lines (constant i)
    for i in 1:size(X, 1)
        plot!(p, X[i, :], Y[i, :]; kwargs...)
    end
    
    return p
end

"""
Plot multiple blocks
"""
function plot_blocks(blocks; kwargs...)
    p = plot(aspect_ratio=:equal, legend=false, grid=false, framestyle=:box)
    
    for block in blocks
        X = block[1, :, :]
        Y = block[2, :, :]
        
        for j in 1:size(Y, 2)
            plot!(p, X[:, j], Y[:, j]; kwargs...)
        end
        
        for i in 1:size(X, 1)
            plot!(p, X[i, :], Y[i, :]; kwargs...)
        end
    end
    
    return p
end

# ============================================================================
# Interactive GUI
# ============================================================================

function create_web_gui()
    
    # ========================================================================
    # Control Widgets
    # ========================================================================
    
    # Domain Selection
    domain_select = dropdown(
        ["Airfoil", "Rectangle"],
        label="Domain Type:",
        value="Airfoil"
    )
    
    metric_problem = slider(1:10, label="Metric Problem:", value=6)
    metric_scale = slider(0.001:0.001:1.0, label="Metric Scale:", value=0.05)
    
    load_button = button("Load Domain")
    
    # Splitting Parameters
    use_splitting = checkbox(true, label="Use Block Splitting")
    i_splits_input = textbox("300,400", label="I-splits (comma-separated):")
    j_splits_input = textbox("30", label="J-splits (comma-separated):")
    
    # Solver Parameters
    use_edge_solver = checkbox(true, label="Use Edge Solver")
    boundary_solver = dropdown(
        ["analytic", "numeric"],
        label="Boundary Solver:",
        value="analytic"
    )
    
    # Smoothing Parameters
    use_smoothing = checkbox(true, label="Use Smoothing")
    smooth_method = dropdown(
        ["ellipticSS"],
        label="Smooth Method:",
        value="ellipticSS"
    )
    
    # Elliptic Parameters
    max_iter = slider(100:100:10000, label="Max Iterations:", value=5000)
    tolerance = slider([1e-8, 1e-7, 1e-6, 1e-5, 1e-4], label="Tolerance:", value=1e-6)
    omega = slider(0.1:0.01:2.0, label="Relaxation ω:", value=0.2)
    decay_param = slider(0.1:0.05:1.0, label="Decay Parameter:", value=0.4)
    
    use_top_wall = checkbox(true, label="Top Wall Forcing")
    use_bottom_wall = checkbox(true, label="Bottom Wall Forcing")
    use_left_wall = checkbox(true, label="Left Wall Forcing")
    use_right_wall = checkbox(true, label="Right Wall Forcing")
    verbose = checkbox(false, label="Verbose Output")
    
    generate_button = button("Generate Grid")
    
    # Status display
    status = Observable("Ready to load domain")
    
    # ========================================================================
    # State Management
    # ========================================================================
    
    initial_grid = Observable{Union{Nothing, Array{Float64, 3}}}(nothing)
    bnd_info = Observable{Vector{Any}}([])
    inter_info = Observable{Vector{Any}}([])
    metric_func = Observable{Function}((x,y) -> [1.0, 1.0])
    
    result_blocks = Observable{Union{Nothing, Vector{Array{Float64, 3}}}}(nothing)
    smooth_blocks = Observable{Union{Nothing, Vector{Array{Float64, 3}}}}(nothing)
    
    # ========================================================================
    # Reactive Plots
    # ========================================================================
    
    # Initial grid plot
    initial_plot = map(initial_grid) do grid
        if isnothing(grid) || size(grid) == (2, 2, 2)
            plot(title="Initial Grid - Not Loaded", legend=false)
        else
            plot_grid(grid, color=:black, linewidth=0.5, title="Initial Grid")
        end
    end
    
    # Result plot
    result_plot = map(result_blocks) do blocks
        if isnothing(blocks)
            plot(title="Generated Grid - Not Generated", legend=false)
        else
            plot_blocks(blocks, color=:blue, linewidth=0.5, title="Generated Grid")
        end
    end
    
    # Smoothed plot
    smooth_plot = map(smooth_blocks) do blocks
        if isnothing(blocks)
            plot(title="Smoothed Grid - Not Generated", legend=false)
        else
            plot_blocks(blocks, color=:darkgreen, linewidth=0.5, title="Smoothed Grid")
        end
    end
    
    # ========================================================================
    # Event Handlers
    # ========================================================================
    
    # Load domain
    on(observe(load_button)) do _
        try
            status[] = "Loading domain..."
            
            case_type = domain_select[] == "Airfoil" ? :airfoil : :rectangle
            
            if case_type == :airfoil
                grid, bnd, inter = GetAirfoilSetup(
                    airfoilPath="examples/airfoil/data/A-airfoil.txt",
                    radius=3,
                    type=:cgrid
                )
                initial_grid[] = grid
                bnd_info[] = bnd
                inter_info[] = inter
                metric_func[] = GetAirfoilMetric(metric_problem[]; scale=metric_scale[])
            else
                domain = GetRectangleDomain()
                initial_grid[] = TFI(domain)
                bnd_info[] = []
                inter_info[] = []
                metric_func[] = GetRectangleMetric(metric_problem[]; scale=10000)
            end
            
            status[] = "Domain loaded successfully!"
        catch e
            status[] = "Error loading domain: $e"
            @error "Failed to load domain" exception=(e, catch_backtrace())
        end
    end
    
    # Generate grid
    on(observe(generate_button)) do _
        try
            if isnothing(initial_grid[])
                status[] = "Please load a domain first!"
                return
            end
            
            status[] = "Generating grid..."
            
            # Parse split locations
            i_splits = if isempty(strip(i_splits_input[]))
                Int[]
            else
                parse.(Int, split(strip(i_splits_input[]), ','))
            end
            
            j_splits = if isempty(strip(j_splits_input[]))
                Int[]
            else
                parse.(Int, split(strip(j_splits_input[]), ','))
            end
            
            # Build parameters
            elliptic_params = EllipticParams(
                max_iter=max_iter[],
                tol=tolerance[],
                ω=omega[],
                useTopWall=use_top_wall[],
                useBottomWall=use_bottom_wall[],
                useLeftWall=use_left_wall[],
                useRightWall=use_right_wall[],
                a_decay_left=decay_param[],
                b_decay_left=decay_param[],
                a_decay_right=decay_param[],
                b_decay_right=decay_param[],
                a_decay_top=decay_param[],
                b_decay_top=decay_param[],
                a_decay_bottom=decay_param[],
                b_decay_bottom=decay_param[],
                verbose=verbose[]
            )
            
            params = SimParams(
                useSplitting=use_splitting[],
                splitLocations=[i_splits, j_splits],
                useEdgeSolver=use_edge_solver[],
                boundarySolver=Symbol(boundary_solver[]),
                useSmoothing=use_smoothing[],
                smoothMethod=Symbol(smooth_method[]),
                elliptic=elliptic_params
            )
            
            # Generate
            result = GenerateGrid(
                initial_grid[],
                bnd_info[],
                inter_info[],
                metric_func[],
                params=params
            )
            
            # Update results
            if length(result) == 6
                smooth_blocks[], result_blocks[], _, _, _, _ = result
            elseif length(result) == 3
                result_blocks[], _, _ = result
                smooth_blocks[] = nothing
            end
            
            status[] = "Grid generated successfully!"
        catch e
            status[] = "Error generating grid: $e"
            @error "Failed to generate grid" exception=(e, catch_backtrace())
        end
    end
    
    # ========================================================================
    # Layout
    # ========================================================================
    
    ui = vbox(
        hbox(
            vbox(  # Left panel - Controls
                html("<h2>GridGeneration Controls</h2>"),
                html("<h3>Domain Setup</h3>"),
                domain_select,
                metric_problem,
                metric_scale,
                load_button,
                html("<hr>"),
                html("<h3>Splitting</h3>"),
                use_splitting,
                i_splits_input,
                j_splits_input,
                html("<hr>"),
                html("<h3>Solvers</h3>"),
                use_edge_solver,
                boundary_solver,
                html("<hr>"),
                html("<h3>Smoothing</h3>"),
                use_smoothing,
                smooth_method,
                html("<hr>"),
                html("<h3>Elliptic Parameters</h3>"),
                max_iter,
                tolerance,
                omega,
                decay_param,
                use_top_wall,
                use_bottom_wall,
                use_left_wall,
                use_right_wall,
                verbose,
                html("<hr>"),
                generate_button,
                html("<hr>"),
                html("<div style='padding: 10px; background-color: #e8f4f8; border-radius: 5px;'>"),
                html("<strong>Status:</strong> "),
                map(x -> html("<span>$x</span>"), status),
                html("</div>")
            ),
            vbox(  # Right panel - Visualizations
                html("<h2>Grid Visualizations</h2>"),
                hbox(initial_plot, result_plot),
                smooth_plot
            )
        )
    )
    
    return ui
end

# ============================================================================
# Launch Application
# ============================================================================

function launch_web_gui()
    @info "Launching web-based GridGeneration GUI..."
    
    ui = create_web_gui()
    
    # Open in Blink window
    w = Window()
    body!(w, ui)
    
    return w
end

# Launch if running as script
if abspath(PROGRAM_FILE) == @__FILE__
    window = launch_web_gui()
end
