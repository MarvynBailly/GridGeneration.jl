include("../../src/GridGeneration.jl")
using .GridGeneration

using GLMakie
using Makie: Polygon
using DelimitedFiles

set_theme!(theme_light())

fig = Figure(size = (1600, 1200))

# --- Main Layout ---
# Plotting area on the left, control panel on the right
vis_layout = fig[1, 2] = GridLayout()
# control_panel = fig[1, 1] = GridLayout(tellwidth = false, width = 350, halign = :center, valign = :top)
control_panel = fig[1, 1] = GridLayout(tellwidth = false, halign = :left, valign = :top)
colsize!(fig.layout, 1, Relative(1/3))

Box(fig[1,1], linestyle = :solid, strokecolor = :black, cornerradius = 20)


# --- Visualization Panels (Axes) ---
# ax_initial = Axis(vis_layout[1, 1], aspect = DataAspect(), title = "Initial Grid (TFI)")
# ax_generated = Axis(vis_layout[2, 1], aspect = DataAspect(), title = "Generated Grid")
# ax_smoothed = Axis(vis_layout[3, 1], aspect = DataAspect(), title = "Smoothed Grid")


# --- Observables for Grid Data ---
# These will hold the (X, Y) tuples for each grid state
initial_grid = Observable{Any}(nothing)
generated_grid = Observable{Any}(nothing)
smoothed_grid = Observable{Any}(nothing)

# Splits Configuration
Label(control_panel[1, 1:2], "Splits Locations", font = :bold, halign = :left, padding=(10,0,10,20), valign=:top)
Label(control_panel[2, 1], "I-splits (e.g., 10,25)", halign=:left, valign=:top, padding=(10,0,0,0))
Label(control_panel[2, 2], "J-splits (e.g., 8)")
i_splits_box = Textbox(control_panel[3, 1], placeholder="none")
j_splits_box = Textbox(control_panel[3, 2], placeholder="none")

split_domain_button = Button(control_panel[4, 1], label = "Split Domain")


# Solver & Smoothing Parameters
Label(control_panel[5, 1:2], "Boundary Solver Parameters", font = :bold, halign = :left, padding=(10,0,10,20), valign=:top)

# Solver Settings
Label(control_panel[6, 1], "Edge Solver Type:", halign=:left, padding=(10,0,0,0))
edge_solver = Menu(control_panel[6, 2], options = ["analytic", "numerical"], default = "analytic")
edge_solve_button = Button(control_panel[7, 1], label = "Solve Edges")


# Smoothing Settings
Label(control_panel[8, 1:2], "Smoothing Parameters", font = :bold, halign = :left, padding=(10,0,10,20), valign=:top)

# Solver Settings
Label(control_panel[9, 1], "Smoothing Type:", halign=:left, padding =(10,0,0,0))
smoothing_type = Menu(control_panel[9, 2], options = ["Elliptic-SS"], default = "Elliptic-SS")

Label(control_panel[10, 1], "Max Iterations:", padding=(10,0,0,0))
max_iter_box = Textbox(control_panel[10, 2], placeholder="5000", validator=Int)
Label(control_panel[11, 1], "Tolerance:", padding=(10,0,0,0))
tolerance_box = Textbox(control_panel[11, 2], placeholder="1e-5", validator=Float64)

Label(control_panel[12, 1], ("Omega: "), halign=:left, padding=(10,0,0,0))
omega_box = Textbox(control_panel[12, 2], placeholder="0.2", validator=Float64)

Label(control_panel[15, 1:2], "Wall Forcing Parameters", font = :bold, halign = :left, padding=(0,0,10,0))
forcing_grid = control_panel[16, 1:2] = GridLayout()
Label(forcing_grid[1, 2], "a", halign=:center)
Label(forcing_grid[1, 3], "b", halign=:center)
Label(forcing_grid[2, 1], "Left:", halign=:right)
wall_forcing_left_a = Textbox(forcing_grid[2, 2], placeholder="0.4", validator=Float64)
wall_forcing_left_b = Textbox(forcing_grid[2, 3], placeholder="0.4", validator=Float64)
Label(forcing_grid[3, 1], "Right:", halign=:right)
wall_forcing_right_a = Textbox(forcing_grid[3, 2], placeholder="0.4", validator=Float64)
wall_forcing_right_b = Textbox(forcing_grid[3, 3], placeholder="0.4", validator=Float64)
Label(forcing_grid[4, 1], "Bottom:", halign=:right)
wall_forcing_bottom_a = Textbox(forcing_grid[4, 2], placeholder="0.4", validator=Float64)
wall_forcing_bottom_b = Textbox(forcing_grid[4, 3], placeholder="0.4", validator=Float64)
Label(forcing_grid[5, 1], "Top:", halign=:right)
wall_forcing_top_a = Textbox(forcing_grid[5, 2], placeholder="0.4", validator=Float64)
wall_forcing_top_b = Textbox(forcing_grid[5, 3], placeholder="0.4", validator=Float64)




# # Section 1: Grid Dimensions
# Label(control_panel[1, 1:2], "Grid Dimensions", font = :bold, tellwidth=false, halign=:left)
# sg_dims = SliderGrid(
#     control_panel[2, 1:2],
#     (label = "Grid Points Nx", range = 5:2:100, startvalue = 20),
#     (label = "Grid Points Ny", range = 5:2:100, startvalue = 15),
# )
# on(sg_dims.sliders[1].value) do val; nx[] = val; end
# on(sg_dims.sliders[2].value) do val; ny[] = val; end

# # Section 2: Deformation Settings
# Label(control_panel[3, 1:2], "Deformation Settings", font = :bold, tellwidth=false, halign=:left, padding=(0,0,10,0))
# # Create a Toggle
# deform_toggle = Toggle(control_panel[4, 1], active = true)
# Label(control_panel[4, 2], "Enable Deformation", tellwidth=false, halign=:left)
# on(deform_toggle.active) do is_active; deformation_active[] = is_active; end
# # Sliders for deformation parameters
# sg_deform = SliderGrid(
#     control_panel[5, 1:2],
#     (label = "Amplitude", range = 0.0:0.01:0.2, startvalue = 0.05),
#     (label = "Frequency", range = 0.5:0.1:10.0, startvalue = 2.0),
# )
# on(sg_deform.sliders[1].value) do val; amplitude[] = val; end
# on(sg_deform.sliders[2].value) do val; frequency[] = val; end

# # Section 3: Algorithm Selection
# Label(control_panel[6, 1:2], "Algorithm", font = :bold, tellwidth=false, halign=:left, padding=(0,0,10,0))
# # Create a Menu
# grid_menu = Menu(control_panel[7, 1:2], options = ["Deformed", "Orthogonal", "Simple Rect"], default = "Deformed")
# on(grid_menu.selection) do selected_type; grid_type[] = selected_type; end

# Section 4: Actions
# Label(control_panel[8, 1:2], "Actions", font = :bold, tellwidth=false, halign=:left, padding=(0,0,10,0))
# button_layout = GridLayout(control_panel[9, 1:2])
# reset_button = Button(button_layout[1, 1], label = "Reset View")
# save_button = Button(button_layout[1, 2], label = "Save Grid")

# on(reset_button.clicks) do _; autolimits!(ax); autolimits!(ax_heatmap); end
# on(save_button.clicks) do _
#     X, Y = grid_coords[]; points_to_save = hcat(vec(X), vec(Y))
#     filename = "grid_output.csv"; writedlm(filename, points_to_save, ',')
#     println("Grid saved to $(filename)!")
# end


# --- MODIFIED Grid Generation Logic ---
function generate_grid(nx, ny, is_active, amp, freq, type)
    x_range = range(0, 1, length=nx)
    y_range = range(0, 1, length=ny)
    X = [x for x in x_range, y in y_range]
    Y = [y for x in x_range, y in y_range]
    amp = 1.0
    freq = 0.4

    if type == "Deformed" && is_active
        # Add deformation only if the type is correct and toggle is on
        Y .+= amp .* sin.(freq * π .* X)
    elseif type == "Orthogonal" # Example for a different grid type
         X .+= amp .* sin.(freq * π .* Y) # swap for a different effect
    end
    # If type is "Simple Rect", no modification is done.
    return X, Y
end

# The @lift now depends on the new observables from the toggle and menu
# grid_coords = @lift(generate_grid($nx, $ny, $deformation_active, $amplitude, $frequency, $grid_type))

# --- Helper Functions and Plotting (Unchanged) ---
# ... (to_linesegments, extract_cell_polygons, calculate_aspect_ratio) ...

# function to_linesegments(X,Y); nx,ny=size(X); points=Point2f[]; for j in 1:ny, i in 1:nx-1; push!(points,Point2f(X[i,j],Y[i,j]),Point2f(X[i+1,j],Y[i+1,j])); end; for i in 1:nx, j in 1:ny-1; push!(points,Point2f(X[i,j],Y[i,j]),Point2f(X[i,j+1],Y[i,j+1])); end; return points; end
# function extract_cell_polygons(X,Y); nx,ny=size(X); polygons=Polygon[]; for j in 1:(ny-1),i in 1:(nx-1); p1=Point2f(X[i,j],Y[i,j]);p2=Point2f(X[i+1,j],Y[i+1,j]);p3=Point2f(X[i+1,j+1],Y[i+1,j+1]);p4=Point2f(X[i,j+1],Y[i,j+1]); push!(polygons,Polygon([p1,p2,p3,p4])); end; return polygons; end
# function calculate_aspect_ratio(X,Y); nx,ny=size(X); quality_values=Float32[]; for j in 1:(ny-1),i in 1:(nx-1); p1=(X[i,j],Y[i,j]);p2=(X[i+1,j],Y[i+1,j]);p3=(X[i,j+1],Y[i,j+1]); e1_sq=(p2[1]-p1[1])^2+(p2[2]-p1[2])^2;e2_sq=(p3[1]-p1[1])^2+(p3[2]-p1[2])^2; ar=(e1_sq<1e-10||e2_sq<1e-10) ? 1000.0 : max(e1_sq,e2_sq)/min(e1_sq,e2_sq); push!(quality_values,ar); end; return quality_values; end

# segments = @lift(to_linesegments($grid_coords...))
# cell_polygons = @lift(extract_cell_polygons($grid_coords...))
# cell_colors = @lift(calculate_aspect_ratio($grid_coords...))

# linesegments!(ax, segments, color = :blue, linewidth=1.5)
# poly_plot = poly!(ax_heatmap, cell_polygons, color = cell_colors, colormap = :thermal, strokecolor = :black, strokewidth = 0.5)
# Colorbar(vis_layout[1:2, 2], poly_plot, label = "Aspect Ratio")

# linkaxes!(ax, ax_heatmap)
# rowsize!(vis_layout, 2, Relative(0.9)) # Give heatmap more vertical space
# autolimits!(ax); autolimits!(ax_heatmap)

fig