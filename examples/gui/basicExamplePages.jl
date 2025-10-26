include("../../src/GridGeneration.jl")
using .GridGeneration

using GLMakie
using Makie: Polygon
using DelimitedFiles

set_theme!(theme_light())

fig = Figure(size = (1400, 900))

# --- Main Layout ---
# Plotting area on the left, control panel on the right
vis_layout = fig[1, 2] = GridLayout()
control_panel = fig[1, 1] = GridLayout(tellwidth = false, width = 350)



# --- Visualization Panels (Axes) ---
ax_initial = Axis(vis_panel[1, 1], aspect = DataAspect(), title = "Initial Grid (TFI)")
ax_generated = Axis(vis_panel[2, 1], aspect = DataAspect(), title = "Generated Grid")
ax_smoothed = Axis(vis_panel[3, 1], aspect = DataAspect(), title = "Smoothed Grid")


# --- Observables for Grid Data ---
# These will hold the (X, Y) tuples for each grid state
initial_grid = Observable{Any}(nothing)
generated_grid = Observable{Any}(nothing)
smoothed_grid = Observable{Any}(nothing)

# --- Control Panel Setup ---
Label(control_panel[1, 1:2], "Control Section", font=:bold)
page_menu = Menu(control_panel[2, 1:2], options=["Splits", "Solver", "Smoothing"])
page_container = control_panel[3, 1:2] = GridLayout()


# Create a separate GridLayout for each page of controls
splits_page = page_container[1,1] = GridLayout(tellwidth=false)
solver_page = page_container[1,1] = GridLayout(tellwidth=false, visible=false)
smoothing_page = page_container[1,1] = GridLayout(tellwidth=false, visible=false)

pages = [splits_page, solver_page, smoothing_page]



# Splits Configuration
Label(splits_page[1, 1:2], "Splits Locations", font = :bold, halign = :left, padding=(0,0,10,0))
Label(splits_page[2, 1], "I-splits (e.g., 10,25)")
Label(splits_page[2, 2], "J-splits (e.g., 8)")
i_splits_box = Textbox(splits_page[3, 1], placeholder="none")
j_splits_box = Textbox(splits_page[3, 2], placeholder="none")

split_domain_button = Button(splits_page[4, 1], label = "Split Domain")


# Solver & Smoothing Parameters
Label(solver_page[1, 1:2], "Boundary Solver Parameters", font = :bold, halign = :left, padding=(0,0,10,0))

# Solver Settings
Label(solver_page[2, 1], "Enable Edge Solver", halign=:left)
use_edge_solver_toggle = Toggle(solver_page[2, 2], active=true)
Label(solver_page[3, 1], "Edge Solver Type", halign=:left)
edge_solver = Menu(solver_page[3, 2], options = ["analytic", "numerical"], default = "analytic")

# Smoothing Settings
Label(smoothing_page[1, 1:2], "Smoothing Parameters", font = :bold, halign = :left, padding=(0,0,10,0))

# Solver Settings
Label(smoothing_page[2, 1], "Enable Smoothing", halign=:left)
use_smoothing_toggle = Toggle(smoothing_page[2, 2], active=true)
Label(smoothing_page[3, 1], "Smoothing Type", halign=:left)
smoothing_type = Menu(smoothing_page[3, 2], options = ["Elliptic-SS"], default = "Elliptic-SS")

Label(smoothing_page[4, 1], "Max Iterations:")
max_iter_box = Textbox(smoothing_page[4, 2], placeholder="5000", validator=Int)
Label(smoothing_page[5, 1], "Tolerance:")
tolerance_box = Textbox(smoothing_page[5, 2], placeholder="1e-5", validator=Float64)

Label(smoothing_page[6, 1], ("Omega: "))
omega_box = Textbox(smoothing_page[6, 2], placeholder="0.2", validator=Float64)

Label(smoothing_page[7, 1:2], "Wall Forcing Parameters", font = :bold, halign = :left, padding=(0,0,10,0))
forcing_grid = smoothing_page[8, 1:2] = GridLayout()
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

# Function to set visibility of all children in a GridLayout
function set_page_visibility(page::GridLayout, is_visible::Bool)
    for child in page.content
        child.content.visible = is_visible
    end
end

# Initial setup: hide all but the first page
for (i, page) in enumerate(pages)
    set_page_visibility(page, i == 1)
end

# Page switching logic
on(page_menu.selection) do selected_page_name
    selected_index = findfirst(==(selected_page_name), page_menu.options[])
    if !isnothing(selected_index)
        for (i, page) in enumerate(pages)
            set_page_visibility(page, i == selected_index)
        end
    end
endd


# # Section 1: Grid Dimensions
# Label(smoothing_page[1, 1:2], "Grid Dimensions", font = :bold, tellwidth=false, halign=:left)
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
grid_coords = @lift(generate_grid($nx, $ny, $deformation_active, $amplitude, $frequency, $grid_type))

# --- Helper Functions and Plotting (Unchanged) ---
# ... (to_linesegments, extract_cell_polygons, calculate_aspect_ratio) ...

function to_linesegments(X,Y); nx,ny=size(X); points=Point2f[]; for j in 1:ny, i in 1:nx-1; push!(points,Point2f(X[i,j],Y[i,j]),Point2f(X[i+1,j],Y[i+1,j])); end; for i in 1:nx, j in 1:ny-1; push!(points,Point2f(X[i,j],Y[i,j]),Point2f(X[i,j+1],Y[i,j+1])); end; return points; end
function extract_cell_polygons(X,Y); nx,ny=size(X); polygons=Polygon[]; for j in 1:(ny-1),i in 1:(nx-1); p1=Point2f(X[i,j],Y[i,j]);p2=Point2f(X[i+1,j],Y[i+1,j]);p3=Point2f(X[i+1,j+1],Y[i+1,j+1]);p4=Point2f(X[i,j+1],Y[i,j+1]); push!(polygons,Polygon([p1,p2,p3,p4])); end; return polygons; end
function calculate_aspect_ratio(X,Y); nx,ny=size(X); quality_values=Float32[]; for j in 1:(ny-1),i in 1:(nx-1); p1=(X[i,j],Y[i,j]);p2=(X[i+1,j],Y[i+1,j]);p3=(X[i,j+1],Y[i,j+1]); e1_sq=(p2[1]-p1[1])^2+(p2[2]-p1[2])^2;e2_sq=(p3[1]-p1[1])^2+(p3[2]-p1[2])^2; ar=(e1_sq<1e-10||e2_sq<1e-10) ? 1000.0 : max(e1_sq,e2_sq)/min(e1_sq,e2_sq); push!(quality_values,ar); end; return quality_values; end

segments = @lift(to_linesegments($grid_coords...))
cell_polygons = @lift(extract_cell_polygons($grid_coords...))
cell_colors = @lift(calculate_aspect_ratio($grid_coords...))

linesegments!(ax, segments, color = :blue, linewidth=1.5)
poly_plot = poly!(ax_heatmap, cell_polygons, color = cell_colors, colormap = :thermal, strokecolor = :black, strokewidth = 0.5)
Colorbar(vis_layout[1:2, 2], poly_plot, label = "Aspect Ratio")

linkaxes!(ax, ax_heatmap)
rowsize!(vis_layout, 2, Relative(0.9)) # Give heatmap more vertical space
autolimits!(ax); autolimits!(ax_heatmap)

fig