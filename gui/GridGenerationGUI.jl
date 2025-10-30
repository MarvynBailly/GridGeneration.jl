"""
GridGeneration.jl Interactive GUI

A GLMakie-based graphical user interface for interactive grid generation, 
supporting domain splitting, boundary solving, and elliptic smoothing.

Usage:
    julia> include("GridGenerationGUI.jl")
    # GUI window will appear
"""

# === Library Imports ===
include("C:\\Users\\admin\\Documents\\GitHub\\GridGeneration\\src\\GridGeneration.jl")
using .GridGeneration

using GLMakie
using Makie: Polygon
using DelimitedFiles
using Dates

# === Include GUI Modules ===
include("utils/helpers.jl")
include("components/ui_components.jl")
include("plotters/plot_blocks.jl")
include("handlers/event_handlers.jl")
include("setup/domain_setup.jl")

# ========================================
# === Domain Setup ===
# ========================================

# Create rectangular domain with default metric
width = 4.0
height = 2.0
num_points_width = 100
num_points_height = 50

initialGrid, initialBndInfo, initialInterfaceInfo = setup_rectangle_domain(
    width, height, num_points_width, num_points_height
)

M = create_default_metric(1000)

# Alternative: Load from Turtle grid files
# Uncomment the following lines to use Turtle grid data instead:

# metricFieldFile = "step/BFstepTest_entropy.metric"
# gridFolder = "step/coarseGrids"
# initialGrid, initialBndInfo, initialInterfaceInfo, M = setup_turtle_grid_domain(
#     metricFieldFile, gridFolder
# )




# ========================================
# === GUI Layout Setup ===
# ========================================

set_theme!(theme_light())
fig = Figure(size = (1800, 1000), backgroundcolor = RGBf(0.85, 0.85, 0.85))  # Light gray background

# Main layout: 3 columns [controls | plots | console]
left_layout = fig[1, 1] = GridLayout(tellwidth = false, halign = :left, valign = :top)
center_layout = fig[1, 2] = GridLayout()
right_layout = fig[1, 3] = GridLayout(tellwidth = false, halign = :left, valign = :top)

colsize!(fig.layout, 1, Relative(1/6))   # Controls - 16%
colsize!(fig.layout, 2, Relative(4/6))   # Plots - 67%
colsize!(fig.layout, 3, Relative(1/6))   # Console - 16%

# Sub-layouts
control_panel = left_layout[1, 1] = GridLayout()

plot_panel = center_layout[1, 1] = GridLayout()
tool_bar = center_layout[2, 1] = GridLayout()
rowsize!(center_layout, 1, Relative(9/10))
rowsize!(center_layout, 2, Relative(1/10))

console_panel = right_layout[1, 1] = GridLayout()

# Decorative boxes
Box(fig[1, 1], linestyle = :solid, strokecolor = :black, cornerradius = 20)
Box(fig[1, 3], linestyle = :solid, strokecolor = :black, cornerradius = 20)
Box(center_layout[1, 1], linestyle = :solid, strokecolor = :black, cornerradius = 20)
Box(center_layout[2, 1], linestyle = :solid, strokecolor = :black, cornerradius = 20)

# ========================================
# === Visualization Panels ===
# ========================================

ax_initial = Axis(plot_panel[1, 1], aspect = DataAspect(), title = "Initial Grid (TFI)")
ax_generated = Axis(plot_panel[2, 1], aspect = DataAspect(), title = "Split/Edge Solved Grid")
ax_final = Axis(plot_panel[3, 1], aspect = DataAspect(), title = "Final Grid")

# Create overlay panel for highlight toggle in top right corner
overlay_panel = GridLayout(plot_panel[1, 1], 
                          tellwidth = false, 
                          tellheight = false,
                          halign = :right,
                          valign = :top)

# Add semi-transparent background box for the overlay
overlay_box = Box(overlay_panel[1, 1:2], 
                  color = RGBAf(1, 1, 1, 0.8),
                  strokecolor = :gray50,
                  cornerradius = 10,
                  )

# Create toggle and label in the overlay
highlight_toggle = Toggle(overlay_panel[1, 1], active = false)
Label(overlay_panel[1, 2], "Highlight", 
      halign = :left, 
      padding = (5, 10, 5, 10))

# ========================================
# === Observables ===
# ========================================

# Console observable
console_text = Observable("")

# Elliptic parameters observable (one per block)
elliptic_params_obs = Observable([GridGeneration.EllipticParams() for _ in 1:length(initialGrid)])

# Grid observables
initial_blocks_obs = Observable(deepcopy(initialGrid))
initial_bndInfo_obs = Observable(deepcopy(initialBndInfo))
initial_interfaceInfo_obs = Observable(deepcopy(initialInterfaceInfo))

split_blocks = Observable(deepcopy(initialGrid))
split_bndInfo = Observable(deepcopy(initialBndInfo))
split_interfaceInfo = Observable(deepcopy(initialInterfaceInfo))

generated_blocks = Observable(deepcopy(initialGrid))
generated_bndInfo = Observable(deepcopy(initialBndInfo))
generated_interfaceInfo = Observable(deepcopy(initialInterfaceInfo))

final_blocks = Observable(deepcopy(initialGrid))
final_bndInfo = Observable(deepcopy(initialBndInfo))
final_interfaceInfo = Observable(deepcopy(initialInterfaceInfo))

# Split line preview observables
i_split_lines = Observable(Point2f[])
j_split_lines = Observable(Point2f[])

# ========================================
# === UI Component Setup ===
# ========================================

controls = populate_control_panel!(control_panel, initialGrid)
tool = populate_button_panel!(tool_bar)

# Create console panel
console = create_console_panel!(console_panel, console_text)

# Log initial message
log_to_console(console_text, "Grid Generation GUI initialized.")

# ========================================
# === Plotting ===
# ========================================

# Initial grid plot with static boundary info
plot_blocks_with_highlighting!(
    ax_initial, 
    initial_blocks_obs, 
    highlight_toggle.active, 
    initial_bndInfo_obs, 
    initial_interfaceInfo_obs
)

# Generated grid plot with reactive boundary info
plot_blocks_with_highlighting!(
    ax_generated, 
    generated_blocks, 
    highlight_toggle.active,
    generated_bndInfo, 
    generated_interfaceInfo
)

# final grid plot with reactive boundary info
plot_blocks_with_highlighting!(
    ax_final, 
    final_blocks, 
    highlight_toggle.active,
    final_bndInfo, 
    final_interfaceInfo
)

# Add split line previews to initial plot
linesegments!(ax_initial, i_split_lines, color = :red, linewidth = 1.5, overdraw=true)
linesegments!(ax_initial, j_split_lines, color = :green, linewidth = 1.5, overdraw=true)

# ========================================
# === Event Handler Setup ===
# ========================================

# Split domain handler
setup_split_domain_handler!(
    tool[:split_domain], controls,
    initialGrid, initialBndInfo, initialInterfaceInfo,
    split_blocks, split_bndInfo, split_interfaceInfo,
    generated_blocks, generated_bndInfo, generated_interfaceInfo,
    final_blocks, final_bndInfo, final_interfaceInfo,
    console_text
)

# Edge solve handler
setup_edge_solve_handler!(
    tool[:edge_solve], controls, M,
    split_blocks, split_bndInfo, split_interfaceInfo,
    generated_blocks, generated_bndInfo, generated_interfaceInfo,
    final_blocks, final_bndInfo, final_interfaceInfo,
    console_text
)

# Smooth grid handler
setup_smooth_grid_handler!(
    tool[:smooth_grid],
    controls,
    generated_blocks,
    final_blocks,
    elliptic_params_obs,
    console_text
)

# Reset view handler
setup_reset_view_handler!(
    tool[:reset_view],
    ax_initial,
    ax_generated,
    ax_final
)

# Clear console handler
setup_clear_console_handler!(tool[:clear_console], console_text)

# Reset all handler
setup_reset_all_handler!(
    tool[:reset_all],
    controls,
    initialGrid, initialBndInfo, initialInterfaceInfo,
    split_blocks, split_bndInfo, split_interfaceInfo,
    generated_blocks, generated_bndInfo, generated_interfaceInfo,
    final_blocks, final_bndInfo, final_interfaceInfo,
    elliptic_params_obs, console_text
)

# Split line preview handlers
setup_split_line_preview!(controls[:i_splits], i_split_lines, initialGrid, :i)
setup_split_line_preview!(controls[:j_splits], j_split_lines, initialGrid, :j)

# Block parameter synchronization handler
setup_block_parameter_sync!(controls, elliptic_params_obs, generated_blocks,console_text)

# Split block selector handler (for initial grid only)
setup_split_block_selector!(controls, initialGrid)

# ========================================
# === Display Figure ===
# ========================================

fig
