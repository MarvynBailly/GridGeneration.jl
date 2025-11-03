# Grid Generation GUI - Code Organization

This directory contains the interactive GUI for GridGeneration.jl, organized into logical modules for better maintainability and readability.

## File Structure

```
gui/
├── GridGenerationGUI.jl         # Main entry point - launches the GUI
├── components/
│   └── ui_components.jl         # UI widget creation functions
├── plotters/
│   └── plot_blocks.jl           # Grid visualization functions
├── handlers/
│   └── event_handlers.jl        # Button click and event callbacks
├── setup/
│   └── domain_setup.jl          # Domain initialization utilities
└── utils/
    └── helpers.jl               # General helper functions
```

## Module Descriptions

### `GridGenerationGUI.jl` (Main File)
The main entry point that:
- Sets up the GLMakie figure and layout
- Creates observables for reactive updates
- Assembles all components
- Displays the GUI window

**Usage:**
```julia
include("GridGenerationGUI.jl")
```

### `components/ui_components.jl`
UI widget creation and layout functions:
- `populate_control_panel!()` - Creates parameter input controls
- `populate_button_panel!()` - Creates action buttons and toggles

### `plotters/plot_blocks.jl`
Visualization functions for grid rendering:
- `plot_blocks_with_highlighting!()` - Reactive grid plotting with boundary/interface highlighting
- `extract_face_xy()` - Helper to extract boundary faces from blocks
- `ordered_range()` - Helper for index range calculation

### `handlers/event_handlers.jl`
Event callback functions:
- `setup_split_domain_handler!()` - Handles domain splitting
- `setup_edge_solve_handler!()` - Handles edge solving
- `setup_reset_view_handler!()` - Handles view reset
- `setup_split_line_preview!()` - Handles split line visualization

### `setup/domain_setup.jl`
Domain initialization utilities:
- `setup_rectangle_domain()` - Creates a simple rectangular domain
- `setup_turtle_grid_domain()` - Loads domain from Turtle grid files
- `create_default_metric()` - Creates a default metric function

### `utils/helpers.jl`
General utility functions:
- `parse_indices()` - Parses comma-separated indices from strings
- `to_linesegments()` - Converts 2D arrays to line segments

## Customization Guide

### Adding a New Button/Control

1. **Add widget in** `components/ui_components.jl`:
   ```julia
   widgets[:my_button] = Button(parent_layout[row, col], label = "My Action")
   ```

2. **Create handler in** `handlers/event_handlers.jl`:
   ```julia
   function setup_my_action_handler!(button, ...)
       on(button.clicks) do _
           # Your action code here
       end
   end
   ```

3. **Wire it up in** `GridGenerationGUI.jl`:
   ```julia
   setup_my_action_handler!(tool[:my_button], ...)
   ```

### Changing the Default Domain

Edit the domain setup section in `GridGenerationGUI.jl`:

```julia
# For custom dimensions
initialGrid, initialBndInfo, initialInterfaceInfo = setup_rectangle_domain(
    8.0,  # width
    4.0,  # height
    200,  # num_points_width
    100   # num_points_height
)

# For custom metric
M = (x, y) -> [1000 * sin(x), 1000 * cos(y)]
```

### Adding New Visualization Panels

1. Add axis in the layout section:
   ```julia
   ax_new = Axis(plot_panel[4, 1], aspect = DataAspect(), title = "My View")
   ```

2. Create observable:
   ```julia
   new_data = Observable(initial_data)
   ```

3. Add plotting call:
   ```julia
   plot_blocks_with_highlighting!(ax_new, new_data, highlight_toggle.active, ...)
   ```

## Development Notes

- All functions that modify the GUI layout use the `!` convention
- Observables are used for reactive updates - modify with `obs[] = new_value`
- Event handlers are set up using `on(observable) do ... end` pattern
- The GUI uses GLMakie's `GridLayout` for responsive positioning

## Backup

A backup of the original monolithic file is saved as `GridGenerationGUI_backup.jl`.
