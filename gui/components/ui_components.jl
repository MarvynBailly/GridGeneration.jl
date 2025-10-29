"""
UI component creation functions for the Grid Generation GUI.
"""

"""
    populate_control_panel!(parent_layout::GridLayout)

Create and populate the control panel with input widgets for grid generation parameters.

Returns a dictionary of widgets with the following keys:
- `:i_splits`, `:j_splits` - Textboxes for split locations
- `:edge_solver` - Menu for boundary solver type
- `:smoothing_type` - Menu for smoothing method
- `:max_iter` - Textbox for maximum iterations
- `:tolerance` - Textbox for convergence tolerance
- `:omega` - Textbox for SOR relaxation parameter
- `:forcing_left_a`, `:forcing_left_b` - Wall forcing parameters for left boundary
- `:forcing_right_a`, `:forcing_right_b` - Wall forcing parameters for right boundary
- `:forcing_bottom_a`, `:forcing_bottom_b` - Wall forcing parameters for bottom boundary
- `:forcing_top_a`, `:forcing_top_b` - Wall forcing parameters for top boundary
"""
function populate_control_panel!(parent_layout::GridLayout)
    widgets = Dict{Symbol, Any}()

    # Title
    Label(parent_layout[1, 1:2], "Parameters", font = :bold, halign = :center, 
          padding=(10,0,0,10), valign=:top)
    
    # === Split Locations Section ===
    Label(parent_layout[2, 1:2], "Splits Locations", font = :bold, halign = :left, 
          padding=(10,0,10,20), valign=:top)
    
    # Block selector for split locations
    Label(parent_layout[3, 1], "Block to Split:", halign=:left, padding=(10,0,0,0))
    widgets[:split_block_selector] = Menu(parent_layout[3, 2], options = ["Block 1"], 
                                           default = "Block 1")
    
    Label(parent_layout[4, 1], "I-splits (e.g., 10,25)", halign=:left, valign=:top, 
          padding=(10,0,0,0))
    Label(parent_layout[4, 2], "J-splits (e.g., 8)")
    widgets[:i_splits] = Textbox(parent_layout[5, 1], placeholder="none", defocus_on_submit=true)
    widgets[:j_splits] = Textbox(parent_layout[5, 2], placeholder="none", defocus_on_submit=true)

    # === Boundary Solver Section ===
    Label(parent_layout[6, 1:2], "Boundary Solver Parameters", font = :bold, halign = :left, 
          padding=(10,0,10,20), valign=:top)
    Label(parent_layout[7, 1], "Edge Solver Type:", halign=:left, padding=(10,0,0,0))
    widgets[:edge_solver] = Menu(parent_layout[7, 2], options = ["analytic", "numerical"], 
                                  default = "analytic")

    # === Smoothing Parameters Section ===
    Label(parent_layout[8, 1:2], "Smoothing Parameters", font = :bold, halign = :left, 
          padding=(10,0,10,20), valign=:top)
    
    # Block selector for multi-block smoothing parameters
    Label(parent_layout[9, 1], "Block:", halign=:left, padding=(10,0,0,0))
    widgets[:block_selector] = Menu(parent_layout[9, 2], options = ["Block 1"], 
                                     default = "Block 1")
    
    Label(parent_layout[10, 1], "Smoothing Type:", halign=:left, padding =(10,0,0,0))
    widgets[:smoothing_type] = Menu(parent_layout[10, 2], options = ["Elliptic-SS"], 
                                     default = "Elliptic-SS")

    Label(parent_layout[11, 1], "Max Iterations:", padding=(10,0,0,0))
    widgets[:max_iter] = Textbox(parent_layout[11, 2], placeholder="5000", validator=Int, stored_string="5000")
    
    Label(parent_layout[12, 1], "Tolerance:", padding=(10,0,0,0))
    widgets[:tolerance] = Textbox(parent_layout[12, 2], placeholder="1e-5", validator=Float64, stored_string="1e-5")

    Label(parent_layout[13, 1], "Omega:", halign=:left, padding=(10,0,0,0))
    widgets[:omega] = Textbox(parent_layout[13, 2], placeholder="0.2", validator=Float64, stored_string="0.2")

    # === Wall Forcing Parameters Section ===
    Label(parent_layout[14, 1:2], "Wall Forcing Parameters", font = :bold, halign = :left, 
          padding=(10,0,10,20))
    
    forcing_grid = parent_layout[14, 1:2] = GridLayout()
    Label(forcing_grid[1, 2], "a", halign=:center)
    Label(forcing_grid[1, 3], "b", halign=:center)
    
    Label(forcing_grid[2, 1], "Left:", halign=:right)
    widgets[:forcing_left_a] = Textbox(forcing_grid[2, 2], placeholder="0.4", validator=Float64, stored_string="0.4")
    widgets[:forcing_left_b] = Textbox(forcing_grid[2, 3], placeholder="0.4", validator=Float64, stored_string="0.4")
    
    Label(forcing_grid[3, 1], "Right:", halign=:right)
    widgets[:forcing_right_a] = Textbox(forcing_grid[3, 2], placeholder="0.4", validator=Float64, stored_string="0.4")
    widgets[:forcing_right_b] = Textbox(forcing_grid[3, 3], placeholder="0.4", validator=Float64, stored_string="0.4")
    
    Label(forcing_grid[4, 1], "Bottom:", halign=:right)
    widgets[:forcing_bottom_a] = Textbox(forcing_grid[4, 2], placeholder="0.4", validator=Float64, stored_string="0.4")
    widgets[:forcing_bottom_b] = Textbox(forcing_grid[4, 3], placeholder="0.4", validator=Float64, stored_string="0.4")
    
    Label(forcing_grid[5, 1], "Top:", halign=:right)
    widgets[:forcing_top_a] = Textbox(forcing_grid[5, 2], placeholder="0.4", validator=Float64, stored_string="0.4")
    widgets[:forcing_top_b] = Textbox(forcing_grid[5, 3], placeholder="0.4", validator=Float64, stored_string="0.4")
    
    Label(forcing_grid[6, 1], " ", halign=:right)

    return widgets
end

"""
    populate_button_panel!(parent_layout)

Create and populate the button panel with action buttons and toggles.

Returns a dictionary of buttons and toggles with the following keys:
- `:split_domain` - Button to trigger domain splitting
- `:edge_solve` - Button to trigger edge solving
- `:smooth_grid` - Button to trigger grid smoothing
- `:save_grid` - Button to save the current grid
- `:reset_view` - Button to reset axis limits
- `:highlight_boundaries` - Toggle to highlight domain boundaries
- `:clear_console` - Button to clear the console output
"""
function populate_button_panel!(parent_layout)
    btns = Dict{Symbol, Any}()

    btns[:split_domain] = Button(parent_layout[1, 1], label = "Split Domain")
    btns[:edge_solve] = Button(parent_layout[1, 2], label = "Solve Edges")
    btns[:smooth_grid] = Button(parent_layout[1, 3], label = "Smooth Grid")
    btns[:save_grid] = Button(parent_layout[1, 4], label = "Save Grid")
    btns[:reset_view] = Button(parent_layout[1, 5], label = "Reset View")

    btns[:highlight_boundaries] = Toggle(parent_layout[1, 6], active = false)
    Label(parent_layout[1, 7], "Highlight Boundaries", halign=:left, padding=(0,0,0,0))
    
    btns[:clear_console] = Button(parent_layout[1, 8], label = "Clear Console")

    return btns
end

"""
    create_console_panel!(parent_layout, console_text_obs)

Create a scrollable console panel for displaying log messages.

# Arguments
- `parent_layout`: The GridLayout to add the console to
- `console_text_obs`: Observable{String} containing the console text

# Returns
- Dictionary with `:textbox` key containing the Label widget
"""
function create_console_panel!(parent_layout, console_text_obs)
    console_dict = Dict{Symbol, Any}()
    
    # Title
    Label(parent_layout[1, 1], "Console Output", font = :bold, halign = :left, 
          padding=(5,10,5,10), valign=:top)
    
    # Scrollable area for console text
    console_scroll = parent_layout[2, 1:2] = GridLayout()
    
    # Text display area - using Axis with text! for scrollable multiline text
    ax = Axis(console_scroll[1, 1], 
              xlabel="", ylabel="",
              xticklabelsvisible=false, yticklabelsvisible=false,
              xticksvisible=false, yticksvisible=false,
              leftspinevisible=false, rightspinevisible=false,
              topspinevisible=false, bottomspinevisible=false,
              xgridvisible=false, ygridvisible=false)
    
    # Add text element that updates with console_text_obs
    txt = text!(ax, 0, 1, text=console_text_obs, 
                fontsize=10, 
                align=(:left, :top),
                space=:relative)
    
    # Set axis limits
    limits!(ax, -0.02, 1, 0, 1.02)
    
    console_dict[:axis] = ax
    console_dict[:text] = txt
    
    return console_dict
end
