"""
Event handler functions for the Grid Generation GUI.
"""

"""
    setup_split_domain_handler!(button, controls, initialGrid, initialBndInfo, initialInterfaceInfo, 
                                generated_blocks, generated_bndInfo, generated_interfaceInfo, console_obs)

Set up the event handler for the Split Domain button.
Parses split locations from textboxes and calls GridGeneration.SplitBlock.
"""
function setup_split_domain_handler!(button, controls, initialGrid, initialBndInfo, initialInterfaceInfo,
                                      split_blocks, split_bndInfo, split_interfaceInfo, 
                                      generated_blocks, generated_bndInfo, generated_interfaceInfo,
                                      final_blocks, final_bndInfo, final_interfaceInfo,
                                      console_obs)
    on(button.clicks) do _
        # Get which block to split
        block_str = string(controls[:split_block_selector].selection[])
        block_num = parse(Int, split(block_str, " ")[end])
        
        i_splits_str = controls[:i_splits].stored_string[]
        j_splits_str = controls[:j_splits].stored_string[]
        i_indices = parse_indices(i_splits_str)
        j_indices = parse_indices(j_splits_str)

        if isempty(i_indices) && isempty(j_indices)
            log_to_console(console_obs, "No splits provided.")
            return
        end

        split_locations = [i_indices, j_indices]
        split_dict = Dict(block_num => split_locations)

        log_to_console(console_obs, "Splitting Block $block_num at i-indices: $(i_indices), j-indices: $(j_indices)...")

        split_block, split_bnd_info, split_inter_info = GridGeneration.SplitMultiBlock(
            deepcopy(initialGrid), 
            split_dict,
            deepcopy(initialBndInfo), 
            deepcopy(initialInterfaceInfo)
        )
        
        # Get the current blocks (or initial grid if first split)
        split_blocks[] = split_block
        split_bndInfo[] = split_bnd_info
        split_interfaceInfo[] = split_inter_info
        generated_blocks[] = split_block
        generated_bndInfo[] = split_bnd_info
        generated_interfaceInfo[] = split_inter_info
        final_blocks[] = split_block
        final_bndInfo[] = split_bnd_info
        final_interfaceInfo[] = split_inter_info

        # clear split inputs
        controls[:i_splits].stored_string[] = ""
        controls[:j_splits].stored_string[] = ""
        controls[:i_splits].displayed_string[] = "none"
        controls[:j_splits].displayed_string[] = "none"



        log_to_console(console_obs, "Domain splitting complete. Generated $(length(generated_blocks[])) block(s).")
    end
end

"""
    setup_edge_solve_handler!(button, controls, M, generated_blocks, generated_bndInfo, generated_interfaceInfo, console_obs)

Set up the event handler for the Solve Edges button.
Solves the grid edges using the selected solver method.
"""
function setup_edge_solve_handler!(button, controls, M, 
    split_blocks, split_bndInfo, split_interfaceInfo,
    generated_blocks, generated_bndInfo, generated_interfaceInfo, 
    final_blocks, final_bndInfo, final_interfaceInfo,
    console_obs)
    on(button.clicks) do _
        edge_type = controls[:edge_solver].selection[]
        solver_sym = edge_type == "analytic" ? :analytic : :numerical
        
        current_blocks = split_blocks[]
        current_bnd_info = split_bndInfo[]
        current_inter_info = split_interfaceInfo[]

        log_to_console(console_obs, "Solving edges for $(length(current_blocks)) block(s) using $solver_sym solver...")
        
        SolvedBlocks, SolvedBndInfo, solvedInterInfo = GridGeneration.SolveAllBlocks(
            M, 
            deepcopy(current_blocks), 
            deepcopy(current_bnd_info), 
            deepcopy(current_inter_info); 
            solver=solver_sym
        )

        generated_blocks[] = SolvedBlocks
        generated_bndInfo[] = SolvedBndInfo
        generated_interfaceInfo[] = solvedInterInfo
        final_blocks[] = SolvedBlocks
        final_bndInfo[] = SolvedBndInfo
        final_interfaceInfo[] = solvedInterInfo

        log_to_console(console_obs, "Edge solving complete.")
    end
end

"""
    smooth_grid_handler!(button, controls, M, generated_blocks, generated_bndInfo, generated_interfaceInfo, console_obs)

Set up the event handler for the Smooth Grid button.
"""
function setup_smooth_grid_handler!(button, controls, generated_blocks, final_blocks, elliptic_params_obs, console_obs)
    on(button.clicks) do _
        smooth_type = controls[:smoothing_type].selection[]
        if smooth_type == "Elliptic-SS"
            solver_sym = :ellipticSS
        end

        current_blocks = generated_blocks[]
        
        # Use the elliptic parameters from the observable (one per block)
        elliptic_params = elliptic_params_obs[]

        log_to_console(console_obs, "Smoothing grid for $(length(current_blocks)) block(s) using $solver_sym smoothing...")

        smoothBlocks, finalErrors, finalIterations = GridGeneration.SmoothBlocks(current_blocks; solver=solver_sym, params=elliptic_params)
        
        final_blocks[] = smoothBlocks

        for i in eachindex(finalErrors)
            log_to_console(console_obs, "Block $i:")
            log_to_console(console_obs, "  Final Error: $(finalErrors[i])")
            log_to_console(console_obs, "  Iterations: $(finalIterations[i])")
        end
        log_to_console(console_obs, "Smoothing complete.")
    end
end




"""
    setup_reset_view_handler!(button, axes...)

Set up the event handler for the Reset View button.
Resets the axis limits for all provided axes.
"""
function setup_reset_view_handler!(button, axes...)
    on(button.clicks) do _
        for ax in axes
            autolimits!(ax)
        end
    end
end

"""
    setup_clear_console_handler!(button, console_obs)

Set up the event handler for the Clear Console button.
Clears all text from the console.
"""
function setup_clear_console_handler!(button, console_obs)
    on(button.clicks) do _
        clear_console(console_obs)
    end
end

"""
    setup_reset_all_handler!(button, controls, initialGrid, initialBndInfo, initialInterfaceInfo,
                             generated_blocks, generated_bndInfo, generated_interfaceInfo,
                             smoothed_blocks, elliptic_params_obs, console_obs)

Set up the event handler for the Reset All button.
Resets all grids to initial state and all parameters to default values.
"""
function setup_reset_all_handler!(button, controls, 
                                    initialGrid, initialBndInfo, initialInterfaceInfo,
                                    split_blocks, split_bndInfo, split_interfaceInfo,
                                   generated_blocks, generated_bndInfo, generated_interfaceInfo,
                                   final_blocks, final_bndInfo, final_interfaceInfo,
                                   elliptic_params_obs, console_obs)
    on(button.clicks) do _
        # Reset all grid observables to initial state
        split_blocks[] = deepcopy(initialGrid)
        split_bndInfo[] = deepcopy(initialBndInfo)
        split_interfaceInfo[] = deepcopy(initialInterfaceInfo)
        generated_blocks[] = deepcopy(initialGrid)
        generated_bndInfo[] = deepcopy(initialBndInfo)
        generated_interfaceInfo[] = deepcopy(initialInterfaceInfo)
        final_blocks[] = deepcopy(initialGrid)
        final_bndInfo[] = deepcopy(initialBndInfo)
        final_interfaceInfo[] = deepcopy(initialInterfaceInfo)



        # Reset elliptic parameters to defaults
        elliptic_params_obs[] = [GridGeneration.EllipticParams() for _ in 1:length(initialGrid)]
        
        # Reset control panel parameters to defaults
        controls[:max_iter].displayed_string[] = "5000"
        controls[:tolerance].displayed_string[] = "1e-5"
        controls[:omega].displayed_string[] = "0.2"
        
        controls[:forcing_left_a].displayed_string[] = "0.4"
        controls[:forcing_left_b].displayed_string[] = "0.4"
        controls[:forcing_right_a].displayed_string[] = "0.4"
        controls[:forcing_right_b].displayed_string[] = "0.4"
        controls[:forcing_bottom_a].displayed_string[] = "0.4"
        controls[:forcing_bottom_b].displayed_string[] = "0.4"
        controls[:forcing_top_a].displayed_string[] = "0.4"
        controls[:forcing_top_b].displayed_string[] = "0.4"
        
        # Reset split textboxes
        controls[:i_splits].displayed_string[] = "none"
        controls[:i_splits].stored_string[] = "none"
        controls[:j_splits].displayed_string[] = "none"
        controls[:j_splits].stored_string[] = "none"
        
        # Reset block selectors to Block 1
        controls[:block_selector].selection[] = "Block 1"
        controls[:split_block_selector].selection[] = "Block 1"
        
        log_to_console(console_obs, "Reset all grids and parameters to initial state.")
    end
end

"""
    setup_split_line_preview!(textbox, observable, initialGrid, direction::Symbol)

Set up reactive preview lines for split locations.

# Arguments
- `textbox`: The textbox containing split indices
- `observable`: Observable to update with line segment points
- `initialGrid`: The initial grid to draw lines on
- `direction`: Either `:i` for vertical lines or `:j` for horizontal lines
"""
function setup_split_line_preview!(textbox, observable, initialGrid, direction::Symbol)
    on(textbox.stored_string) do s
        indices = parse_indices(s)
        points = Point2f[]
        
        if !isempty(indices)
            X, Y = initialGrid[1][1, :, :], initialGrid[1][2, :, :]
            nx, ny = size(X)
            
            if direction == :i
                # Vertical lines (constant i)
                for i in indices
                    if 1 <= i <= nx
                        for j in 1:(ny-1)
                            push!(points, Point2f(X[i, j], Y[i, j]))
                            push!(points, Point2f(X[i, j+1], Y[i, j+1]))
                        end
                    end
                end
            elseif direction == :j
                # Horizontal lines (constant j)
                for j in indices
                    if 1 <= j <= ny
                        for i in 1:(nx-1)
                            push!(points, Point2f(X[i, j], Y[i, j]))
                            push!(points, Point2f(X[i+1, j], Y[i, j]))
                        end
                    end
                end
            end
        end
        
        observable[] = points
    end
end

"""
    setup_block_parameter_sync!(controls, elliptic_params_obs, generated_blocks_obs, console_obs)

Set up synchronization between GUI controls and per-block elliptic parameters.
When block selector changes, load that block's parameters into the GUI.
When parameters change, update the selected block's parameters in the list.
"""
function setup_block_parameter_sync!(controls, elliptic_params_obs, generated_blocks_obs, console_obs)
    # Track which block is currently selected
    current_block_idx = Ref(1)
    
    # When block selector changes, load parameters for that block
    on(controls[:block_selector].selection) do selection
        block_str = string(selection)
        block_num = parse(Int, split(block_str, " ")[end])
        current_block_idx[] = block_num
        
        # Load parameters from elliptic_params for this block
        if block_num <= length(elliptic_params_obs[])
            params = elliptic_params_obs[][block_num]
            
            # Update GUI controls with this block's parameters
            # Use displayed_string to update the visible text
            controls[:max_iter].displayed_string[] = string(params.max_iter)
            controls[:tolerance].displayed_string[] = string(params.tol)
            controls[:omega].displayed_string[] = string(params.ω)
            
            controls[:forcing_left_a].displayed_string[] = string(params.a_decay_left)
            controls[:forcing_left_b].displayed_string[] = string(params.b_decay_left)
            controls[:forcing_right_a].displayed_string[] = string(params.a_decay_right)
            controls[:forcing_right_b].displayed_string[] = string(params.b_decay_right)
            controls[:forcing_bottom_a].displayed_string[] = string(params.a_decay_bottom)
            controls[:forcing_bottom_b].displayed_string[] = string(params.b_decay_bottom)
            controls[:forcing_top_a].displayed_string[] = string(params.a_decay_top)
            controls[:forcing_top_b].displayed_string[] = string(params.b_decay_top)
            
            log_to_console(console_obs, "Loaded parameters for Block $block_num")
        end
    end
    
    # When parameters change, update the elliptic_params list for current block
    function update_params()
        block_idx = current_block_idx[]
        if block_idx <= length(elliptic_params_obs[])
            params_list = elliptic_params_obs[]
            
            # Parse values from GUI
            max_iter = tryparse(Int, controls[:max_iter].stored_string[])
            tol = tryparse(Float64, controls[:tolerance].stored_string[])
            omega = tryparse(Float64, controls[:omega].stored_string[])
            
            left_a = tryparse(Float64, controls[:forcing_left_a].stored_string[])
            left_b = tryparse(Float64, controls[:forcing_left_b].stored_string[])
            right_a = tryparse(Float64, controls[:forcing_right_a].stored_string[])
            right_b = tryparse(Float64, controls[:forcing_right_b].stored_string[])
            bottom_a = tryparse(Float64, controls[:forcing_bottom_a].stored_string[])
            bottom_b = tryparse(Float64, controls[:forcing_bottom_b].stored_string[])
            top_a = tryparse(Float64, controls[:forcing_top_a].stored_string[])
            top_b = tryparse(Float64, controls[:forcing_top_b].stored_string[])
            
            # Create updated parameters if all parsed successfully
            if !isnothing(max_iter) && !isnothing(tol) && !isnothing(omega) &&
               !isnothing(left_a) && !isnothing(left_b) &&
               !isnothing(right_a) && !isnothing(right_b) &&
               !isnothing(bottom_a) && !isnothing(bottom_b) &&
               !isnothing(top_a) && !isnothing(top_b)
                
                params_list[block_idx] = GridGeneration.EllipticParams(
                    max_iter = max_iter,
                    tol = tol,
                    ω = omega,
                    a_decay_left = left_a,
                    b_decay_left = left_b,
                    a_decay_right = right_a,
                    b_decay_right = right_b,
                    a_decay_bottom = bottom_a,
                    b_decay_bottom = bottom_b,
                    a_decay_top = top_a,
                    b_decay_top = top_b
                )
                
                elliptic_params_obs[] = params_list
            end
        end
    end
    
    # Attach update_params to all parameter textboxes
    on(controls[:max_iter].stored_string) do _; update_params(); end
    on(controls[:tolerance].stored_string) do _; update_params(); end
    on(controls[:omega].stored_string) do _; update_params(); end
    on(controls[:forcing_left_a].stored_string) do _; update_params(); end
    on(controls[:forcing_left_b].stored_string) do _; update_params(); end
    on(controls[:forcing_right_a].stored_string) do _; update_params(); end
    on(controls[:forcing_right_b].stored_string) do _; update_params(); end
    on(controls[:forcing_bottom_a].stored_string) do _; update_params(); end
    on(controls[:forcing_bottom_b].stored_string) do _; update_params(); end
    on(controls[:forcing_top_a].stored_string) do _; update_params(); end
    on(controls[:forcing_top_b].stored_string) do _; update_params(); end
    
    # When blocks change (after split/solve), update smoothing block selector options
    on(generated_blocks_obs) do blocks
        if !isnothing(blocks) && isa(blocks, Vector) && !isempty(blocks)
            num_blocks = length(blocks)
            block_options = ["Block $i" for i in 1:num_blocks]
            controls[:block_selector].options[] = block_options
            
            # Ensure elliptic_params has the right number of entries
            current_params = elliptic_params_obs[]
            if length(current_params) < num_blocks
                # Add default params for new blocks
                new_params = copy(current_params)
                for i in (length(current_params)+1):num_blocks
                    push!(new_params, GridGeneration.EllipticParams())
                end
                elliptic_params_obs[] = new_params
            elseif length(current_params) > num_blocks
                # Trim extra params
                elliptic_params_obs[] = current_params[1:num_blocks]
            end
        end
    end
end

"""
    setup_split_block_selector!(controls, initialGrid)

Set up the split block selector to show initial grid blocks and update
the I-splits and J-splits labels with the index ranges of the selected block.
"""
function setup_split_block_selector!(controls, initialGrid)
    # When split block selector changes, update the labels with index ranges
    on(controls[:split_block_selector].selection) do selection
        block_str = string(selection)
        block_num = parse(Int, split(block_str, " ")[end])
        
        if block_num <= length(initialGrid)
            block = initialGrid[block_num]
            ni = size(block, 2)  # i-direction size
            nj = size(block, 3)  # j-direction size
            
            # Update labels with valid index ranges (1 to ni-1 for splits, etc.)
            controls[:i_splits_label][] = "I-splits (2 to $(ni-1))"
            controls[:j_splits_label][] = "J-splits (2 to $(nj-1))"
        end
    end
    
    # Initialize with Block 1 ranges
    if !isempty(initialGrid)
        block = initialGrid[1]
        ni = size(block, 2)
        nj = size(block, 3)
        controls[:i_splits_label][] = "I-splits (2 to $(ni-1))"
        controls[:j_splits_label][] = "J-splits (2 to $(nj-1))"
    end
end
