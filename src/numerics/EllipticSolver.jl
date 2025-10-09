using Plots


"""
Calculates the metric terms α, β, γ and the Jacobian J.
Uses central differences for the interior and one-sided for boundaries.
"""
function calculate_metrics(x::Matrix, y::Matrix)
    Ni, Nj = size(x)
    
    x_xi  = zeros(Ni, Nj); y_xi  = zeros(Ni, Nj)
    x_eta = zeros(Ni, Nj); y_eta = zeros(Ni, Nj)

    for j in 1:Nj
        for i in 1:Ni
            # ξ-derivatives
            if i == 1
                x_xi[i, j] = x[i+1, j] - x[i, j]
                y_xi[i, j] = y[i+1, j] - y[i, j]
            elseif i == Ni
                x_xi[i, j] = x[i, j] - x[i-1, j]
                y_xi[i, j] = y[i, j] - y[i-1, j]
            else
                x_xi[i, j] = (x[i+1, j] - x[i-1, j]) / 2.0
                y_xi[i, j] = (y[i+1, j] - y[i-1, j]) / 2.0
            end
        end
    end

    for i in 1:Ni
        for j in 1:Nj
             # η-derivatives
            if j == 1
                x_eta[i, j] = x[i, j+1] - x[i, j]
                y_eta[i, j] = y[i, j+1] - y[i, j]
            elseif j == Nj
                x_eta[i, j] = x[i, j] - x[i, j-1]
                y_eta[i, j] = y[i, j] - y[i, j-1]
            else
                x_eta[i, j] = (x[i, j+1] - x[i, j-1]) / 2.0
                y_eta[i, j] = (y[i, j+1] - y[i, j-1]) / 2.0
            end
        end
    end
    
    # Metric terms
    alpha = x_eta.^2 + y_eta.^2
    beta  = x_xi .* x_eta + y_xi .* y_eta
    gamma = x_xi.^2 + y_xi.^2
    J     = x_xi .* y_eta - x_eta .* y_xi

    return alpha, beta, gamma, J
end

"""
Computes the forcing terms at the wall boundary (j=N)
to enforce orthogonality and spacing. 
"""
function compute_forcing_eta(x::Matrix, y::Matrix, a_decay::Float64, b_decay::Float64; wall::Int = 1)
    @assert wall == size(x,2) || wall == 1 "wall must be either 1 (bottom) or Nj (top)"
    
    Ni, Nj = size(x)
    RHS_x = zeros(Ni)
    RHS_y = zeros(Ni)

    RHS_x_full = zeros(Ni, Nj)
    RHS_y_full = zeros(Ni, Nj)

    # Determine the sign for the desired η-derivative based on the wall
    dir = wall == 1 ? 1 : -1
    
    
    # 1. Calculate the true initial spacing at both ends of the boundary (Euclidean distance)
    s1 = sqrt((x[1, wall] - x[1, wall + dir*1])^2 + (y[1, wall] - y[1, wall + dir*1])^2)
    s2 = sqrt((x[Ni, wall] - x[Ni, wall + dir*1])^2 + (y[Ni, wall] - y[Ni, wall + dir*1])^2)
    
    # 2. Linearly interpolate between s1 and s2 to get a smooth spacing vector
    s_vec = LinRange(s1, s2, Ni)
    

    x_xi = (x[3:Ni, wall] - x[1:Ni-2, wall]) / 2.0
    y_xi = (y[3:Ni, wall] - y[1:Ni-2, wall]) / 2.0
    x_xixi = x[3:Ni, wall] - 2*x[2:Ni-1, wall] + x[1:Ni-2, wall]
    y_xixi = y[3:Ni, wall] - 2*y[2:Ni-1, wall] + y[1:Ni-2, wall]

    # 2. Impose orthogonality and spacing to find desired η-derivatives
    # let's use x - and y + coordinates and account for it in the next computation
    ds_inv = 1.0 ./ sqrt.(x_xi.^2 + y_xi.^2)
    x_eta_desired = -s_vec[2:Ni-1] .* y_xi .* ds_inv
    y_eta_desired = s_vec[2:Ni-1] .* x_xi .* ds_inv

    # 3. Compute second derivatives normal to the wall using the initial grid
    # and the one-sided, second-order accurate formula.
    # note that we account for the direction of the wall (top or bottom) using `dir` here
    x_etaeta = 0.5*(-7*x[2:Ni-1, wall] + 8*x[2:Ni-1, wall+dir*1] - x[2:Ni-1, wall+dir*2]) - dir*3*x_eta_desired
    y_etaeta = 0.5*(-7*y[2:Ni-1, wall] + 8*y[2:Ni-1, wall+dir*1] - y[2:Ni-1, wall+dir*2]) - dir*3*y_eta_desired

    # define boundary metrics
    alpha_b = x_eta_desired.^2 + y_eta_desired.^2 #s_vec[2:Ni-1].^2 # Since alpha = x_eta^2 + y_eta^2
    gamma_b = x_xi.^2 + y_xi.^2

    # 5. Compute the final RHS for the governing equations at the boundary
    RHS_x[2:Ni-1] = -(alpha_b .* x_xixi + gamma_b .* x_etaeta)
    RHS_y[2:Ni-1] = -(alpha_b .* y_xixi + gamma_b .* y_etaeta)


    # plot the compoments of the RHS_x
    # p = plot(1:Ni, RHS_x, label="RHS_x", title="RHS_x = -(α x_xixi + γ x_etaeta) along wall j=$wall", xlabel="i", ylabel="RHS_x =")
    # plot!(alpha_b, label="alpha_b")
    # plot!(x_xixi, label="x_xixi")
    # plot!(gamma_b, label="gamma_b")
    # plot!(x_etaeta, label="x_etaeta")
    # display(p)
    # readline()


    # for i in 2:Ni-1
    #     s = s_vec[i]
    #     # s = min(s1, s2)

    #     # 1. Derivatives along the boundary (ξ-derivatives)
    #     x_xi   = (x[i+1, wall] - x[i-1, wall]) / 2.0
    #     y_xi   = (y[i+1, wall] - y[i-1, wall]) / 2.0
    #     x_xixi = x[i+1, wall] - 2*x[i, wall] + x[i-1, wall]
    #     y_xixi = y[i+1, wall] - 2*y[i, wall] + y[i-1, wall]

    #     # 2. Impose orthogonality and spacing to find desired η-derivatives
    #     # let's use x - and y + coordinates and account for it in the next computation
    #     ds_inv = 1.0 / sqrt(x_xi^2 + y_xi^2)
    #     x_eta_desired = -s * y_xi * ds_inv
    #     y_eta_desired = s * x_xi * ds_inv

    #     # 3. Compute second derivatives normal to the wall using the initial grid
    #     # and the one-sided, second-order accurate formula.
    #     # note that we account for the direction of the wall (top or bottom) using `dir` here
    #     x_etaeta = 0.5*(-7*x[i,wall] + 8*x[i,wall+dir*1] - x[i,wall+dir*2]) - dir*3*x_eta_desired
    #     y_etaeta = 0.5*(-7*y[i,wall] + 8*y[i,wall+dir*1] - y[i,wall+dir*2]) - dir*3*y_eta_desired


    #     # 4. Define boundary metrics. By construction, β=0.
    #     alpha_b = s^2 # Since alpha = x_eta^2 + y_eta^2
    #     gamma_b = x_xi^2 + y_xi^2

    #     # 5. Compute the final RHS for the governing equations at the boundary
    #     # This is the term that forces the grid to the desired state.
    #     RHS_x[i] = -(alpha_b * x_xixi + gamma_b * x_etaeta)
    #     RHS_y[i] = -(alpha_b * y_xixi + gamma_b * y_etaeta)
    # end



    # Propagate the boundary forcing terms into the domain with exponential decay
    for j in 1:Nj
        eta_dist = abs(j - wall) # Distance from the wall
        RHS_x_full[:, j] = 1000 * RHS_x .* exp(-a_decay * eta_dist)
        RHS_y_full[:, j] = 1000 * RHS_y .* exp(-b_decay * eta_dist)
    end

    return RHS_x_full, RHS_y_full
end

"""
Computes the forcing terms at the wall boundary (i=1 or N)
to enforce orthogonality and spacing. 
"""
function compute_forcing_xi(x::Matrix, y::Matrix, a_decay::Float64, b_decay::Float64; wall::Int = 1)
    @assert wall == 1 || wall == size(x, 1) "wall must be either 1 (left) or Ni (right)"

    # struggles when the wall is close to orthogonal. Look into cost method?
    
    Ni, Nj = size(x)
    RHS_x_boundary = zeros(Nj)
    RHS_y_boundary = zeros(Nj)
    RHS_x_full = zeros(Ni, Nj)
    RHS_y_full = zeros(Ni, Nj)

    # Determine the sign for the desired η-derivative based on the wall
    dir = wall == 1 ? 1 : -1
    
    
    # compute s from the edge of the wall on both sides
    s1 = abs(x[wall, 1] - x[wall + dir*1, 1])
    s2 = abs(x[wall, Nj] - x[wall + dir*1, Nj])
    # take s as the min
    s = min(s1, s2)


    s1 = sqrt((x[wall, 1] - x[wall + dir*1, 1])^2 + (y[wall, 1] - y[wall + dir*1, 1])^2)
    s2 = sqrt((x[wall, Nj] - x[wall + dir*1, Nj])^2 + (y[wall, Nj] - y[wall + dir*1, Nj])^2)
    
    # 2. Linearly interpolate between s1 and s2 to get a smooth spacing vector
    s_vec = LinRange(s1, s2, Nj)


    # interpolate s1 and s2. Then look up the s value at each i location
    for j in 2:Nj-1
        s = s_vec[j]
        # s = min(s1, s2)

        # 1. Derivatives along the boundary (η-derivatives)
        x_eta   = (x[wall, j+1] - x[wall, j-1]) / 2.0
        y_eta   = (y[wall, j+1] - y[wall, j-1]) / 2.0
        x_etaeta = x[wall, j+1] - 2*x[wall, j] + x[wall, j-1]
        y_etaeta = y[wall, j+1] - 2*y[wall, j] + y[wall, j-1]

        # 2. Impose orthogonality and spacing to find desired η-derivatives
        # let's use x - and y + coordinates and account for it in the next computation
        ds_inv = 1.0 / sqrt(x_eta^2 + y_eta^2)
        x_xi_desired = -s * y_eta * ds_inv
        y_xi_desired = s * x_eta * ds_inv

        # 3. Compute second derivatives normal to the wall using the initial grid
        # and the one-sided, second-order accurate formula.
        # note that we account for the direction of the wall (top or bottom) using `dir` here
        x_xixi = 0.5*(-7*x[wall, j] + 8*x[wall+dir*1, j] - x[wall+dir*2, j]) + dir*3*x_xi_desired
        y_xixi = 0.5*(-7*y[wall, j] + 8*y[wall+dir*1, j] - y[wall+dir*2, j]) + dir*3*y_xi_desired


        # 4. Define boundary metrics. By construction, β=0.
        alpha_b = s^2 # Since alpha = x_eta^2 + y_eta^2
        gamma_b = x_eta^2 + y_eta^2

        # 5. Compute the final RHS for the governing equations at the boundary
        # This is the term that forces the grid to the desired state.
        RHS_x_boundary[j] = -(alpha_b * x_xixi + gamma_b * x_etaeta)
        RHS_y_boundary[j] = -(alpha_b * y_xixi + gamma_b * y_etaeta)
    end

    # Propagate the boundary forcing terms into the domain with exponential decay
    for i in 1:Ni
        xi_dist = abs(i - wall) # Distance from the wall
        RHS_x_full[i, :] = RHS_x_boundary .* exp(-a_decay * xi_dist)
        RHS_y_full[i, :] = RHS_y_boundary .* exp(-b_decay * xi_dist)
    end

    return RHS_x_full, RHS_y_full
end



"""
Main iterative solver for the elliptic grid generation equations using SOR.
"""
function EllipticSolver(x::Matrix, y::Matrix; params)
    function verbose_print(msg)
        if params.verbose
            println(msg)
        end
    end

    # @info "Applying forcing terms at walls: Top: $useTopWall, Bottom: $useBottomWall, Left: $useLeftWall, Right: $useRightWall"

    Ni, Nj = size(x)
    error = 0
    finalIter = 0


    verbose_print("Starting SOR iterations...")
    for iter in 1:params.max_iter
        finalIter = iter
        RHS_x_full = zeros(Ni, Nj)
        RHS_y_full = zeros(Ni, Nj)

        if params.useTopWall
            RHS_x_top, RHS_y_top = compute_forcing_eta(x, y, params.a_decay_top, params.b_decay_top, wall = Nj)
            RHS_x_full .+= RHS_x_top
            RHS_y_full .+= RHS_y_top
        end
        if params.useBottomWall
            RHS_x_bottom, RHS_y_bottom = compute_forcing_eta(x, y, params.a_decay_bottom, params.b_decay_bottom, wall = 1)
            RHS_x_full .+= RHS_x_bottom
            RHS_y_full .+= RHS_y_bottom
            
        end
        if params.useLeftWall
            RHS_x_left, RHS_y_left = compute_forcing_xi(x, y, params.a_decay_left, params.b_decay_left, wall = 1)
            RHS_x_full .+= RHS_x_left
            RHS_y_full .+= RHS_y_left

        end
        if params.useRightWall
            RHS_x_right, RHS_y_right = compute_forcing_xi(x, y, params.a_decay_right, params.b_decay_right, wall = Ni)
            RHS_x_full .+= RHS_x_right
            RHS_y_full .+= RHS_y_right
        end
        
        # Create two heatmap plots
        # p1 = heatmap(RHS_x_full, title="Heatmap of RHS_x_full", c=:viridis, aspect_ratio=:equal)
        # p2 = heatmap(RHS_y_full, title="Heatmap of RHS_y_full", c=:viridis, aspect_ratio=:equal)

        # # Display them in a single figure with a 1x2 layout
        # p = plot(p1, p2, layout = (1, 2), size = (800, 400))
        # println("Displaying RHS heatmaps. Close the plot window to continue...")
        # println("Max of RHS_x_full: ", maximum(abs.(RHS_x_full)))
        # println("Max of RHS_y_full: ", maximum(abs.(RHS_y_full)))
        # display(p)
        # readline()

        x_old = copy(x)
        
        # Calculate metric terms based on the current grid
        alpha, beta, gamma, J = calculate_metrics(x, y)
        
        # Perform one SOR iteration on the interior points
        for j in 2:Nj-1
            for i in 2:Ni-1
                denom = 2.0 * (alpha[i,j] + gamma[i,j])
                
                # --- X Equation ---
                term_xixi_x = alpha[i,j] * (x[i+1,j] + x[i-1,j])
                term_etaeta_x = gamma[i,j] * (x[i,j+1] + x[i,j-1])
                term_xieta_x = 0.5 * beta[i,j] * (x[i+1,j+1] - x[i-1,j+1] - x[i+1,j-1] + x[i-1,j-1])
                
                x_new = (term_xixi_x + term_etaeta_x - term_xieta_x + RHS_x_full[i,j]) / denom
                x[i, j] = (1 - params.ω) * x[i, j] + params.ω * x_new

                # --- Y Equation ---
                term_xixi_y = alpha[i,j] * (y[i+1,j] + y[i-1,j])
                term_etaeta_y = gamma[i,j] * (y[i,j+1] + y[i,j-1])
                term_xieta_y = 0.5 * beta[i,j] * (y[i+1,j+1] - y[i-1,j+1] - y[i+1,j-1] + y[i-1,j-1])
                
                y_new = (term_xixi_y + term_etaeta_y - term_xieta_y + RHS_y_full[i,j]) / denom
                y[i, j] = (1 - params.ω) * y[i, j] + params.ω * y_new
            end
        end

        # Check for convergence
        error = norm(x - x_old)

        if iter % 500 == 0
            verbose_print("Iter: $iter, Error: $error")
        end
        if error < params.tol
            verbose_print("Convergence reached at iteration $iter with error $error.")
            break
        end
        if iter == params.max_iter
            verbose_print("Warning: Maximum iterations reached without convergence.")
        end
    end
    
    return x, y, error, finalIter
end


