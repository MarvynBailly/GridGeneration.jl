using Plots
using LinearAlgebra

# Note on Indexing Convention:
# This code uses (i, j) indexing for matrices, where:
# i corresponds to the ξ-direction (along the channel, size Ni)
# j corresponds to the η-direction (across the channel, size Nj)
# So, arrays are dimensioned as x[Ni, Nj].

"""
Helper function to calculate the y-coordinates of a symmetric NACA 4-digit airfoil.
For NACA 0012, t=0.12.
"""
function naca_0012(x, c=1.0)
    # Thickness for NACA 0012 is 12% of the chord
    t = 0.12
    # Standard formula for NACA 4-digit airfoils (half-thickness)
    y_t = 5 * t * c * (0.2969 * sqrt(x/c) 
                     - 0.1260 * (x/c) 
                     - 0.3516 * (x/c)^2 
                     + 0.2843 * (x/c)^3 
                     - 0.1015 * (x/c)^4) # Note: -0.1036, -0.1015
    return y_t
end

function setup_airfoil(Ni::Int, Nj::Int)
    @assert isodd(Ni) "Ni must be an odd number to have a center point at the leading edge."

    x = zeros(Float64, Ni, Nj)
    y = zeros(Float64, Ni, Nj)

    # --- Parameters ---
    chord_length = 1.0
    far_field_radius = 5.0

    # --- 1. Define the boundaries ---
    
    # Bottom boundary (j=1): The airfoil surface
    # We use cosine clustering to get more points at the leading and trailing edges.
    Ni_half = (Ni + 1) ÷ 2
    x_dist = 0.5 * chord_length * (1.0 .- cos.(range(0, π, length=Ni_half)))
    
    # Upper surface (from leading edge to trailing edge)
    x[Ni_half:Ni, 1] .= x_dist
    y[Ni_half:Ni, 1] .= naca_0012.(x_dist, chord_length)
    
    # Lower surface (from trailing edge to leading edge)
    # Note: We reverse the arrays and skip the leading edge point which is already defined.
    x[1:Ni_half-1, 1] .= reverse(x_dist[1:end-1])
    y[1:Ni_half-1, 1] .= -naca_0012.(reverse(x_dist[1:end-1]), chord_length)

    # Top boundary (j=Nj): Far-field boundary (a semi-circle)
    # The semi-circle is centered around the airfoil's mid-chord
    theta = range(2*π, 0, length=Ni)
    x[:, Nj] .= chord_length / 2.0 .+ far_field_radius .* cos.(theta)
    y[:, Nj] .= far_field_radius .* sin.(theta)
    
    # Left (i=1) and Right (i=Ni) boundaries connect the trailing edge to the far-field.
    # We define them via linear interpolation.
    y[1, :] = range(y[1, 1], y[1, Nj], length=Nj)
    x[1, :] = range(x[1, 1], x[1, Nj], length=Nj)
    
    y[Ni, :] = range(y[Ni, 1], y[Ni, Nj], length=Nj)
    x[Ni, :] = x[1, :] # In a C-grid, the wake cut lines have the same x-coords

    # --- 2. Make an initial guess for the interior grid (Linear Interpolation) ---
    # This interpolates between the airfoil surface and the far-field boundary.
    for i in 2:Ni-1
        for j in 2:Nj-1
            # Linear interpolation factor
            eta_j = (j - 1) / (Nj - 1)
            
            # Interpolate between the airfoil (j=1) and far-field (j=Nj) boundaries
            x[i, j] = (1 - eta_j) * x[i, 1] + eta_j * x[i, Nj]
            y[i, j] = (1 - eta_j) * y[i, 1] + eta_j * y[i, Nj]
        end
    end
    
    return x, y
end


"""
Sets up the physical domain and creates an initial algebraic grid.
The domain is a channel with a sinusoidal bump on the bottom wall.
"""
function setup_geometry_and_initial_grid(Ni::Int, Nj::Int)
    x = zeros(Float64, Ni, Nj)
    y = zeros(Float64, Ni, Nj)


    # Physical domain extents
    x_min, x_max = 0.0, 10.0
    y_max = 5.0 # Channel height

    # 1. Define the boundaries
    xi_pts = range(x_min, x_max, length=Ni)
    
    # Bottom boundary (j=1): Flat plate with a sinusoidal bump
    for i in 1:Ni
        x_val = xi_pts[i]
        x[i, 1] = x_val
        # A smooth bump in the middle of the domain
        if 1.0 <= x_val <= 9.0
            y[i, 1] = 1.0 * (1 + cos(pi * (x_val - 2.0)))
        else
            y[i, 1] = 0.0
        end
    end

    # Top boundary (j=Nj)
    # x[:, Nj] = xi_pts
    for i in 1:Ni
        x_val = xi_pts[i]
        x[i, Nj] = x_val
        #a smooth bump
        if 1.0 <= x_val <= 9.0
            y[i, Nj] = y_max - 1.0 * (1 + cos(pi * (x_val - 2.0)))
        else
            y[i, Nj] = y_max
        end
    end

    # Left (i=1) and Right (i=Ni) boundaries are straight vertical lines
    y[1, :] = range(y[1, 1], y_max, length=Nj)
    x[1, :] .= x_min
    y[Ni, :] = range(y[Ni, 1], y_max, length=Nj)
    x[Ni, :] .= x_max

    # 2. Make an initial guess for the interior grid (Linear Interpolation)
    for i in 2:Ni-1
        for j in 2:Nj-1
            eta_j = (j - 1) / (Nj - 1)
            x[i, j] = (1 - eta_j) * x[i, 1] + eta_j * x[i, Nj]
            y[i, j] = (1 - eta_j) * y[i, 1] + eta_j * y[i, Nj]
        end
    end
    
    return x, y
end


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
function compute_forcing_eta(x::Matrix, y::Matrix, sold::Float64, a_decay::Float64, b_decay::Float64; wall::Int = 1)
    @assert wall == size(x,2) || wall == 1 "wall must be either 1 (bottom) or Nj (top)"
    
    Ni, Nj = size(x)
    RHS_x = zeros(Ni)
    RHS_y = zeros(Ni)

    RHS_x_full = zeros(Ni, Nj)
    RHS_y_full = zeros(Ni, Nj)

    # Determine the sign for the desired η-derivative based on the wall
    dir = wall == 1 ? 1 : -1
    
    
    # compute s from the edge of the wall on both sides
    s1 = abs(y[1, wall] - y[1, wall + dir*1])
    s2 = abs(y[Ni, wall] - y[Ni, wall + dir*1])
    # take s as the min 
    s = min(s1, s2)
    # println("Computing forcing at wall j=$wall with dir=$dir using s=$s" )
    

    # interpolate s1 and s2. Then look up the s value at each i location

    for i in 2:Ni-1
        # 1. Derivatives along the boundary (ξ-derivatives)
        x_xi   = (x[i+1, wall] - x[i-1, wall]) / 2.0
        y_xi   = (y[i+1, wall] - y[i-1, wall]) / 2.0
        x_xixi = x[i+1, wall] - 2*x[i, wall] + x[i-1, wall]
        y_xixi = y[i+1, wall] - 2*y[i, wall] + y[i-1, wall]

        # 2. Impose orthogonality and spacing to find desired η-derivatives
        # let's use x - and y + coordinates and account for it in the next computation
        ds_inv = 1.0 / sqrt(x_xi^2 + y_xi^2)
        x_eta_desired = -s * y_xi * ds_inv
        y_eta_desired = s * x_xi * ds_inv

        # 3. Compute second derivatives normal to the wall using the initial grid
        # and the one-sided, second-order accurate formula.
        # note that we account for the direction of the wall (top or bottom) using `dir` here
        x_etaeta = 0.5*(-7*x[i,wall] + 8*x[i,wall+dir*1] - x[i,wall+dir*2]) - dir*3*x_eta_desired
        y_etaeta = 0.5*(-7*y[i,wall] + 8*y[i,wall+dir*1] - y[i,wall+dir*2]) - dir*3*y_eta_desired


        # 4. Define boundary metrics. By construction, β=0.
        alpha_b = s^2 # Since alpha = x_eta^2 + y_eta^2
        gamma_b = x_xi^2 + y_xi^2

        # 5. Compute the final RHS for the governing equations at the boundary
        # This is the term that forces the grid to the desired state.
        RHS_x[i] = -(alpha_b * x_xixi + gamma_b * x_etaeta)
        RHS_y[i] = -(alpha_b * y_xixi + gamma_b * y_etaeta)
    end

 

    # Propagate the boundary forcing terms into the domain with exponential decay
    for j in 1:Nj
        eta_dist = abs(j - wall) # Distance from the wall
        RHS_x_full[:, j] = RHS_x .* exp(-a_decay * eta_dist)
        RHS_y_full[:, j] = RHS_y .* exp(-b_decay * eta_dist)
    end

    return RHS_x_full, RHS_y_full
end

"""
Computes the forcing terms at the wall boundary (i=1 or N)
to enforce orthogonality and spacing. 
"""
function compute_forcing_xi(x::Matrix, y::Matrix, sold::Float64, a_decay::Float64, b_decay::Float64; wall::Int = 1)
    @assert wall == 1 || wall == size(x, 1) "wall must be either 1 (left) or Ni (right)"
    
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
    # println("Computing forcing at wall j=$wall with dir=$dir using s=$s" )
    

    # interpolate s1 and s2. Then look up the s value at each i location

    for j in 2:Nj-1
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
function EllipticSolver(; x::Matrix, y::Matrix,
                        max_iter::Int=2000, tol::Float64=1e-6, ω::Float64=1.7, 
                        s_left::Float64=0.02, a_decay_left::Float64=0.5, b_decay_left::Float64=0.5, 
                        s_right::Float64=0.02, a_decay_right::Float64=0.5, b_decay_right::Float64=0.5,
                        s_top::Float64=0.02, a_decay_top::Float64=0.5, b_decay_top::Float64=0.5,
                        s_bottom::Float64=0.02, a_decay_bottom::Float64=0.5, b_decay_bottom::Float64=0.5,
                        use_top_wall::Bool=true, use_bottom_wall::Bool=true, use_left_wall::Bool=true, use_right_wall::Bool=true,
                        verbose::Bool=false
                              )
    function verbose_print(msg)
        if verbose
            println(msg)
        end
    end

    @info "Applying forcing terms at walls: Top: $use_top_wall, Bottom: $use_bottom_wall, Left: $use_left_wall, Right: $use_right_wall"

    Ni, Nj = size(x)
    error = 0
    finalIter = 0


    verbose_print("Starting SOR iterations...")
    for iter in 1:max_iter
        finalIter = iter
        RHS_x_full = zeros(Ni, Nj)
        RHS_y_full = zeros(Ni, Nj)
        
        if use_top_wall
            RHS_x_top, RHS_y_top = compute_forcing_eta(x, y, s_top, a_decay_top, b_decay_top, wall = Nj)
            RHS_x_full .+= RHS_x_top
            RHS_y_full .+= RHS_y_top
        end
        if use_bottom_wall
            RHS_x_bottom, RHS_y_bottom = compute_forcing_eta(x, y, s_bottom, a_decay_bottom, b_decay_bottom, wall = 1)
            RHS_x_full .+= RHS_x_bottom
            RHS_y_full .+= RHS_y_bottom
        end
        if use_left_wall
            RHS_x_left, RHS_y_left = compute_forcing_xi(x, y, s_left, a_decay_left, b_decay_left, wall = 1)
            RHS_x_full .+= RHS_x_left
            RHS_y_full .+= RHS_y_left

        end
        if use_right_wall
            RHS_x_right, RHS_y_right = compute_forcing_xi(x, y, s_right, a_decay_right, b_decay_right, wall = Ni)
            RHS_x_full .+= RHS_x_right
            RHS_y_full .+= RHS_y_right
        end
        
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
                x[i, j] = (1 - ω) * x[i, j] + ω * x_new

                # --- Y Equation ---
                term_xixi_y = alpha[i,j] * (y[i+1,j] + y[i-1,j])
                term_etaeta_y = gamma[i,j] * (y[i,j+1] + y[i,j-1])
                term_xieta_y = 0.5 * beta[i,j] * (y[i+1,j+1] - y[i-1,j+1] - y[i+1,j-1] + y[i-1,j-1])
                
                y_new = (term_xixi_y + term_etaeta_y - term_xieta_y + RHS_y_full[i,j]) / denom
                y[i, j] = (1 - ω) * y[i, j] + ω * y_new
            end
        end

        # Check for convergence
        error = norm(x - x_old)

        if iter % 50 == 0
            verbose_print("Iter: $iter, Error: $error")
        end
        if error < tol
            verbose_print("Convergence reached at iteration $iter with error $error.")
            break
        end
        if iter == max_iter
            verbose_print("Warning: Maximum iterations reached without convergence.")
        end
    end
    
    return x, y, error, finalIter
end


