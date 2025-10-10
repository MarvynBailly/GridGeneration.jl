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
Compute forcing terms for bottom/top wall to enforce orthogonality and spacing.
wall=1 for bottom (j=1), wall=Nj for top (j=Nj).
"""
function compute_forcing_eta(x::Matrix, y::Matrix, a_decay::Float64, b_decay::Float64; wall::Int = 1)
    @assert wall == size(x,2) || wall == 1 "wall must be either 1 (bottom) or Nj (top)"
    
    Ni, Nj = size(x)
    RHS_x = zeros(Ni)
    RHS_y = zeros(Ni)
    RHS_x_full = zeros(Ni, Nj)
    RHS_y_full = zeros(Ni, Nj)

    dir = wall == 1 ? 1 : -1
    
    # Calculate initial spacing at wall ends
    s1 = sqrt((x[1, wall] - x[1, wall + dir*1])^2 + (y[1, wall] - y[1, wall + dir*1])^2)
    s2 = sqrt((x[Ni, wall] - x[Ni, wall + dir*1])^2 + (y[Ni, wall] - y[Ni, wall + dir*1])^2)
    s_vec = LinRange(s1, s2, Ni)
    
    # ξ-derivatives along boundary
    x_xi = (x[3:Ni, wall] - x[1:Ni-2, wall]) / 2.0
    y_xi = (y[3:Ni, wall] - y[1:Ni-2, wall]) / 2.0
    x_xixi = x[3:Ni, wall] - 2*x[2:Ni-1, wall] + x[1:Ni-2, wall]
    y_xixi = y[3:Ni, wall] - 2*y[2:Ni-1, wall] + y[1:Ni-2, wall]

    # Impose orthogonality and spacing to find desired η-derivatives
    ds_inv = 1.0 ./ sqrt.(x_xi.^2 + y_xi.^2)
    x_eta_desired = -s_vec[2:Ni-1] .* y_xi .* ds_inv
    y_eta_desired = s_vec[2:Ni-1] .* x_xi .* ds_inv

    # Compute η second derivatives using one-sided second-order formula
    x_etaeta = 0.5*(-7*x[2:Ni-1, wall] + 8*x[2:Ni-1, wall+dir*1] - x[2:Ni-1, wall+dir*2]) - dir*3*x_eta_desired
    y_etaeta = 0.5*(-7*y[2:Ni-1, wall] + 8*y[2:Ni-1, wall+dir*1] - y[2:Ni-1, wall+dir*2]) - dir*3*y_eta_desired

    # Boundary metrics
    alpha_b = x_eta_desired.^2 + y_eta_desired.^2
    gamma_b = x_xi.^2 + y_xi.^2

    # Compute final RHS for governing equations at boundary
    RHS_x[2:Ni-1] = -(alpha_b .* x_xixi + gamma_b .* x_etaeta)
    RHS_y[2:Ni-1] = -(alpha_b .* y_xixi + gamma_b .* y_etaeta)

    # Propagate boundary forcing terms into domain with exponential decay
    for j in 1:Nj
        eta_dist = abs(j - wall) # Distance from the wall
        RHS_x_full[:, j] = 1000 * RHS_x .* exp(-a_decay * eta_dist)
        RHS_y_full[:, j] = 1000 * RHS_y .* exp(-b_decay * eta_dist)
    end

    return RHS_x_full, RHS_y_full
end

"""
Compute forcing terms for left/right wall to enforce orthogonality and spacing.
wall=1 for left (i=1), wall=Ni for right (i=Ni).
"""
function compute_forcing_xi(x::Matrix, y::Matrix, a_decay::Float64, b_decay::Float64; wall::Int = 1)
    @assert wall == 1 || wall == size(x, 1) "wall must be either 1 (left) or Ni (right)"
    
    Ni, Nj = size(x)
    RHS_x_boundary = zeros(Nj)
    RHS_y_boundary = zeros(Nj)
    RHS_x_full = zeros(Ni, Nj)
    RHS_y_full = zeros(Ni, Nj)

    dir = wall == 1 ? 1 : -1
    
    # Calculate initial spacing at wall ends
    s1 = sqrt((x[wall, 1] - x[wall + dir*1, 1])^2 + (y[wall, 1] - y[wall + dir*1, 1])^2)
    s2 = sqrt((x[wall, Nj] - x[wall + dir*1, Nj])^2 + (y[wall, Nj] - y[wall + dir*1, Nj])^2)
    s_vec = LinRange(s1, s2, Nj)

    for j in 2:Nj-1
        s = s_vec[j]

        # η-derivatives along boundary
        x_eta   = (x[wall, j+1] - x[wall, j-1]) / 2.0
        y_eta   = (y[wall, j+1] - y[wall, j-1]) / 2.0
        x_etaeta = x[wall, j+1] - 2*x[wall, j] + x[wall, j-1]
        y_etaeta = y[wall, j+1] - 2*y[wall, j] + y[wall, j-1]

        # Impose orthogonality and spacing to find desired ξ-derivatives
        ds_inv = 1.0 / sqrt(x_eta^2 + y_eta^2)
        x_xi_desired = -s * y_eta * ds_inv
        y_xi_desired = s * x_eta * ds_inv

        # Compute ξ second derivatives using one-sided second-order formula
        x_xixi = 0.5*(-7*x[wall, j] + 8*x[wall+dir*1, j] - x[wall+dir*2, j]) + dir*3*x_xi_desired
        y_xixi = 0.5*(-7*y[wall, j] + 8*y[wall+dir*1, j] - y[wall+dir*2, j]) + dir*3*y_xi_desired

        # Boundary metrics
        alpha_b = s^2
        gamma_b = x_eta^2 + y_eta^2

        # Compute final RHS for governing equations at boundary
        RHS_x_boundary[j] = -(alpha_b * x_xixi + gamma_b * x_etaeta)
        RHS_y_boundary[j] = -(alpha_b * y_xixi + gamma_b * y_etaeta)
    end

    # Propagate boundary forcing terms into domain with exponential decay
    for i in 1:Ni
        xi_dist = abs(i - wall)
        RHS_x_full[i, :] = RHS_x_boundary .* exp(-a_decay * xi_dist)
        RHS_y_full[i, :] = RHS_y_boundary .* exp(-b_decay * xi_dist)
    end

    return RHS_x_full, RHS_y_full
end



"""
Main iterative solver for elliptic grid generation using SOR.
"""
function EllipticSolver(x::Matrix, y::Matrix; params)
    function verbose_print(msg)
        if params.verbose
            println(msg)
        end
    end

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

        x_old = copy(x)
        
        # Calculate metric terms based on current grid
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


