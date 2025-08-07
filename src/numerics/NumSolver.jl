using DifferentialEquations, BoundaryValueDiffEq


"""
 Get the optimal solution for the ODE grid spacing problem.
 This function solves the ODE system for grid spacing and computes the optimal number of points based on the metric values.
 Then recomputes the solution with the optimal number of points.
"""
function GetOptimalSolution(m, mx, N, xs; method = :numeric)    

    sol = SolveODE(m, mx, N, xs; method=method)

    @assert length(sol[1, :]) == N "Solution length ($(length(sol[1, :]))) does not match expected number of points (N = $N)"

    N_opt = ComputeOptimalNumberofPoints(sol[1, :], m, xs)
    sol_opt = SolveODE(m, mx, N_opt, xs; method=method)


    @info("Optimal number of points: ", N_opt)
    return sol_opt, sol
end



"""
    SolveODE(M, Mx, N, x0, x1; method = :numeric, verbose = false)

Numerical solver for the ODE grid spacing problem using DifferentialEquations.jl.
"""
function SolveODE(M, Mx, N, xs; method = :numeric, verbose = false)
    if verbose @info("Solving ODE for grid spacing using method: $method...") end

    if method == :numeric
        m_func = GridGeneration.build_interps_linear(xs, M)
        mx_func = GridGeneration.build_interps_linear(xs, Mx)


        x0 = xs[1]
        x1 = xs[end]
        sol = SolveLinearSystem(m_func, mx_func, N, x0, x1)
        if verbose println("Solution found.") end
    end

    if method == :analytic
        s_vals, x_sol = semi_analytic_solution_from_data(xs, M, N)
        sol = zeros(2, N)
        sol[1, :] = x_sol
        if verbose println("Analytic solution found.") end
    end

    if verbose
        @info("ODE solution details:")
        @info("  Number of points: $N")
        @info("  Start point: $x0")
        @info("  End point: $x1")
        @info("  Method used: $method")
    end

    return sol
end

"""
    ComputeOptimalSpacing(x, M, s)
Compute the optimal grid spacing based on the metric values `M` at points `x` and spacing `s`.
return Ceil(Int, 1 / sigma_opt) as the optimal number of points.
"""
function ComputeOptimalNumberofPoints(x, M, s)
    @assert length(x) == length(M) == length(s) "x, M, s must have same length (x: $(length(x)), M: $(length(M)), s: $(length(s)))"
    N = length(x)
    @assert N ≥ 2 "need at least two points"

    numer = 0.0
    denom = 0.0
    @inbounds for i in 2:N-1
        Δs  = s[i+1] - s[i]
        @assert Δs > 0 "s must be strictly increasing"
        x_s = (x[i+1] - x[i-1]) / (2*Δs)
        Mc  = 0.5*(M[i] + M[i+1])          # cell-avg M
        p   = Mc * x_s^2
        numer += p * Δs
        denom += (p^2) * Δs
    end

    sigma_opt = sqrt(numer / denom)
    N_opt = ceil(Int, 1 / sigma_opt)

    return N_opt
end



#####################
# NUMERICAL SOLVER - System of ODEs
#####################


# solve First Order System
function SolveLinearSystem(M, Mx, N, x0, x1)

    
    # Define the ODE system
    function SpacingODE!(du, u, p, s)
        M_func, M_u1_func, sigma, _, _ = p
        u1, u2 = u

        M_u1 = M_u1_func(u1)
        M = M_func(u1)
        
        du[1] = u2
        du[2] = - (M_u1 * u2^2) / (2 * M)    
    end
    
    
    # Set the boundary conditions
    function BoundaryConditions!(residual, u, p, s)
        _, _, _, x0, x1 = p
        residual[1] = u[1][1] - x0
        residual[2] = u[end][1] - x1
    end
    
    
    sspan = (0.0, 1.0)
    
    uGuess(s) = [x0 + s * (x1 - x0), x1 - x0]
    
    s_grid = range(0, 1, length=N)
    sigma = 1/N
    
    p = (M, Mx, sigma, x0, x1)
    
    bvp = BVProblem(SpacingODE!, BoundaryConditions!, uGuess, sspan, p)

    sol = solve(bvp, Shooting(Tsit5()), saveat=s_grid)

    return sol
end



#####################
# NUMERICAL SOLVER - Semi-analytic solution
#####################

function make_linear_inverse_interpolator(x_vals, y_vals)
    @assert length(x_vals) == length(y_vals)
    @assert all(diff(y_vals) .> 0) "y_vals must be strictly increasing for invertibility"

    function Iinv(y_query)
        @assert y_query ≥ y_vals[1] && y_query ≤ y_vals[end] "Query out of bounds"

        # Binary search to find the interval [i, i+1] such that y_vals[i] ≤ y_query ≤ y_vals[i+1]
        low, high = 1, length(y_vals) - 1
        while low ≤ high
            mid = div(low + high, 2)
            if y_vals[mid] ≤ y_query ≤ y_vals[mid+1]
                # linear interpolation
                y1, y2 = y_vals[mid], y_vals[mid+1]
                x1, x2 = x_vals[mid], x_vals[mid+1]
                t = (y_query - y1) / (y2 - y1)
                return x1 + t * (x2 - x1)
            elseif y_query < y_vals[mid]
                high = mid - 1
            else
                low = mid + 1
            end
        end

        error("Value not found in data range.")
    end

    return Iinv
end

function semi_analytic_solution_from_data(x_vals, M_vals, N)
    @assert length(x_vals) == length(M_vals) "x_vals and M_vals must be the same length"
    @assert all(diff(x_vals) .>= 0) "x_vals must be increasing"

    # Step 1: Compute I(x) via trapezoidal integration of sqrt(M)
    sqrtM = sqrt.(M_vals)
    I_vals = zeros(length(x_vals))
    for i in 2:length(x_vals)
        Δx = x_vals[i] - x_vals[i-1]
        I_vals[i] = I_vals[i-1] + 0.5 * Δx * (sqrtM[i] + sqrtM[i-1])
    end

    I_total = I_vals[end]

    # Step 2: Interpolation of I(x) to get inverse: I(x) → x
    # Itp = LinearInterpolation(I_vals, x_vals, extrapolation_bc=Throw())
    Iinv = make_linear_inverse_interpolator(x_vals, I_vals)

    # Step 3: Compute x(s) from I(x(s)) = s * I_total
    s_vals = range(0, 1; length=N)
    x_of_s = [Iinv(s * I_total) for s in s_vals]

    return s_vals, x_of_s
end
