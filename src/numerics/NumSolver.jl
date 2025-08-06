"""
 Get the optimal solution for the ODE grid spacing problem.
 This function solves the ODE system for grid spacing and computes the optimal number of points based on the metric values.
 Then recomputes the solution with the optimal number of points.
"""

function GetOptimalSolution(m, mx, N, xs; method = "system of odes")
    m_func = GridGeneration.build_interps_linear(xs, m)
    mx_func = GridGeneration.build_interps_linear(xs, mx)


    if method == "system of odes"
        sol = GridGeneration.SolveODE(m_func, mx_func, N, xs[1], xs[end])

        @assert length(sol[1, :]) == N "Solution length ($(length(sol[1, :]))) does not match expected number of points (N = $N)"

        N_opt = GridGeneration.ComputeOptimalNumberofPoints(sol[1, :], m, xs)
        sol_opt = GridGeneration.SolveODE(m_func, mx_func, N_opt, xs[1], xs[end])
    end
    
    @info("Optimal number of points: ", N_opt)
    return sol_opt, sol
end



"""
    SolveODE(M, Mx, N, x0, x1; method = :numeric, verbose = false)

Numerical solver for the ODE grid spacing problem using DifferentialEquations.jl.
"""
function SolveODE(M, Mx, N, x0, x1; method = :numeric, verbose = false)
    if verbose @info("Solving ODE for grid spacing using method: $method...") end

    if method == :numeric
        sol = SolveLinearSystem(M, Mx, N, x0, x1)
        if verbose println("Solution found.") end
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

#####################
# NUMERICAL SOLVER
#####################
using DifferentialEquations, BoundaryValueDiffEq


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