include("FirstOrderSolver.jl")
include("AnalyticSolver.jl")
include("SecondOrderSolver.jl")


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
Solve the ODE system for grid spacing using the metric values `M` and `Mx` at points `x0` and `x1`.
The method can be either `:numeric` for numerical solution or `:analytic` for semi-analytic solution.
Returns the solution as a 2D array where the first row is the x-coordinates and

"""
function SolveODE(M, Mx, N, xs; method = :numeric, verbose = false)
    if verbose @info("Solving ODE for grid spacing using method: $method...") end

    if method == :numeric
        m_func = GridGeneration.build_interps_linear(xs, M)
        mx_func = GridGeneration.build_interps_linear(xs, Mx)


        x0 = xs[1]
        x1 = xs[end]
        sol = SolveFirstOrderLinearSystem(m_func, mx_func, N, x0, x1)
        if verbose println("Solution found.") end
    end

    if method == :analytic
        s_vals, x_sol = SolveAnalytic(xs, M, N)
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

