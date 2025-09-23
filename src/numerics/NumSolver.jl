include("AnalyticSolver.jl")
include("SecondOrderSolver.jl")
include("CentralDiff.jl")
include("OptPointsSolver.jl")

"""
 Get the optimal solution for the ODE grid spacing problem.
 This function solves the ODE system for grid spacing and computes the optimal number of points based on the metric values.
 Then recomputes the solution with the optimal number of points.
"""
# function GetOptimalSolution(m, mx, N, xs; dir=1, method = "analytic")    

#     sol = SolveODE(m, mx, N, xs; dir=dir, method=method)

#     @assert length(sol[1, :]) == N "Solution length ($(length(sol[1, :]))) does not match expected number of points (N = $N)"

    
#     N_opt = ComputeOptimalNumberofPoints(xs, m)
#     @info("Optimal number of points: ", N_opt)

#     # if N_opt == 1 
#     #     N_opt = 3
#     # end

#     # if N_opt % 2 == 0
#     #     N_opt += 1
#     # end

#     sol_opt = SolveODE(m, mx, N_opt, xs; dir=dir, method=method)

#     return sol_opt, sol
# end



"""
    SolveODE(M, xs; method = :numeric, verbose = false)
Solve the ODE x'' + M_x/(2M) x' = 0 for grid spacing using function M(x) on discrete points xs. 
method can be :analytic or :numeric.
Returns the solution and the residual norm (res only for numeric).
"""

function SolveODE(m_func, xs; solver=:analytic)
    if solver == :analytic
        sol = GridGeneration.AnalyticalSolution(xs, m_func)
        return sol
    elseif solver == :numeric
        N = length(xs)
        mx = GridGeneration.CentralDiff(m_func, xs)
        mx_func = GridGeneration.LinearInterpolate(xs, mx)

        f = x -> mx_func(x) ./ (2 * m_func(x))

        sol, resNorm = GridGeneration.SolveSecondOrder(f, xs; N=N, omega=0.02, max_iter=10, tol=1e-5, verbose=false)
        return sol, resNorm
    end
end

"""
    SolveODE(M, xs, N; method = :numeric, verbose = false)
Solve the ODE x'' + M_x/(2M) x' = 0 for grid spacing using function M(x) on discrete points xs. Resulting solution will have N points. 
method can be :analytic or :numeric.
Returns the solution and the residual norm (res only for numeric).
"""
function SolveODEFixedN(m_func, xs, N; solver=:analytic)
    xsFixed = range(0, xs[end], length=N)
    m_funcFixed = GridGeneration.LinearInterpolate(xsFixed, m_func.(xsFixed))
    
    if solver == :analytic
        sol = GridGeneration.AnalyticalSolution(xsFixed, m_funcFixed)
        return sol
    elseif solver == :numeric
        mx = GridGeneration.CentralDiff(m_funcFixed, xsFixed)
        mx_func = GridGeneration.LinearInterpolate(xsFixed, mx)

        f = x -> mx_func(x) ./ (2 * m_funcFixed(x))

        sol, resNorm = GridGeneration.SolveSecondOrder(f, xsFixed; N=N, omega=0.02, max_iter=10, tol=1e-5, verbose=false)
        return sol, resNorm
    end
end