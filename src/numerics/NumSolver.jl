include("AnalyticSolver.jl")
include("SecondOrderSolver.jl")
include("CentralDiff.jl")
include("OptPointsSolver.jl")


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