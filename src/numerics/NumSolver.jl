include("AnalyticSolver.jl")
include("SecondOrderSolver.jl")
include("CentralDiff.jl")
include("OptPointsSolver.jl")


"""
    SolveODE(m_func, xs; solver=:analytic, omega=0.5, max_iter=500, tol=1e-8)

Solve the equidistribution ODE x'' + M_x/(2M) x'^2 = 0 for the grid point locations, using the
metric function `m_func(x)` on the discrete points `xs`.
`solver` can be `:analytic` or `:numeric`; `omega`, `max_iter` and `tol` only apply to `:numeric`.
Returns the vector of point locations.
"""
function SolveODE(m_func, xs; solver=:analytic, omega=0.5, max_iter=500, tol=1e-8)
    if solver == :analytic
        return AnalyticalSolution(xs, m_func)
    elseif solver == :numeric
        return solve_numeric(m_func, xs; omega=omega, max_iter=max_iter, tol=tol)
    else
        throw(ArgumentError("Unknown solver :$solver (expected :analytic or :numeric)"))
    end
end

"""
    SolveODEFixedN(m_func, xs, N; solver=:analytic, omega=0.5, max_iter=500, tol=1e-8)

Same as [`SolveODE`](@ref), but the resulting solution has `N` points.
"""
function SolveODEFixedN(m_func, xs, N; solver=:analytic, omega=0.5, max_iter=500, tol=1e-8)
    xsFixed = range(0, xs[end], length=N)
    m_funcFixed = LinearInterpolate(xsFixed, m_func.(xsFixed))

    if solver == :analytic
        return AnalyticalSolution(xsFixed, m_funcFixed)
    elseif solver == :numeric
        return solve_numeric(m_funcFixed, xsFixed; omega=omega, max_iter=max_iter, tol=tol)
    else
        throw(ArgumentError("Unknown solver :$solver (expected :analytic or :numeric)"))
    end
end

function solve_numeric(m_func, xs; omega, max_iter, tol)
    mx = CentralDiff(m_func, xs)
    mx_func = LinearInterpolate(xs, mx)

    f = x -> mx_func(x) ./ (2 * m_func(x))

    sol, _ = SolveSecondOrder(f, xs; N=length(xs), omega=omega, max_iter=max_iter, tol=tol, verbose=false)
    return sol
end
