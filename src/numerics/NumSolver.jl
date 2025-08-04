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