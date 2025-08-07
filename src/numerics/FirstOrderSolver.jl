using DifferentialEquations, BoundaryValueDiffEq


# solve First Order System
function SolveFirstOrderLinearSystem(M, Mx, N, x0, x1)

    
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


