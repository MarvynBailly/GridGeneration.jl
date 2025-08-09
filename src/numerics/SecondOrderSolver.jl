function SolveSecondOrder(f, x0, x1; N=100, omega=0.5, max_iter=100, tol=1e-8, verbose=false)
    function normInf(x)
        return maximum(abs, x)
    end
    

    h = 1.0 / (N + 1)
    x = range(0, 1, length=N+2)




    u = [x0; collect(x[2:end-1]); x1]  # initial guess: linear


    u_new = similar(u)

    # Storage for tridiagonal matrix
    a = fill(1 / h^2, N)    # subdiagonal
    b = fill(-2 / h^2, N)   # diagonal
    c = fill(1 / h^2, N)    # superdiagonal
    rhs = zeros(N)
    resNorm = 0

    for iter in 1:max_iter
        for i in 1:N
            fi = f(u[i+1])

            du = u[i+2] - u[i]
            
            rhs[i] = - (fi / (4h^2)) * du^2
        end
        rhs[1] += - x0 / h^2
        rhs[end] += - x1 / h^2


        # Solve linear system
        δu = ThomasAlg(a, b, c, rhs)

        # Update interior points
        u_new[1] = x0
        u_new[end] = x1
        u_new[2:N+1] = δu

        # Check convergence
        resNorm = normInf(u_new - u)
        if resNorm < tol
            if verbose @info("Converged in $iter iterations with norm $(resNorm)") end
            return x, u_new, resNorm
        end

        # if verbose println("Iteration $iter: norm = $(resNorm)") end

        # u .= u_new

        # add relaxation factor
        # omega = 0.5
        u = (1 - omega) * u + omega * u_new
    end

    @info("Did not converge after $max_iter iterations with norm $(resNorm).")
    return x, u_new, resNorm
end

function ThomasAlg(a, b, c, d)
    n = length(d)
    cp = similar(c)
    dp = similar(d)

    # forward sweep
    cp[1] = c[1] / b[1]
    dp[1] = d[1] / b[1]

    for i in 2:n
        denom = b[i] - a[i] * cp[i-1]
        cp[i] = c[i] / denom
        dp[i] = (d[i] - a[i] * dp[i-1]) / denom
    end

    # backward substitution
    x = zeros(n)
    x[end] = dp[end]
    for i = n-1:-1:1
        x[i] = dp[i] - cp[i] * x[i+1]
    end

    return x
end
