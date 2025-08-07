function solve_bvp_fixed_point(f, x0, x1; N=100, max_iter=100, tol=1e-8)
    h = 1.0 / (N + 1)
    x = range(0, 1, length=N+2)
    u = [x0; collect(x[2:end-1]); x1]  # initial guess: linear

    u_new = similar(u)

    for iter in 1:max_iter
        rhs = zeros(N)

        for i in 2:N+1
            fi = f(u[i])
            du = u[i+1] - u[i-1]
            rhs[i-1] = - (fi / (4h^2)) * du^2
        end

        # Build tridiagonal matrix
        a = fill(1/h^2, N)    # subdiagonal
        b = fill(-2/h^2, N)   # diagonal
        c = fill(1/h^2, N)    # superdiagonal

        # Solve linear system
        δu = thomas_algorithm(a, b, c, rhs)

        # Update interior points
        u_new[1] = x0
        u_new[end] = x1
        u_new[2:N+1] = δu

        # Check convergence
        if norm(u_new - u, Inf) < tol
            println("Converged in $iter iterations")
            return x, u_new
        end

        u .= u_new
    end

    error("Did not converge after $max_iter iterations.")
end

function thomas_algorithm(a, b, c, d)
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
