using LinearAlgebra

function residual(u, f, h)
    ui = u[2:end-1]          # interior nodes
    up = u[3:end]            # i+1
    um = u[1:end-2]          # i-1

    term1 = (up .- 2 .* ui .+ um) ./ h^2
    term2 = f.(ui) .* (up .- um).^2 ./ (4h^2)

    return term1 .+ term2
end

function SolveSecondOrder(f, xs; N=100, omega=0.5, max_iter=100, tol=1e-8, verbose=false)
    x0 = xs[1]
    x1 = xs[end]

    h = 1.0 / N

    u = collect(xs)  # initial guess: linear
    u_new = similar(u)

    # Storage for tridiagonal matrix
    Nint = N - 2
    a = fill(1 / h^2, Nint)    # subdiagonal
    b = fill(-2 / h^2, Nint)   # diagonal
    c = fill(1 / h^2, Nint)    # superdiagonal
    rhs = zeros(Nint)
    resNorm = zeros(max_iter)

    for iter in 1:max_iter
        rhs = - f.(u[2:end-1]) .* ((u[3:end] - u[1:end-2]).^2) ./ (4 * h^2)

        rhs[1] += - x0 / h^2
        rhs[end] += - x1 / h^2


        # Solve linear system
        δu = GridGeneration.ThomasAlg(a, b, c, rhs)

        # Update interior points
        u_new[1] = x0
        u_new[end] = x1
        u_new[2:end-1] = δu

        # Check convergence
        res = residual(u, f, h)
        resNorm[iter] = norm(res)

        if verbose @info("Iteration $iter: norm = $(resNorm[iter])") end

        if resNorm[iter] < tol
            @info("Converged in $iter iterations with norm $(resNorm[iter])")
            resNorm = resNorm[1:iter]
            return u_new, resNorm
        end

        u = (1 - omega) * u + omega * u_new
    end

    @info("Didn't converge with last norm $(resNorm[end])")
    return u, resNorm
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
