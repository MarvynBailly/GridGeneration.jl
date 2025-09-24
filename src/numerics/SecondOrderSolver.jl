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
    h = 1.0 / N

    u = collect(xs)  # initial guess: linear
    u_new = similar(u)
    
    resNorm = zeros(max_iter)

    # Storage for tridiagonal matrix
    a = zeros(N) 
    b = zeros(N)
    c = zeros(N)
    rhs = zeros(N)
    u_s = similar(u)
    
    # main diagonal
    b[2:end-1] .= -2 / h^2      
    b[1] = 1.0
    b[end] = 1.0

    rhs[1] = xs[1]
    rhs[end] = xs[end]

    for iter in 1:max_iter
        u_s[2:end-1] = (u[3:end] - u[1:end-2]) / (2 * h)
        u_s[1] = (u[2] - u[1]) / h
        u_s[end] = (u[end] - u[end-1]) / h


        a[2:end-1] .= (1 / h^2) .- f.(u[2:end-1]) .* u_s[2:end-1] * (1 / (2 * h))   

        c[2:end-1] .= (1 / h^2) .+ f.(u[2:end-1]) .* u_s[2:end-1] * (1 / (2 * h))   

        # Solve linear system
        δu = GridGeneration.ThomasAlg(a, b, c, rhs)

        # Update interior points
        u_new = δu

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
