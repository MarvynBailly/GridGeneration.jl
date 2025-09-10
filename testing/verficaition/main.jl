using Plots
using LinearAlgebra

include("../../src/GridGeneration.jl")

function ComputeOptimalNumberofPoints(x, M)
    Nn = length(x)
    # set up computational coordinate
    s = range(0, 1, length=Nn)

    # compute x_s
    x_s = zeros(Nn)
    @inbounds for i in 2:Nn-1
        Δs = s[i+1] - s[i-1]
        x_s[i] = (x[i+1] - x[i-1]) / Δs
    end
    x_s[1] = (x[2] - x[1]) / (s[2] - s[1])               # forward difference
    x_s[end] = (x[end] - x[end-1]) / (s[end] - s[end-1]) # backward difference

    # compute integral with trapezoidal rule
    numer = 0.0
    denom = 0.0
    @inbounds for i in 2:Nn
        Δs   = s[i+1] - s[i]
        p    = 0.5*(M(x[i-1])*x_s[i-1]^2 + M(x[i])*x_s[i]^2)
        numer += p * Δs
        denom += (p^2) * Δs
    end
    # println("numer: $numer, denom: $denom")

    sigma_opt = sqrt(numer / denom)

    N_opt = floor(Int, 1/(sigma_opt))
    return N_opt
end

function prolongate(u_coarse, Nfine)
    Nc = length(u_coarse)
    x_coarse = range(0, 1; length=Nc)
    x_fine   = range(0, 1; length=Nfine)

    u_fine = similar(x_fine)

    j = 1
    for (i, xf) in enumerate(x_fine)
        # advance j so that x_coarse[j] <= xf <= x_coarse[j+1]
        while j < Nc && x_coarse[j+1] < xf
            j += 1
        end
        if xf ≤ x_coarse[1]
            u_fine[i] = u_coarse[1]
        elseif xf ≥ x_coarse[end]
            u_fine[i] = u_coarse[end]
        else
            xL, xR = x_coarse[j], x_coarse[j+1]
            uL, uR = u_coarse[j], u_coarse[j+1]
            θ = (xf - xL) / (xR - xL)
            u_fine[i] = (1-θ)*uL + θ*uR
        end
    end

    return u_fine
end



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
            if verbose @info("Converged in $iter iterations with norm $(resNorm[iter])") end
            resNorm = resNorm[1:iter]
            return u_new, resNorm
        end

        u = (1 - omega) * u + omega * u_new
    end

    return u, resNorm
end

function SolveFirstOrder(M, xs; N=100, omega=0.5, max_iter=100, tol=1e-8, verbose=false)
    # Boundary conditions
    x0, x1 = xs[1], xs[end]
    h = 1.0 / N
    Nint = N - 2
    s = range(0, 1, length=N)

    # Initial guess: just copy xs (linear or given)
    u = collect(xs)
    u_new = similar(u)

    resNorm = zeros(max_iter)

    for iter in 1:max_iter
        # compute M_s from current u
        M_x = zeros(N)
        x_s = zeros(N)
        @inbounds for i in 2:N-1
            Δs = s[i+1] - s[i-1]
            x_s[i] = (u[i+1] - u[i-1]) / Δs
            M_x[i] = (M(u[i+1]) - M(u[i-1])) / (u[i+1] - u[i-1])
        end
        M_s = M_x .* x_s
        f = M_s ./ (2 * M.(u))

        # set up tridiagonal matrix
        a = zeros(Nint)
        b = zeros(Nint)
        c = zeros(Nint)
        rhs = zeros(Nint)

        @inbounds for i in 1:Nint
            a[i] = 1 / h^2 - f[i] / (2 * h)
            b[i] = -2 / h^2
            c[i] = 1 / h^2 + f[i] / (2 * h)
        end

        rhs[1] += (f[2] / (2 * h) - 1 / h^2) * x0
        rhs[end] += (-1/h^2 - f[end-1] / (2 * h)) * x1

        # Solve linear system
        δu = GridGeneration.ThomasAlg(a, b, c, rhs)
        u_new[1], u_new[end] = x0, x1
        u_new[2:end-1] = δu

        # Convergence check
        resNorm[iter] = norm(u_new - u)
        if verbose
            @info "Iteration $iter: norm = $(resNorm[iter])"
        end
        if resNorm[iter] < tol
            return u_new, resNorm[1:iter]
        end

        # Relax update
        u .= (1 - omega) .* u .+ omega .* u_new
    end

    return u, resNorm
end



problem = 3

# domain
L = 5.0
N = 100
xs = range(0, L, length=N)
s = range(0, 1, length=N)

showTrue = true


if problem == 1
    k = 200^2
    M = x -> k
    Mx = x -> 0.0
    numOptTrue = sqrt(k) * L
    trueSol = s -> s * L
elseif problem == 2
    a = 10
    b = 200
    M = x -> (a + b * x)^2
    Mx = x -> 2 * b * (a + b * x)
    
    numOptTrue = 1/(2 * sqrt(1 / (L^2 * (2 * a + b * L)^2 ) ) )

    C = a * L + 0.5 * b * L^2
    trueSol = s -> (-a + sqrt(a^2 + 2 * b * C* s)) / b
elseif problem == 3
    b = 300
    e = 2.5
    M = x -> 100 + (b * (e - x))^2
    Mx = x -> -2 * (b^2) * (e - x)

    showTrue = false
end


f = x -> Mx(x) / (2 * M(x))

sol1, _ = SolveSecondOrder(f, xs; N=N, omega=0.02, max_iter=10000, tol=1e-8, verbose=false)

nOpt = ComputeOptimalNumberofPoints(sol1, M)
@info "Computed Optimal number of points: $nOpt"

# xs = prolongate(sol1, nOpt)

sol1, resNorm = SolveSecondOrder(f, xs; N=nOpt, omega=0.02, max_iter=10000, tol=1e-8, verbose=false)

nOpt = ComputeOptimalNumberofPoints(sol1, M)
@info "Computed Optimal number of points: $nOpt"




pltRes = plot(resNorm, yscale=:log10, xlabel="Iteration", ylabel="Residual Norm", title="Convergence History")
display(pltRes)


# sol2 = SolveFirstOrder(M, sol1; N=100, omega=0.02, max_iter=10000, tol=1e-8, verbose=true)



pM = plot(xs, M.(xs), label="M(x)", xlabel="x", ylabel="M", title="Metric Function")

pXvS = plot(s, sol1, label="Numerical Method 1", lw=2, xlabel="s", ylabel="x", title="Optimal Mapping x(s)")
# plot!(pXvS, s, sol2, label="Numerical Method 2", lw=2, legend=:topleft)

pPoints = scatter(
    sol1, zero.(sol1);
    markershape = :vline,
    markersize  = 10,         
    markerstrokewidth = 1.5,
    markerstrokecolor = :black,
    color = :black,
    legend = false,
    xlabel = "x", ylabel = "", yticks = false,
    title = "Generated Grid Points", 
    label="Numerical Method 1"
)
# scatter!(sol2, zero.(sol2); markershape=:vline, markersize=10, markerstrokewidth=1.5, markerstrokecolor=:red, color=:red, label="Numerical Method 2")



if showTrue
    @info "True Optimal number of points: $numOptTrue"

    plot!(pXvS, s, trueSol.(s), ls=:dash, label="True Solution", lw=2, legend=:topleft)
end

p = plot(pM, pXvS, pPoints, layout=(3,1), size=(800,600))
display(p)


