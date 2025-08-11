include("../src/GridGeneration.jl")
using Plots

function M_func(x, scale, base, problem)
    if problem == 1
        return scale * x^2
    elseif problem == 2
        return scale * (1 + 15 * x)^(-2)
    elseif problem == 3
        return scale * (1 + 15 * (1-x))^(-2)
    elseif problem == 4
        return scale * exp(-(x - 0.5)^2 / base)
    elseif problem == 5
        return scale * (1 + 15 * (x))^(-2) + scale * (1 + 15 * (1-x))^(-2)
    elseif problem == 6
        return scale * (1 + 15 * (x))^(-2) + scale * (1 + 15 * (1-x))^(-2) + scale * exp(-(x - 0.5)^2 / base)
    end
end

function true_solution(x, scale, problem)
    if problem == 1
        return sqrt.(x)
    else
        return nothing
    end
end

function equidistribute_linear_inverse(x, m, N)
    @assert length(x) == length(m) "x and m must have same length"
    @assert isapprox(first(x), 0.0; atol=1e-12) && isapprox(last(x), 1.0; atol=1e-12)
    @assert issorted(x) "x must be strictly increasing"
    @assert minimum(m) >= 0 "m must be positive"

    n = length(x) - 1
    Δx = x[2:end] .- x[1:end-1]
    w  = sqrt.(max.(m))  # √m

    # cumulative trapezoid
    I = similar(x, Float64)
    I[1] = 0.0
    for i in 1:n
        I[i+1] = I[i] + 0.5*(w[i] + w[i+1]) * Δx[i]
    end
    Itot = I[end]

    # uniform s-grid and targets
    s = range(0.0, 1.0; length=N)
    t = Itot .* s

    # invert by linear interpolation on (I, x)
    x_nodes = similar(s, Float64)
    for (j, tj) in enumerate(t)
        # find i with I[i] <= tj <= I[i+1]
        i = searchsortedlast(I, tj)
        if i == length(I)         # tj == I[end]
            x_nodes[j] = x[end]
        elseif I[i+1] == I[i]     # flat segment (m ~ 0)
            println("here")
            x_nodes[j] = x[i]     # or x[i] + θ*Δx[i] if you want uniform spread
        else
            θ = (tj - I[i]) / (I[i+1] - I[i])
            x_nodes[j] = x[i] + θ * (x[i+1] - x[i])
        end
    end
    return x_nodes
end


######################## 
######################## 
######################## 

# problem = 6
# problems = [1,2,3,4,5,6]

scale = 4000
base = 0.01
N = 100
xs = range(0.0, 1.0; length=N)

saveFig = true




for problem in 2:6
m_vals = M_func.(xs, scale, base, problem)
true_sol = true_solution.(xs, scale, problem)

# sol = equidistribute_linear_inverse(xs, m_vals, N)

# N_opt = GridGeneration.ComputeOptimalNumberofPoints(xs, m_vals, sol)

# sol_opt1 = equidistribute_linear_inverse(xs, m_vals, N_opt)

# @info("Optimal N found $N_opt")


sol_opt2, sol = GridGeneration.GetOptimalSolution(m_vals, m_vals, N, xs; method="analytic")
sol = sol[1, :]
sol_opt = sol_opt2[1, :]

N_opt = size(sol_opt, 2)

p1 = plot(xs, m_vals, label="m(x)", xlabel="x", ylabel="m", title="Mass Distribution")
p2 = scatter(sol, zero.(sol), label="Non-optimal x(s) (N = $N)", xlabel="x", ylabel="s", title="Point Distribution", markerstrokewidth = 0) 
scatter!(p2, sol_opt, zero.(sol_opt), label="Optimal x(s)(N = $(N_opt))", xlabel="x", ylabel="s", title="Point Distribution", markerstrokewidth = 0, color =:red, markersize = 2) 

if true_sol[1] != nothing
    scatter!(p2, true_sol, zero.(true_sol), label="True Solution", color=:black, marker=:diamond, markerstrokewidth = 0, markersize = 3)

    err = abs.(true_sol .- sol)
    p3 = scatter(sol, err, label="Abs Error", xlabel="x", ylabel="Error", title="Error in Equidistributed Grid Points", linewidth=2)

    plots = [p1, p2, p3]
    plotFinal = plot(plots..., layout=@layout([a; b; c]), size=(800, 600))
else
    plots = [p1, p2]
    plotFinal = plot(plots..., layout=@layout([a; b]), size=(800, 600))
end



if saveFig
    folder = "ODENumericalMethods"
    path = "docs/src/assets/images/$folder/"
    filename = "semianalytic=$(problem)_N=$N.png"
    savefig(plotFinal, joinpath(path, filename))
else
    display(plotFinal)
end

end