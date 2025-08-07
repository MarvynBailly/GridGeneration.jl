include("../src/GridGeneration.jl")

using Plots

# println the directory of the script

folder = "ODENumericalMethods"
path = "docs/src/assets/images/$folder/"


scale = 1000
base = 0.05


function M_func(x, problem)
    
    if problem == 1
        return scale 
    elseif problem == 2
        return scale * (1 + 15 * (x))^(-2)
    elseif problem == 3
        return scale * (1 + 15 * (1-x))^(-2)
    elseif problem == 4
        return scale * exp(-(x - 0.5)^2 / base)
    elseif problem == 5
        return scale * (1 + 15 * (x))^(-2) + scale * (1 + 15 * (1-x))^(-2)
    end
end

function M_u1_func(x, problem)
    
    if problem == 1
        return 0
    elseif problem == 2
        return -2 * scale * (1 + 15 * (x))^(-3) * (15)
    elseif problem == 3
        return -2 * scale * (1 + 15 * (1-x))^(-3) * (-15)
    elseif problem == 4
        return scale * exp(-(x - 0.5)^2 / base ) * (-2 * (x - 0.5) / base)
    elseif problem == 5
        return -2 * scale * (1 + 15 * (x))^(-3) * (15) + -2 * scale * (1 + 15 * (1-x))^(-3) * (-15)
    end
end


problems = [5];
# problem = 3 # 1: uniform, 2: clustering at x=0, 3: clustering at x=1, 4: clustering at x=0.5

method = :numeric # :numeric or :analytic

names = ["Uniform", "x=0", "x=1", "x=0.5", "Edges"]
scale = 1000

for problem in problems
    N = 500
    x0 = 0.0
    x1 = 1.0



    M(x) = M_func(x, problem)
    M_u1(x) = M_u1_func(x, problem)

    xs = range(x0, x1, length=N)

    m = M_func.(xs, problem)
    mx = M_u1_func.(xs, problem)

    sol = GridGeneration.SolveODE(m, mx, N, xs, method=method);

    imageName = "$(name)_N$(N)_$(method).svg"
    imagePath = "$(path)$(imageName)"

    x_solution = sol[1, :]
    M_solution = M.(x_solution)

    # plot the solution
    p1 = plot(range(0, 1, length=N), x_solution, label="x(s)", xlabel="s", ylabel="x(s)", title="Solution for x(s) vs s", linewidth=2)
    plot!(p1, range(0, 1, length=N), range(0, 1, length=N), label="x=s", xlabel="s", ylabel="x(s)", title="Solution for x(s) vs s", linewidth=2)

    p2 = plot(range(0, 1, length=N), M_solution, label="M(x(s))", xlabel="s", ylabel="M(x(s))", title="M(x(s)) vs s", linewidth=2)

    p3 = scatter(x_solution, ones(length(x_solution)), label="grid points", xlabel="x(s)", ylabel="s", title="Distribution of x(s)", linewidth=2)

    titlepanel = plot(title = "$name Clustering Metric (N = $N, method = $method)",
                    framestyle = :none, grid = false, ticks = nothing)

                    
    p4 = plot(titlepanel, p1, p2, p3,
            layout = @layout([A{0.001h}; B C; D]),
            size = (800, 650), legend = :topleft)

    display(p4)
    # savebig(p4, imagePath)
end