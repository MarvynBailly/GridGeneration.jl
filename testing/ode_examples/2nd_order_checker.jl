using Plots

include("../../src/GridGeneration.jl")




function M_func(x, scale, problem)
    base = 0.05
    if problem == 1
        return scale 
    elseif problem == 2
        return scale * (1 + 15 * (x))^(-2)
    elseif problem == 3
        return scale * (1 + 15 * (1-x))^(-2)
    elseif problem == 4
        return scale * exp(-(x -0.5)^2 / base)
    elseif problem == 5
        return scale * (1 + 15 * (x))^(-2) + scale * (1 + 15 * (1-x))^(-2)
    elseif problem == 6
        return scale * exp(-(x) / base)
    end
end

function M_u1_func(x, scale, problem)
    base = 0.05
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
    elseif problem == 6
        return  scale * exp(-(x) / base ) * (-1/base)
    end
end


function f(x, scale, problem)
    # return (M_u1_func(x, scale, problem) / (2 * M_func(x, scale, problem)))
    # return ( / (2 * M_func(x, scale, problem)))
    # return x^2 #(M_u1_func(x, scale, problem) / (2 * M_func(x, scale, problem)))
end



######################

problem = 4

name = ["Uniform", "x=0", "x=1", "x=0.5", "Edges", "Other"][problem]

folder = "ODENumericalMethods"
path = "docs/src/assets/images/$folder/"
saveFig = false

x0 = 0.0
x1 = 1.0
scale = 40000
omega = 0.5


method = "2ndOrder-omega=$omega"

forcing = x -> f(x, scale, problem)


# Ns = [5, 50, 100, 500]
Ns = [100]
norms = zeros(length(Ns))

for (i, N) in enumerate(Ns)
    x, u, resNorm = GridGeneration.SolveSecondOrder(forcing, x0, x1; N=N, omega=omega, tol=1e-15, max_iter=100, verbose=true)
    norms[i] = resNorm
end


N = Ns[end]
x, u, resNorm = GridGeneration.SolveSecondOrder(forcing, x0, x1; N=N, omega=omega, tol=1e-15, max_iter=1000, verbose=true)


x_solution = u
M_solution = M_func.(x, scale, problem)
M_u_solution = M_u1_func.(x, scale, problem)

p1 = plot(x, x_solution, label="x(s)", xlabel="s", ylabel="x(s)", title="Solution for x(s) vs s", linewidth=2)
plot!(p1, range(0, 1, length=N), range(0, 1, length=N), label="x=s", xlabel="s", ylabel="x(s)", title="Solution for x(s) vs s", linewidth=2)

p2 = plot(x, M_solution, label="M(x(s))", xlabel="s", ylabel="M(x(s))", title="M(x(s)) vs s", linewidth=2)

# plot!(p2, x, M_u_solution, label="M'(x(s))", xlabel="s", ylabel="M'(x(s))", title="M'(x(s)) vs s", linewidth=2)

p3 = scatter(x_solution, ones(length(x_solution)), label="grid points", xlabel="x(s)", ylabel="s", title="Distribution of x(s)", linewidth=2)

titlepanel = plot(title = "$name Clustering Metric (N = $N, method = $method)",
                framestyle = :none, grid = false, ticks = nothing)


                
p5 = plot(title = "Convergence Norms",
        xlabel = "N", ylabel = "Norm",
        legend = :topleft,
        linewidth = 2)
plot!(p5, Ns, norms, label="Norms", marker=:circle, yscale = :log10)

p6 = plot(title = "Forcing Function",
        xlabel = "x", ylabel = "f(x)",
        legend = :topleft,
        linewidth = 2)
    plot!(p6, x, forcing.(x), label="f(x)", xlabel="x", ylabel="f(x)", title="Forcing Function", linewidth=2)

p4 = plot(titlepanel, p1, p2, p3, p5, p6,
        layout = @layout([A{0.001h}; B C; D; E F]),
        size = (800, 650), legend = :topleft)

imageName = "$(name)_N$(N)_$(method).svg"
imagePath = "$(path)$(imageName)"

if saveFig
    savefig(p4, imagePath)
else 
    display(p4)
end