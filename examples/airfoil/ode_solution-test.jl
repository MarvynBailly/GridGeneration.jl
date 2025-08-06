include("../../src/GridGeneration.jl")
include("GetAirfoilGrid.jl")
# include("metric/Metric.jl")

using Plots, MAT, DelimitedFiles

function PlotProjectedPoints(projected_x, projected_y, projected_x_opt, projected_y_opt, x_bnd, y_bnd, metric, method, name)
    # 5. Plot the projected points on the airfoil
    p1 = plot(title = "Projected Points on Section (no ar)" )

    plot!(p1, x_bnd, y_bnd, 
        seriestype=:scatter, 
        label="Original Airfoil ($(length(x_bnd)))", 
        color=:black, 
        marker=:diamond,
        markersize=5, 
        markerstrokewidth=0)


    plot!(p1, x_bnd, y_bnd, 
        label="Original Airfoil", 
        color=:black, 

    )

    plot!(p1, projected_x, projected_y, 
        seriestype=:scatter, 
        label="Projected Points non-Optimal ($(length(projected_x)))", 
        color=:green, 
        markersize=3, 
        markerstrokewidth=0,
        xlabel="x",
        ylabel="y")

    plot!(p1, projected_x_opt, projected_y_opt, 
        seriestype=:scatter, 
        label="Projected Points Optimal ($(length(projected_x_opt)))", 
        color=:red, 
        marker=:rect,
        markersize=3, 
        markerstrokewidth=0,
        title="Projected Points on Airfoil Section",
        xlabel="x",
        ylabel="y")


    p2 = plot(projected_x_opt, projected_y_opt, 
        seriestype=:scatter, 
        label="Projected Points Optimal ($(length(projected_x_opt)))", 
        marker=:rect,
        color=:red, 
        markersize=2, 
        markerstrokewidth=0,
        title="Projected Points on Airfoil Section",
        xlabel="x",
        ylabel="y",
        aspect_ratio = 1)

    plot!(p2, x_bnd, y_bnd, 
        label="Original Airfoil", 
        color=:black, 
    )

    p3 = plot( title = "Discrete Boundary Metric Values", aspect_ratio = 1)
    plot!(p3, x_bnd, y_bnd, label="Boundary", color=:black, linewidth=1)
    
    scatter!(p3, x_bnd, y_bnd, label="Discrete Metric Values ($(length(x_bnd)))", marker_z = metric, markersize=4, markerstrokewidth=0, c = :brg)
    
    scatter!(p3, projected_x_opt, projected_y_opt, label="Projected Points Optimal ($(length(projected_x_opt)))", color = :red, marker=:rect, markersize=2, markerstrokewidth=0)


    titlepanel = plot(title = "1D Metric (method = $method) with $name clustering",
                framestyle = :none, grid = false, ticks = nothing)

    p5 = plot(titlepanel, p1, p2, p3, layout = @layout([A{0.001h}; b; c d]), size = (1000, 800))


    return p5
end



"""
 Get the optimal solution for the ODE grid spacing problem.
 This function solves the ODE system for grid spacing and computes the optimal number of points based on the metric values.
 Then recomputes the solution with the optimal number of points.
"""

function GetOptimalSolution(m, mx, N, xs; method = "system of odes")
    m_func = GridGeneration.build_interps_linear(xs, m)
    mx_func = GridGeneration.build_interps_linear(xs, mx)
    sol = SolveODE(m_func, mx_func, N, xs[1], xs[end])





    return sol, sol
end



"""
    SolveODE(M, Mx, N, x0, x1; method = :numeric, verbose = false)

Numerical solver for the ODE grid spacing problem using DifferentialEquations.jl.
"""
function SolveODE(M, Mx, N, x0, x1; method = :numeric, verbose = false)
    sol = SolveLinearSystem(M, Mx, N, x0, x1)

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
    
    # scatter!(ptest, xs, m, markershape=:circle, markersize=4, markerstrokewidth=0, c = :black, label="Boundary Values")
    
    bvp = BVProblem(SpacingODE!, BoundaryConditions!, uGuess, sspan, p)
    
    sol = solve(bvp, Shooting(Tsit5()), saveat=s_grid)
    
    return sol
end



"""
    ComputeOptimalSpacing(x, M, s)
Compute the optimal grid spacing based on the metric values `M` at points `x` and spacing `s`.
return Ceil(Int, 1 / sigma_opt) as the optimal number of points.
"""
function ComputeOptimalNumberofPoints(x, M, s)
    @assert length(x) == length(M) == length(s) "x, M, s must have same length (x: $(length(x)), M: $(length(M)), s: $(length(s)))"
    N = length(x)
    @assert N ≥ 2 "need at least two points"

    numer = 0.0
    denom = 0.0
    @inbounds for i in 2:N-1
        Δs  = s[i+1] - s[i]
        @assert Δs > 0 "s must be strictly increasing"
        x_s = (x[i+1] - x[i-1]) / (2*Δs)
        Mc  = 0.5*(M[i] + M[i+1])          # cell-avg M
        p   = Mc * x_s^2
        numer += p * Δs
        denom += (p^2) * Δs
    end

    sigma_opt = sqrt(numer / denom)
    N_opt = ceil(Int, 1 / sigma_opt)

    return N_opt
end


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


#####################
# SET UP DOMAIN
#####################

initialGrid = GetAirfoilGrid(airfoilPath = "examples/airfoil/data/A-airfoil.txt", radius = 3)

bottom = initialGrid[:,:,1]
airfoil = bottom[:, 101:end-100]  
SectionIndices = 400:500
# SectionIndices = 1:length(airfoil[1, :])
boundarySection = airfoil[:, SectionIndices]

N = length(boundarySection[1, :])

# xs = GridGeneration.Boundary2Dto1D(boundarySection)

# build the metric
saveFig = false
method = "local"
folder = "PointProjection"
path = "docs/src/assets/images/$folder/"

scale = 40000

problems = [1] # 1: x=0, 2: x=1, 3: uniform
names = ["x=0", "x=1", "uniform"]

# for (name, problem) in zip(names, problems) 
problem = 1
name = names[problem]

# M_func = (x,y) -> Metric(x, y, scale, problem)
# M_u1_func = (x,y) -> MetricDerivative(x, y, scale, problem)

# m = GridGeneration.Get1DMetric(boundarySection, M_func, method = method)
# mx = GridGeneration.Get1DMetric(boundarySection, M_u1_func, method = method)
xs = range(0, 1, length=N)

m = M_func.(boundarySection[1,:],  problem + 1 )
mx = M_u1_func.(boundarySection[1,:], problem + 1)



p = plot(xs, m, title = "1D Metric (method = $method) with $name clustering",
        xlabel = "s", ylabel = "m(x(s))", label = "m(x(s))",
        legend = :topright, linewidth=2)
display(p)
readline()









sol_opt, sol = GetOptimalSolution(m, mx, N, xs)

x_sol, x_sol_opt = sol[1, :], sol_opt[1, :]

projected_x, projected_y = GridGeneration.InterpolateBoundaryManually(x_sol, boundarySection)
projected_x_opt, projected_y_opt = GridGeneration.InterpolateBoundaryManually(x_sol_opt, boundarySection)



minM = minimum(m)

p1 = plot(xs, m, title = "ODE Solution for 1D Metric (method = $method) with $name clustering",
        xlabel = "s", ylabel = "m(x(s))", label = "m(x(s))",
        legend = :topright, linewidth=2)
scatter!(p1, xs, m, markershape=:circle, markersize=4, markerstrokewidth=0, c = :black, label="Boundary Values")
scatter!(p1, xs, zeros(length(xs)), markershape=:circle, markersize=4, markerstrokewidth=0, c = :black, label="")
scatter!(p1, xs, minM/2 * ones(length(xs)), markershape=:circle, markersize=4, markerstrokewidth=0, c = :black, label="")

scatter!(p1, sol[1,:], minM/2 * ones(length(sol[1,:])), markershape=:circle, markersize=2.5, markerstrokewidth=0, c = :red, label="Non-Optimal Solution (N = $(length(sol[1,:])))")
scatter!(p1, sol_opt[1,:], minM/2 * zeros(length(sol_opt[1,:])), markershape=:circle, markersize=2, markerstrokewidth=0, c = :green, label="Optimal Solution (N = $(length(sol_opt[1,:])))")


p2 = plot( title = "Optimal Points Projected on Boundary", aspect_ratio = 1)

plot!(p2, boundarySection[1, :], boundarySection[2, :], label="Boundary", color=:black, linewidth=1)
scatter!(p2, boundarySection[1, :], boundarySection[2, :], label="1D Metric", marker_z = m, markersize=4, markerstrokewidth=0)
scatter!(p2, projected_x_opt, projected_y_opt, label="Projected Points Optimal ($(length(projected_x_opt)))", color = :green, marker=:circle, markersize=2, markerstrokewidth=0)
# scatter!(p3, boundarySection[1, :], boundarySection[2, :], label="Discrete Metric Values ($(length(boundarySection[1,:])))", marker_z = m_vals_section, markersize=2, markerstrokewidth=0, c = :brg)

# p3 = plot(p1, p2, layout = @layout([A{0.001h}; b c]), size = (1000, 600))
p3 = plot(p1, p2, layout = @layout([a ; b]), size = (1000, 600))
display(p3)
# end