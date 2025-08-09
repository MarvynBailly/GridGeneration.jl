include("../src/GridGeneration.jl")
include("airfoil/GetAirfoilGrid.jl")
include("airfoil/metric/Metric.jl")

using Plots, MAT, DelimitedFiles

function M_func(x, scale, problem)
    base = 0.1
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

function M_u1_func(x, scale, problem)
    base = 0.1
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


function f(x, scale, problem)
    return (M_u1_func(x, scale, problem) / (2 * M_func(x, scale, problem)))
    # return x^2 #(M_u1_func(x, scale, problem) / (2 * M_func(x, scale, problem)))
end

problem = 3
scale = 40000

met1 = x -> f(x, scale, problem) # M_func(x, scale, problem)




M_func_test = (x, y) -> Metric(x, 0, scale, problem - 1 )
M_u1_func_test = (x, y) -> MetricDerivative(x, 0, scale, problem - 1)

met2 = (x) -> (M_u1_func_test(x, 0) ./ (2 * M_func_test(x, 0)))[1]




s = range(0,1,length=100)

initialGrid = GetAirfoilGrid(airfoilPath = "examples/airfoil/data/A-airfoil.txt", radius = 3)
bottom = initialGrid[:,:,1]
airfoil = bottom[:, 101:end-100]  
sectionIndices = 100:300
boundarySection = airfoil[:, sectionIndices]
boundarySectionReverse = reverse(boundarySection, dims=2)


xs1 = GridGeneration.ProjectBoundary2Dto1D(boundarySection)
xs2 = GridGeneration.ProjectBoundary2Dto1D(boundarySectionReverse)
xs = xs1


# boom the issue is with get1dmetric
m = GridGeneration.Get1DMetric(reverse(boundarySection, dims=2), M_func_test)



mx = GridGeneration.Get1DMetric(boundarySection, M_u1_func_test)

m_func = GridGeneration.build_interps_linear(xs, m)
mx_func = GridGeneration.build_interps_linear(xs, mx)

#################
# TODO: The issue is here. This means that it must be either coming from build_inters 
met3 = x -> mx_func(x) ./ (2 * m_func(x))
#################

N = 100
omega = 0.5
x0 = 0
x1 = 1

x, u, resNorm = GridGeneration.SolveSecondOrder(met3, x0, x1; N=N, omega=omega, tol=1e-15, max_iter=200, verbose=true)

p = plot(range(0,1,length=length(m)), m, label="m", xlabel="x", ylabel="y", title="Metric Function", legend=:topleft)


first = [v[1] for v in M_func_test.(s, 0)]

plot!(p, s, first, label="M_func_test", xlabel="x", ylabel="y", title="Metric Function", linewidth=2)
display(p)


# p1 = plot(s, met2.(s), label="met2", xlabel="x", ylabel="y", title="Metric Function", legend=:topleft)
# plot!(p1, s, met3.(s), label="met3", xlabel="s", ylabel="y", title="Metric Function", linewidth=2)
# p2 = scatter(u, zeros(length(x)), label="Grid", xlabel="x", ylabel="y", title="Grid Generation", legend=:topleft, markershape=:circle, markersize=2)

# p3 = plot(x, u, label="x(s)", xlabel="s", ylabel="x(s)", title="Solution for x(s) vs s", linewidth=2)
# plot!(p3, x, x, label="x=s", xlabel="s", ylabel="x(s)", title="Solution for x(s) vs s", linewidth=2)


# finalPlot = plot(p1, p2, p3, layout=@layout([a;b;c]), size=(800, 600), title="Airfoil Grid Generation with Metric Function")

# display(finalPlot)