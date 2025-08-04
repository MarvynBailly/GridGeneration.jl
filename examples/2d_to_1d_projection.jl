include("../src/GridGeneration.jl")

using Plots, MAT, DelimitedFiles


function plotGrid!(p, X, Y, clr; sz = 1.0, skipV=1, skipH=1)
    for j in 1:skipV:size(X, 1)
        plot!(p, X[j, :], Y[j, :], color=clr, lw=sz, label=false)
    end

    for i in 1:skipH:size(X, 2)
        plot!(p, X[:, i], Y[:, i], color=clr, lw=sz, label=false)
    end
end


function metric(x,y, scale, problem)
    if problem == 1
        M11 = scale * (1 .+ 15 * (x)).^(-2)
        M22 = 0 #scale * (1 .+ 15 * (y)).^(-2)
        return [M11, M22]
    elseif problem == 2
        M11 = scale * (1 .+ 15 * (1 .- x)).^(-2)
        M22 = 0 #scale * (1 .+ 15 * (1 .- y)).^(-2)
        return [M11, M22]
    elseif problem == 3
        M11 = scale
        M22 = scale
        return [M11, M22]
    end
end

function metricDerivative(x, y, scale, problem)
    if problem == 1
        M11 = -2 * scale * (1 .+ 15 * (x)).^(-3) * (15)
        M22 = 0 # -2 * scale * (1 .+ 15 * (y)).^(-3) * (15)
        return [M11, M22]
    elseif problem == 2
        M11 = -2 * scale * (1 .+ 15 * (1 .- x)).^(-3) * (-15)
        M22 = 0 # -2 * scale * (1 .+ 15 * (1 .- y)).^(-3) * (-15)
        return [M11, M22]
    elseif problem == 3
        M11 = 0
        M22 = 0 # 0
        return [M11, M22]
    end
end



function GetMetricValues(points, getMetric; method = "local")
    #---------------------------
    # Description: Compute metric between points depending on method 
    # Input: 2xn array of points 
    # Output: 1xn array of metric 
    #---------------------------

    n = size(points, 2)
    m_vals = zeros(Float64, n)
    diff = zeros(Float64, 2, n)
    diff[:, 2:n-1] = points[:, 3:n] - points[:, 1:n-2]
    diff[:, n] = points[:, n] - points[:, n-1]
    diff[:,1] = points[:, 2] - points[:, 1]
    
    for i in 1:n
        # get metric value for the points M
        metricValues = getMetric(points[1, i], points[2, i])
        M = zeros(Float64, 2, 2)
        M[1, 1] = metricValues[1]
        M[2, 2] = metricValues[2]
        
        
        localDiff = diff[:, i]
        if method == "local"
            m_vals[i] = localDiff' * M * localDiff                        
        elseif method == "nuclear"
            m_vals[i] = M[1, 1] + M[2, 2]
        end
    end


    return m_vals
end


# Build M(x) as linear interpolation, Mx(x) as piecewise-constant slope.
function build_interps_linear(xv, Mv)
    @assert length(xv) == length(Mv) ≥ 2
    @assert issorted(xv)

    # precompute slopes on each interval
    Δx = diff(xv)

    # linear interpolation with clamping to [xv[1], xv[end]]
    M_of = function (x::Real)
        if x ≤ xv[1]; return Mv[1] end
        if x ≥ xv[end]; return Mv[end] end
        i = searchsortedlast(xv, x)       # i such that xv[i] ≤ x < xv[i+1]
        t = (x - xv[i]) / Δx[i]
        return (1 - t)*Mv[i] + t*Mv[i+1]
    end

    return M_of
end



function getPlots(m_vals_section, x_sol, x_sol_opt, method, name)    
    diff = x_sol[2:end] - x_sol[1:end-1]
    diff_opt = x_sol_opt[2:end] - x_sol_opt[1:end-1]
    s = LinRange(0, 1, length(diff))
    s_opt = LinRange(0,1, length(diff_opt))

    p5 = scatter(x_sol, zeros(size(x_sol)) .+ 0.01, label="Non-opt Solution", color=:blue, linewidth=1.5,  title= "Grid Spacing with real metric")
    scatter!(p5, x_sol_opt, zeros(size(x_sol_opt)), label="opt Solution ($(length(x_sol_opt)))", color=:black, linewidth=1.5)
    scatter!(p5, s, diff, label="Non-opt Spacing", color=:red, markersize=2, markerstrokewidth=0)
    scatter!(p5, s_opt, diff_opt, label="opt Spacing", color=:green, markersize=2, markerstrokewidth=0)

    p6 = scatter(s, title = "M values", m_vals_section, label="M(s)", color=:black, markersize=2, markerstrokewidth=0)

        titlepanel = plot(title = "1D Metric (method = $method) with $name clustering",
                    framestyle = :none, grid = false, ticks = nothing)

    p7 = plot(titlepanel, p5, p6, layout = @layout([A{0.001h}; b c]), size = (1000, 600),)

    return p7
end

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




function ComputeOptimalSpacing(x, M, s)
    @assert length(x) == length(M) == length(s) "x, M, s must have same length"
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
    return sqrt(numer / denom)
end

function plotBoundary(airfoil, boundarySection, m_vals, m_vals_section, method, name)
    p1 = plot(title = "Discrete Metric Values", aspect_ratio = 1)
    plot!(p1, airfoil[1, :], airfoil[2, :], label="Boundary", color=:black, linewidth=0.5)
    scatter!(p1, airfoil[1, :], airfoil[2, :], label="Airfoil", marker_z = m_vals, markersize=0.5, markerstrokewidth=0, c = :brg)

    p2 = plot( title = "Zoomed Leading Edge (no aspect_ratio)",
            xlims = (boundarySection[1, 1] - 0.1, boundarySection[1, end] + 0.1),
            ylims = (minimum(boundarySection[2, :]) - 0.1, maximum(boundarySection[2, :]) + 0.1))

    plot!(p2, airfoil[1, :], airfoil[2, :], label="Boundary", color=:black, linewidth=0.5)
    scatter!(p2, airfoil[1, :], airfoil[2, :], label="Airfoil", marker_z = m_vals, markersize=0.5, markerstrokewidth=0, c = :brg)


    p3 = plot( title = "Discrete Boundary Metric Values (no aspect_ratio)")
    plot!(p3, boundarySection[1, :], boundarySection[2, :], label="Boundary", color=:black, linewidth=1)
    scatter!(p3, boundarySection[1, :], boundarySection[2, :], label="Discrete Metric Values ($(length(boundarySection[1,:])))", marker_z = m_vals_section, markersize=2, markerstrokewidth=0, c = :brg)

    titlepanel = plot(title = "1D Metric (method = $method) with $name clustering",
                    framestyle = :none, grid = false, ticks = nothing)

    p4 = plot(titlepanel, p1, p2,  p3, layout = @layout([A{0.001h}; b c; d]), size = (1000, 600))

    return p4
end

function InterpolateBoundaryManually(x_spacing, boundary; align = :auto)
    # boundary: 2×N (x row, y row)
    x_bnd = collect(boundary[1, :])
    y_bnd = collect(boundary[2, :])
    n = length(x_bnd)
    @assert size(boundary, 1) == 2 && n ≥ 2 "boundary must be 2×N with N≥2"

    # Step 1: cumulative arc length along boundary order
    arc_lengths = zeros(eltype(x_bnd), n)
    @inbounds for i in 2:n
        dx = x_bnd[i] - x_bnd[i-1]
        dy = y_bnd[i] - y_bnd[i-1]
        arc_lengths[i] = arc_lengths[i-1] + hypot(dx, dy)
    end
    L = arc_lengths[end]
    @assert L > 0 "boundary has zero total length"

    # Step 2: map x_spacing to [0, L]
    xs = (x_spacing .- first(x_spacing)) ./ (last(x_spacing) - first(x_spacing)) .* L

    # Step 2.1: decide orientation
    # If xs starts nearer to L than 0, flip parameterization (s -> L - s).
    do_flip = align == :reverse ? true :
              align == :forward ? false :
              abs(xs[1] - L) < abs(xs[1] - 0)

    if do_flip
        xs = L .- xs
    end

    # Clamp to [0,L] to avoid tiny numerical overshoots
    @inbounds for j in eachindex(xs)
        xs[j] = clamp(xs[j], 0, L)
    end

    # Step 3: piecewise linear interpolation along arc length
    projected_x = similar(xs, eltype(x_bnd))
    projected_y = similar(xs, eltype(y_bnd))

    @inbounds for j in eachindex(xs)
        s = xs[j]
        i = searchsortedlast(arc_lengths, s)
        if i ≥ n
            projected_x[j] = x_bnd[end]
            projected_y[j] = y_bnd[end]
            continue
        elseif i == 0
            projected_x[j] = x_bnd[1]
            projected_y[j] = y_bnd[1]
            continue
        end

        s0 = arc_lengths[i]
        s1 = arc_lengths[i+1]
        denom = s1 - s0
        t = denom == 0 ? 0.0 : (s - s0) / denom 

        projected_x[j] = (1 - t)*x_bnd[i] + t*x_bnd[i+1]
        projected_y[j] = (1 - t)*y_bnd[i] + t*y_bnd[i+1]
    end

    return projected_x, projected_y
end


#####################
# SET UP DOMAIN
#####################

airfoilPath = "examples/airfoil/A-airfoil.txt"

# read the airfoil data
airfoilData = readdlm(airfoilPath, '\t', skipstart=1)

RADIUS = 3
vertN = 100
horzN = 100

# number of points in the airfoil
airfoilN = length(airfoilData[:, 1])

# set up a c grid around the provided inner boundary
boundary = GridGeneration.SetupDomain(
    airfoilData, 
    RADIUS, 
    vertN, 
    horzN;
    type = "cgrid"
)


initialGrid = GridGeneration.TFI_2D(boundary)
initialGrid = permutedims(cat(initialGrid...; dims=3), (3, 1, 2));

bottom = initialGrid[:,:,1]
airfoil = bottom[:, 101:end-100]  
SectionIndices = 100:250
# SectionIndices = 1:length(airfoil[1, :])
boundarySection = airfoil[:, SectionIndices]

arclength = sum(sqrt.(diff(boundarySection[1, :]).^2 .+ diff(boundarySection[2, :]).^2))
N = length(boundarySection[1, :])

# boundarySection is 2×N (rows: x,y; columns: points in order)
x = boundarySection[1, :]
y = boundarySection[2, :]

# segment lengths (N-1)
Δx = diff(x)
Δy = diff(y)
Δs = sqrt.(Δx.^2 .+ Δy.^2)

xs = [0.0; cumsum(Δs)]   # length N, xs[1]=0, xs[end]=arclength

# normalize 
xs = xs ./ xs[end]  # now xs is in [0, 1]


# build the metric
problem = 3
saveFig = true
method = "local" # "local" or "nuclear"



folder = "Mapping2Dto1D"
path = "docs/src/assets/images/$folder/"


scale = 80000

M_func = (x,y) -> metric(x, y, scale, problem)
M_u1_func = (x,y) -> metricDerivative(x, y, scale, problem)

M_func_values = M_func.(boundarySection[1, :], boundarySection[2, :])


if problem == 1
    name = "x=0"
elseif problem == 2
    name = "x=1"
elseif problem == 3
    name = "uniform"
end


m_vals = GetMetricValues(airfoil, M_func, method = method)
m_vals_section = GetMetricValues(boundarySection, M_func, method = method)

p = plotBoundary(airfoil, boundarySection, m_vals, m_vals_section, method, name)
imageName = "metricboundary_$(name)_$(method).svg"
imagePath = "$path$imageName"


if saveFig
    savefig(p, imagePath)
else
    display(p)
    readline()
end




m = GetMetricValues(boundarySection, M_func, method = method)
mx = GetMetricValues(boundarySection, M_u1_func, method = method)
m_func = build_interps_linear(xs, m)
mx_func = build_interps_linear(xs, mx)

# pass the functions to the solver
sol = GridGeneration.SolveODE(m_func, mx_func, N, xs[1], xs[end])

sigma_opt = ComputeOptimalSpacing(sol[1, :], m, xs)
N_opt = ceil(Int, 1 / sigma_opt)

@info("Optimal number of points: ", N_opt)

if N_opt < 2
    @info("Optimal number of points is less than 2, using N=3")
    N_opt = 10
end

sol_opt = GridGeneration.SolveODE(m_func, mx_func, N_opt, xs[1], xs[end])

x_sol = sol[1, :]
x_sol_opt = sol_opt[1, :]

p = getPlots(m, x_sol, x_sol_opt, method, name)

imageName = "pointsmetric_$(name)_$(method).svg"
imagePath = "$path$imageName"


if saveFig
    savefig(p, imagePath)
else
    display(p)
    readline()

end

# # scale the solutions back to the correct size
# x_sol = x_sol * xs[end]
# x_sol_opt = x_sol_opt * xs[end]

# projected_x, projected_y = InterpolateBoundaryManually(x_sol, boundarySection)
# projected_x_opt, projected_y_opt = InterpolateBoundaryManually(x_sol_opt, boundarySection)

# p = PlotProjectedPoints(projected_x, projected_y, projected_x_opt, projected_y_opt, boundarySection[1, :], boundarySection[2, :], m_vals_section, method, name)

# imageName = "result_$(name)_$(method).svg"
# imagePath = "$path$imageName"


# if saveFig
#     savefig(p, imagePath)
# else
#     display(p)
# end