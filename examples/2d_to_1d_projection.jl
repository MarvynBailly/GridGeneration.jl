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
        M22 = scale * (1 .+ 15 * (y)).^(-2)
        return [M11, M22]
    end
end

function metricDerivative(x, y, scale, problem)
    if problem == 1
        M11 = -2 * scale * (1 .+ 15 * (x)).^(-3) * (15)
        M22 = -2 * scale * (1 .+ 15 * (y)).^(-3) * (15)
        return [M11, M22]
    end
end



function GetMetricValues(points, getMetric; method = "central")
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
    
    if method == "central"
        for i in 1:n
            # get metric value for the points M
            metricValues = getMetric(points[1, i], points[2, i])
            M = zeros(Float64, 2, 2)
            M[1, 1] = metricValues[1]
            M[2, 2] = metricValues[2]

            localDiff = diff[:, i]

            m_vals[i] = localDiff' * M * localDiff                        
        end
    end


    return m_vals
end


# Build M(x) as linear interpolation, Mx(x) as piecewise-constant slope.
function build_interps_linear(xv, Mv, Mxv)
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

    Mx_of = function (x::Real)
        if x ≤ xv[1]; return Mxv[1] end
        if x ≥ xv[end]; return Mxv[end] end
        i = searchsortedlast(xv, x)       # i such that xv[i] ≤ x < xv[i+1]
        t = (x - xv[i]) / Δx[i]
        return (1 - t)*Mxv[i] + t*Mxv[i+1]
    end

    return M_of, Mx_of
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
SectionIndices = 300:350
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

# build the metric
problem = 1
scale = 40000

M_func = (x,y) -> metric(x, y, scale, problem)
M_u1_func = (x,y) -> metricDerivative(x, y, scale, problem)
m = GetMetricValues(boundarySection, M_func)
mx = GetMetricValues(boundarySection, M_u1_func);


m_func, mx_func = build_interps_linear(xs, m, mx)

# pass the functions to the solver
sol = GridGeneration.SolveODE(m_func, mx_func, N, xs[1], xs[end])




p = scatter(sol[1, :], ones(size(sol[1, :])), label="x", color=:blue, markersize=2.5)
scatter!(p, xs, ones(size(xs)), label="x", color=:red, markersize=2)
display(p)