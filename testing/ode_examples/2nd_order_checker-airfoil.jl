using Plots, MAT

include("../../src/GridGeneration.jl")



include("../ode_examples/airfoil/metric/CustomMetric.jl")
include("airfoil/GetAirfoilGrid.jl")

include("../../plotter/metric_grid_plotter.jl")
include("../../plotter/blocks_interfaces_boundaries.jl")


function compute_Ms_over_M(xs, m_vals)
    N = length(xs)
    ratio = zeros(N)

    for i in 2:N-1
        ds = xs[i+1] - xs[i-1]
        ms = (m_vals[i+1] - m_vals[i-1]) / ds
        ratio[i] = ms / m_vals[i]
    end

    # one-sided at boundaries
    ratio[1]   = (m_vals[2] - m_vals[1]) / ( (xs[2]-xs[1]) * m_vals[1] )
    ratio[end] = (m_vals[end] - m_vals[end-1]) / ( (xs[end]-xs[end-1]) * m_vals[end] )

    return ratio
end


function SolveSecondOrder_MS(xs, f; x0=0.0, x1=1.0)
    n = length(xs)
    @assert n ≥ 2
    # ratio = compute_Ms_over_M(xs, m_vals)          # ≈ (M_s/M) at nodes

    a = zeros(n); b = zeros(n); c = zeros(n); rhs = zeros(n)

    # Dirichlet BCs
    b[1] = 1.0; rhs[1] = x0
    b[n] = 1.0; rhs[n] = x1

    # Interior rows (non-uniform)
    for i in 2:n-1
        hi   = xs[i]   - xs[i-1]
        hip1 = xs[i+1] - xs[i]
        Hi   = hi + hip1
        ai   =  f[i]#0.5 * ratio[i]                 # a_i = (1/2) * (M_s/M)

        a[i] = 2/(hi*Hi)   - ai/Hi
        b[i] = -2/(hi*hip1)
        c[i] = 2/(hip1*Hi) + ai/Hi
    end

    u = GridGeneration.ThomasAlg(a, b, c, rhs)
    return u
end

function Metric(x,y; problem = 1, scale = 4000)
    """
    - problem 1: Constant metric
    - problem 2: leading edge
    - problem 3: trailing edge
    - problem 4: leading edge and trailing edge
    - problem 5: for fun
    - problem 6: real metric
    """
    if problem == 1
        M11 = scale
        M22 = scale
        return [M11, M22]
    elseif problem == 2
        metricFunction =  make_getMetric(airfoil;
                A_airfoil = 100,  ℓ_airfoil = 0.5, p_airfoil = 2,   
                A_origin  = scale,  ℓ_origin  = 0.1, p_origin  = 10,   
                floor     = 1e-4,  origin_center=(0, 0),
                profile   = :rational
            ) 
        return metricFunction(x,y)
    elseif problem == 3
        metricFunction =  make_getMetric(airfoil;
                A_airfoil = 100.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
                A_origin  = 500.0,  ℓ_origin  = 0.1, p_origin  = 10,   
                floor     = 1e-4,  origin_center=(1, 0),
                profile   = :rational
            ) 
        return metricFunction(x,y)
    elseif problem == 4
        metricFunction1 =  make_getMetric(airfoil;
                A_airfoil = 100.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
                A_origin  = 500.0,  ℓ_origin  = 0.1, p_origin  = 10,   
                floor     = 1e-4,  origin_center=(0, 0),
                profile   = :rational
            ) 
        metricFunction2 =  make_getMetric(airfoil;
                A_airfoil = 100.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
                A_origin  = 500.0,  ℓ_origin  = 0.1, p_origin  = 10,   
                floor     = 1e-4,  origin_center=(1, 0),
                profile   = :rational
            )

        return metricFunction1(x,y) .+ metricFunction2(x,y)
    elseif problem == 5
        metricFunction1 =  make_getMetric(airfoil;
                A_airfoil = 100.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
                A_origin  = 1000.0,  ℓ_origin  = 0.1, p_origin  = 10,   
                floor     = 1e-4,  origin_center=(0, 0),
                profile   = :rational
            ) 
        metricFunction2 =  make_getMetric(airfoil;
                A_airfoil = 100.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
                A_origin  = 1000.0,  ℓ_origin  = 0.1, p_origin  = 10,   
                floor     = 1e-4,  origin_center=(1, 0),
                profile   = :rational
            )
        metricFunction3 =  make_getMetric(airfoil;
                A_airfoil = 100.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
                A_origin  = 1000.0,  ℓ_origin  = 0.1, p_origin  = 10,   
                floor     = 1e-4,  origin_center=(0.30, 0.13),
                profile   = :rational
            ) 
        metricFunction4 =  make_getMetric(airfoil;
                A_airfoil = 100.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
                A_origin  = 1000.0,  ℓ_origin  = 0.1, p_origin  = 10,   
                floor     = 1e-4,  origin_center=(0.50, -0.09),
                profile   = :rational
            )

        return metricFunction1(x,y) .+ metricFunction2(x,y) .+ metricFunction3(x,y) .+ metricFunction4(x,y)
    elseif problem == 6
        M_func = (x,y) -> scale * GridGeneration.find_nearest_kd(metricData, tree, refs, x, y)
        return M_func(x,y)
    else
        error("Unknown problem type: $problem")
    end
end

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
end




function central_diff_nonuniform(coords,
                                 fvals; dir=:x)
    N = length(coords)
    @assert length(fvals) == N
    df = similar(fvals, Float64)

    # left boundary (forward difference)
    df[1] = (fvals[2] - fvals[1]) / (coords[2] - coords[1])

    # interior points (3-point nonuniform central)
    for i in 2:N-1
        xm, xi, xp = coords[i-1], coords[i], coords[i+1]
        fm, fi, fp = fvals[i-1], fvals[i], fvals[i+1]
        h0 = xi - xm
        h1 = xp - xi
        df[i] = ( -h1^2 * fm + (h1^2 - h0^2) * fi + h0^2 * fp ) / (h0*h1*(h0+h1))
    end

    # right boundary (backward difference)
    df[N] = (fvals[N] - fvals[N-1]) / (coords[N] - coords[N-1])

    return df
end


# set up real metric data
metricPath = "examples/airfoil/metric/A-airfoil_grid_data.mat"
# load metric
metricData = matread(metricPath)
# set up metric tree for fast nearest neighbor search
tree, refs = GridGeneration.setup_metric_tree(metricData)


# metricFunc = (x,y) -> Metric(x,y; problem = 1, scale = 100)

# metric_fg = make_getMetric_with_gradient(airfoil;
#     A_airfoil=100.0, ℓ_airfoil=0.5, p_airfoil=2,
#     A_origin=500.0,  ℓ_origin=0.1, p_origin=10,
#     origin_center=(0.0, 0.0),
#     floor=1e-4, profile=:rational,
# )

# metricFunc  = (x,y) -> metric_fg(x,y)[1]  # => [M, M]

# # Gradient duplicated per component:
# gradM = function (x,y)
#     _, (gx, gy) = metric_fg(x,y)
#     return ([gx, gx], [gy, gy])     # gradM[1] = [dMdx, dMdx], gradM[2] = [dMdy, dMdy]
# end

# gradMx = (x,y) -> gradM(x,y)[1]  # => [dMdx, dMdx]
# gradMy = (x,y) -> gradM(x,y)[2]  # => [dMdy, dMdy]

# gradMy = (x,y) -> metric_fg(x, y)[2][2]

# M11, M22 = metricFunc
# dMdx, dMdy = gradM

# # cast 2d points to 1d


# xs = GridGeneration.ProjectBoundary2Dto1D(boundarySection)

# m = GridGeneration.Get1DMetric(boundarySection, metricFunc)
# mx = GridGeneration.Get1DMetric(boundarySection, gradMx)


# function solver(


# get boundary
initialGrid = GetAirfoilGrid(airfoilPath = "examples/airfoil/data/A-airfoil.txt", radius = 3)
bottom = initialGrid[:,:,1]
airfoil = bottom[:, 101:end-100]  
boundarySection = airfoil




#[:, 100:500]#[:, sectionIndices]
# sectionIndices = 100:600

xs = GridGeneration.ProjectBoundary2Dto1D(boundarySection)


# define metric
metricFunc1 = make_getMetric(airfoil;
    A_airfoil = 50.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
    A_origin  = 5000.0,  ℓ_origin  = 0.1, p_origin  = 10,   
    floor     = 1e-4,  origin_center=(0, 0),
    profile   = :rational
)  

metricFunc2 = make_getMetric(airfoil;
    A_airfoil = 0.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
    A_origin  = 1000.0,  ℓ_origin  = 0.1, p_origin  = 10,   
    floor     = 1e-4,  origin_center=(0.5, 0.1),
    profile   = :rational
)  

metricFunc = (x,y) -> metricFunc1(x,y) .+  metricFunc2(x,y)


m = GridGeneration.Get1DMetric(boundarySection, metricFunc)
mx = central_diff_nonuniform(xs, m)
    
m_func = GridGeneration.build_interps_linear(xs, m)
mx_func = GridGeneration.build_interps_linear(xs, mx)

forcing = x -> mx_func(x) ./ (2 * m_func(x))
forcing_vals = forcing.(xs)

solS = SolveSecondOrder_MS(xs, forcing_vals; x0=0.0, x1=1.0)




mxAirfoil = central_diff_nonuniform(airfoil[1,:], getindex.(metricFunc.(airfoil[1,:], airfoil[2,:]),1))

m_func1 = GridGeneration.build_interps_linear(xs, m)
mx_func1 = GridGeneration.build_interps_linear(xs, mxAirfoil)

forcing1 = x -> mx_func1(x) ./ (2 * m_func1(x))

N = length(xs) - 2
omega = 0.5

x, solX, resNorm = GridGeneration.SolveSecondOrder(forcing1, xs[1], xs[end]; N=N, omega=omega, tol=1e-15, max_iter=1000, verbose=false)

# iters = 1:length(resNorm)
# pltRes = plot(iters,resNorm,title="L2 Norm of the Residual")
# display(pltRes)




p2 = plot(xs, solS, label="M_s", xlabel="s", ylabel="x(s)", title="Grid Generation along airfoil boundary", legend=:topleft)
plot!(p2, xs, solX, label="M_x", xlabel="s", ylabel="x(s)", title="Grid Generation along airfoil boundary", legend=:topleft)


# projectedSol = GridGeneration.ProjectBoundary1Dto2D(boundarySection, sol, xs)

# pltProjected = scatter(projectedSol[1, :], projectedSol[2, :], label="Projected Grid", xlabel="x", ylabel="y", title="Projected Grid onto Airfoil Boundary", aspect_ratio=1)







########### plots
# plot the metric and metric derivative values along the airfoil boundary



# mxvals = gradMx.(boundarySection[1,:], boundarySection[2,:])
# mxvalsx = getindex.(mxvals,1)

# N = length(mvalsx)
# xr = 1:N

# mplt = plot(xr, mvalsx, label="Mx", title="Mx along airfoil")
# # plot!(mplt, xr, mxvalsx, label="dMdx_x")
# display(mplt)

plt1DMetric = plot(xs, m, label="M(s)", xlabel="s", ylabel="M(s)", title="1D metric along airfoil boundary")
plot!(plt1DMetric, xs, mx, label="dM/ds", xlabel="s", ylabel="dM/ds", title="1D metric along airfoil boundary")


w = (x,y) -> first(metricFunc(x,y))
xl,xr = -0.5, 1.5
yl, yr = -1.0, 1.0
xs = range(xl, xr, length=450)
ys = range(yl, yr, length=300)

pltMetricField, _ = plot_scalar_field(w, xs, ys; boundary=airfoil,
                        title="Metric x field with airfoil",
                        cb_label="M(x,y)", equal_aspect=true, colormap =cgrad([RGB(1,1,1), RGB(1,0,0)]))#colormap=cgrd(:imola, rev=true))

plot!(pltMetricField, xlims=(xl, xr), ylims=(yl, yr))


display(pltMetricField)
display(plt1DMetric)
display(p2)
display(pltProjected)