using Plots
using LinearAlgebra


"""
Helper function to calculate the y-coordinates of a symmetric NACA 4-digit airfoil.
For NACA 0012, t=0.12.
"""
function naca_0012(x, c=1.0)
    # Thickness for NACA 0012 is 12% of the chord
    t = 0.12
    # Standard formula for NACA 4-digit airfoils (half-thickness)
    y_t = 5 * t * c * (0.2969 * sqrt(x/c) 
                     - 0.1260 * (x/c) 
                     - 0.3516 * (x/c)^2 
                     + 0.2843 * (x/c)^3 
                     - 0.1015 * (x/c)^4) # Note: -0.1036, -0.1015
    return y_t
end

function setup_airfoil(Ni::Int, Nj::Int)
    @assert isodd(Ni) "Ni must be an odd number to have a center point at the leading edge."

    x = zeros(Float64, Ni, Nj)
    y = zeros(Float64, Ni, Nj)

    # --- Parameters ---
    chord_length = 1.0
    far_field_radius = 5.0

    # --- 1. Define the boundaries ---
    
    # Bottom boundary (j=1): The airfoil surface
    # We use cosine clustering to get more points at the leading and trailing edges.
    Ni_half = (Ni + 1) ÷ 2
    x_dist = 0.5 * chord_length * (1.0 .- cos.(range(0, π, length=Ni_half)))
    
    # Upper surface (from leading edge to trailing edge)
    x[Ni_half:Ni, 1] .= x_dist
    y[Ni_half:Ni, 1] .= naca_0012.(x_dist, chord_length)
    
    # Lower surface (from trailing edge to leading edge)
    # Note: We reverse the arrays and skip the leading edge point which is already defined.
    x[1:Ni_half-1, 1] .= reverse(x_dist[1:end-1])
    y[1:Ni_half-1, 1] .= -naca_0012.(reverse(x_dist[1:end-1]), chord_length)

    # Top boundary (j=Nj): Far-field boundary (a semi-circle)
    # The semi-circle is centered around the airfoil's mid-chord
    theta = range(2*π, 0, length=Ni)
    x[:, Nj] .= chord_length / 2.0 .+ far_field_radius .* cos.(theta)
    y[:, Nj] .= far_field_radius .* sin.(theta)
    
    # Left (i=1) and Right (i=Ni) boundaries connect the trailing edge to the far-field.
    # We define them via linear interpolation.
    y[1, :] = range(y[1, 1], y[1, Nj], length=Nj)
    x[1, :] = range(x[1, 1], x[1, Nj], length=Nj)
    
    y[Ni, :] = range(y[Ni, 1], y[Ni, Nj], length=Nj)
    x[Ni, :] = x[1, :] # In a C-grid, the wake cut lines have the same x-coords

    # --- 2. Make an initial guess for the interior grid (Linear Interpolation) ---
    # This interpolates between the airfoil surface and the far-field boundary.
    for i in 2:Ni-1
        for j in 2:Nj-1
            # Linear interpolation factor
            eta_j = (j - 1) / (Nj - 1)
            
            # Interpolate between the airfoil (j=1) and far-field (j=Nj) boundaries
            x[i, j] = (1 - eta_j) * x[i, 1] + eta_j * x[i, Nj]
            y[i, j] = (1 - eta_j) * y[i, 1] + eta_j * y[i, Nj]
        end
    end
    
    return x, y
end


"""
Sets up the physical domain and creates an initial algebraic grid.
The domain is a channel with a sinusoidal bump on the bottom wall.
"""
function setup_geometry_and_initial_grid(Ni::Int, Nj::Int)
    x = zeros(Float64, Ni, Nj)
    y = zeros(Float64, Ni, Nj)


    # Physical domain extents
    x_min, x_max = 0.0, 10.0
    y_max = 5.0 # Channel height

    # 1. Define the boundaries
    xi_pts = range(x_min, x_max, length=Ni)
    
    # Bottom boundary (j=1): Flat plate with a sinusoidal bump
    for i in 1:Ni
        x_val = xi_pts[i]
        x[i, 1] = x_val
        # A smooth bump in the middle of the domain
        if 1.0 <= x_val <= 9.0
            y[i, 1] = 1.0 * (1 + cos(pi * (x_val - 2.0)))
        else
            y[i, 1] = 0.0
        end
    end

    # Top boundary (j=Nj)
    # x[:, Nj] = xi_pts
    for i in 1:Ni
        x_val = xi_pts[i]
        x[i, Nj] = x_val
        #a smooth bump
        if 1.0 <= x_val <= 9.0
            y[i, Nj] = y_max - 1.0 * (1 + cos(pi * (x_val - 2.0)))
        else
            y[i, Nj] = y_max
        end
    end

    # Left (i=1) and Right (i=Ni) boundaries are straight vertical lines
    y[1, :] = range(y[1, 1], y_max, length=Nj)
    x[1, :] .= x_min
    y[Ni, :] = range(y[Ni, 1], y_max, length=Nj)
    x[Ni, :] .= x_max

    # 2. Make an initial guess for the interior grid (Linear Interpolation)
    for i in 2:Ni-1
        for j in 2:Nj-1
            eta_j = (j - 1) / (Nj - 1)
            x[i, j] = (1 - eta_j) * x[i, 1] + eta_j * x[i, Nj]
            y[i, j] = (1 - eta_j) * y[i, 1] + eta_j * y[i, Nj]
        end
    end
    
    return x, y
end


function plot_grid(x, y, title_str; plt = nothing, c = nothing, lw = 0.1 )
    p = plt == nothing ? plot() : plt
    plot!(p, title=title_str, aspect_ratio=:equal, legend=false, framestyle=:box, grid = false)
    # Plot η-lines (lines of constant j)
    for j in 1:1:size(x, 2)
        plot!(p, x[:, j], y[:, j], lw=lw, color=c)
    end
    # Plot ξ-lines (lines of constant i)
    for i in 1:5:size(x, 1)
        plot!(p, x[i, :], y[i, :], lw=lw, color=c)
    end
    return p
end


struct EllipticParams
    Ni::Int
    Nj::Int
    max_iter::Int
    tol::Float64
    ω::Float64
    use_top_wall::Bool
    s_top::Float64
    a_decay_top::Float64
    b_decay_top::Float64
    use_left_wall::Bool
    s_left::Float64
    a_decay_left::Float64
    b_decay_left::Float64
    use_right_wall::Bool
    s_right::Float64
    a_decay_right::Float64
    b_decay_right::Float64
    use_bottom_wall::Bool
    s_bottom::Float64
    a_decay_bottom::Float64
    b_decay_bottom::Float64
    EllipticSolver_verbose::Bool
    showPlots::Bool
end


##############################################
############ Elliptic Parameters #############
##############################################
Ni = 101
Nj = 41

max_iter::Int=1000
tol::Float64=1e-6
ω::Float64= 0.5

# note user defined `s_side` is being ignored in favor of automatic calculation based on grid spacing

use_top_wall::Bool=true 
s_top::Float64=-1
a_decay_top::Float64=0.9
b_decay_top::Float64=0.9

use_left_wall::Bool=false 
s_left::Float64=-1
a_decay_left::Float64=0.8
b_decay_left::Float64=0.8


use_right_wall::Bool=false
s_right::Float64=-1
a_decay_right::Float64=0.8
b_decay_right::Float64=0.8

use_bottom_wall::Bool=true 
s_bottom::Float64=-1
a_decay_bottom::Float64=0.1
b_decay_bottom::Float64=0.1

EllipticSolver_verbose::Bool=false

showPlots::Bool=false

EllipticParams(; Ni = Ni, Nj = Nj,
    max_iter = max_iter,
    tol = tol,
    ω = ω,
    use_top_wall = use_top_wall,
    s_top = s_top,
    a_decay_top = a_decay_top,
    b_decay_top = b_decay_top,
    use_left_wall = use_left_wall,
    s_left = s_left,
    a_decay_left = a_decay_left,
    b_decay_left = b_decay_left,
    use_right_wall = use_right_wall,
    s_right = s_right,
    a_decay_right = a_decay_right,
    b_decay_right = b_decay_right,
    use_bottom_wall = use_bottom_wall,
    s_bottom = s_bottom,
    a_decay_bottom = a_decay_bottom,
    b_decay_bottom = b_decay_bottom,
    EllipticSolver_verbose = EllipticSolver_verbose,
    showPlots = showPlots
) = EllipticParams(Ni, Nj, max_iter, tol, ω, use_top_wall, s_top, a_decay_top, b_decay_top,
    use_left_wall, s_left, a_decay_left, b_decay_left,
    use_right_wall, s_right, a_decay_right, b_decay_right,
    use_bottom_wall, s_bottom, a_decay_bottom, b_decay_bottom,
    EllipticSolver_verbose, showPlots
)


    



function main(params)
    # x_init, y_init = setup_geometry_and_initial_grid(params.Ni, params.Nj)
    # x_init, y_init = setup_airfoil(params.Ni, params.Nj)

    airfoilOuterRadius::Float64 = 0.8
    initialGrid = GetAirfoilGrid(airfoilPath="airfoil/data/A-airfoil.txt", radius = airfoilOuterRadius)
    airfoilGrid = initialGrid[:, 101:end-100,:]
    x_init = airfoilGrid[1, :, :]
    y_init = airfoilGrid[2, :, :]

    xr, yr, finalError, finalIter = EllipticSolver(x = copy(x_init), y = copy(y_init),
                            max_iter = params.max_iter, tol = params.tol, ω = params.ω,
                            s_left = params.s_left, a_decay_left = params.a_decay_left, b_decay_left = params.b_decay_left,
                            s_right = params.s_right, a_decay_right = params.a_decay_right, b_decay_right = params.b_decay_right,
                            s_top = params.s_top, a_decay_top = params.a_decay_top, b_decay_top = params.b_decay_top,
                            s_bottom = params.s_bottom, a_decay_bottom = params.a_decay_bottom, b_decay_bottom = params.b_decay_bottom,
                            use_top_wall = params.use_top_wall, use_bottom_wall = params.use_bottom_wall, use_left_wall = params.use_left_wall, use_right_wall = params.use_right_wall,
                            verbose = params.EllipticSolver_verbose
                            )
    return xr, yr, finalError, finalIter, x_init, y_init
end


runSolver = true
saveFig = false
case = 4

if case == 1
    params = EllipticParams(
        use_left_wall = false, 
        use_right_wall = false, 
        use_bottom_wall = false, 
        use_top_wall = false)
elseif case == 2
    params = EllipticParams(
        use_left_wall = false, 
        use_right_wall = false, 
        use_bottom_wall = true,
        use_top_wall = false)
elseif case == 3
    params = EllipticParams(
        use_left_wall = false, 
        use_right_wall = false, 
        use_bottom_wall = true,
        use_top_wall = true)
elseif case == 4
    params = EllipticParams(
        use_left_wall = true, 
        use_right_wall = true, 
        use_bottom_wall = true,
        use_top_wall = true)
else
    error("Invalid case number. Please choose between 1 and 4.")
end

@info "runSolver = $runSolver..."

if runSolver
    x, y, finalError, finalIter, x_init, y_init = main(params)
end

@info "Final error: $finalError"
@info "Final iterations: $finalIter"


pInit = plot_grid(x_init, y_init, "Initial Grid", c = RGB(0.0, 0.0, 0.0), lw=0.2)
pFinal = plot_grid(x, y, "Final Grid After Elliptic Case $case", c = RGB(0.0, 0.0, 0.0), lw=0.2)
p = plot(pInit, pFinal, layout=(1,2), size=(900,400))
display(p)

if saveFig
    savefig(p, "images/elliptic_grid_$case.pdf")
end