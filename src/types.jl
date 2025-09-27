struct EllipticParams
    max_iter::Int
    tol::Float64
    ω::Float64
    useTopWall::Bool
    useBottomWall::Bool
    useLeftWall::Bool
    useRightWall::Bool
    a_decay_left::Float64
    b_decay_left::Float64
    a_decay_right::Float64
    b_decay_right::Float64
    a_decay_top::Float64
    b_decay_top::Float64
    a_decay_bottom::Float64
    b_decay_bottom::Float64
    verbose::Bool
end


struct SimParams
    useSplitting::Bool
    splitLocations::Vector{Vector{Int}}
    useEdgeSolver::Bool
    boundarySolver::Symbol
    useSmoothing::Bool
    smoothMethod::Symbol
    elliptic::EllipticParams
end


```
GridGeneration.EllipticParams
- max_iter::Int # maximum number of iterations
- tol::Float64 # tolerance for convergence
- ω::Float64 # relaxation factor
- useTopWall::Bool # whether to apply forcing on the top wall
- useBottomWall::Bool # whether to apply forcing on the bottom wall
- useLeftWall::Bool # whether to apply forcing on the left wall
- useRightWall::Bool # whether to apply forcing on the right wall
- a_decay_left::Float64 # decay parameter for left wall forcing
- b_decay_left::Float64 # decay parameter for left wall forcing
- a_decay_right::Float64 # decay parameter for right wall forcing
- b_decay_right::Float64 # decay parameter for right wall forcing
- a_decay_top::Float64 # decay parameter for top wall forcing
- b_decay_top::Float64 # decay parameter for top wall forcing
- a_decay_bottom::Float64 # decay parameter for bottom wall forcing
- b_decay_bottom::Float64 # decay parameter for bottom wall forcing
- verbose::Bool # whether to print convergence information
```
function EllipticParams(; max_iter::Int=5000, tol::Float64=1e-6, ω::Float64=0.2,
                 useTopWall::Bool=false, useBottomWall::Bool=true, useLeftWall::Bool=false, useRightWall::Bool=false,
                 a_decay_left::Float64=0.4, b_decay_left::Float64=0.4,
                 a_decay_right::Float64=0.4, b_decay_right::Float64=0.4,
                 a_decay_top::Float64=0.4, b_decay_top::Float64=0.4,
                 a_decay_bottom::Float64=0.4, b_decay_bottom::Float64=0.4,
                 verbose::Bool=false) 
                 return EllipticParams(
                    max_iter, tol, ω,
                    useTopWall, useBottomWall, useLeftWall, useRightWall,
                    a_decay_left, b_decay_left,
                    a_decay_right, b_decay_right,
                    a_decay_top, b_decay_top,
                    a_decay_bottom, b_decay_bottom,
                    verbose)
end

```
GridGeneration.SimParams
- useSplitting::Bool # whether to use block splitting
- splitLocations::Vector{Vector{Int}} # locations to split the grid, defined using indices of the initial grid
- useEdgeSolver::Bool # whether to use the edge solver on each block after splitting
- boundarySolver::Symbol # :analytic or :numeric, which boundary solver to use
- useSmoothing::Bool # whether to use smoothing on the final grid
- smoothMethod::Symbol # :ellipticSS, which smoothing method to use
- elliptic::EllipticParams # parameters for the elliptic solver if using elliptic smoothing
```
function SimParams(; useSplitting::Bool=true, 
            splitLocations::Vector{Vector{Int}}=Vector{Vector{Int}}(), 
            useEdgeSolver::Bool=true, boundarySolver::Symbol=:none,
            useSmoothing::Bool=true, smoothMethod::Symbol=:none,
            elliptic::EllipticParams=EllipticParams())
            return SimParams(useSplitting, splitLocations, useEdgeSolver, boundarySolver, useSmoothing, smoothMethod, elliptic)
end