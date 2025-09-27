# In types.jl

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


# EllipticParams(; max_iter::Int=5000, tol::Float64=1e-6, ω::Float64=0.2,
#                  useTopWall::Bool=false, useBottomWall::Bool=true, useLeftWall::Bool=false, useRightWall::Bool=false,
#                  a_decay_left::Float64=0.4, b_decay_left::Float64=0.4,
#                  a_decay_right::Float64=0.4, b_decay_right::Float64=0.4,
#                  a_decay_top::Float64=0.4, b_decay_top::Float64=0.4,
#                  a_decay_bottom::Float64=0.4, b_decay_bottom::Float64=0.4,
#                  verbose::Bool=false) = EllipticParams(max_iter, tol, ω,
#                                                         useTopWall, useBottomWall, useLeftWall, useRightWall,
#                                                         a_decay_left, b_decay_left,
#                                                         a_decay_right, b_decay_right,
#                                                         a_decay_top, b_decay_top,
#                                                         a_decay_bottom, b_decay_bottom,
#                                                         verbose)

# SimParams(; useSplitting::Bool=true, 
#             splitLocations::Vector{Vector{Int}}=Vector{Vector{Int}}(), 
#             useEdgeSolver::Bool=true, boundarySolver::Symbol=:none,
#             useSmoothing::Bool=true, smoothMethod::Symbol=:none,
#             elliptic::EllipticParams=EllipticParams()) = SimParams(useSplitting, splitLocations, useEdgeSolver, boundarySolver, useSmoothing, smoothMethod, elliptic)


