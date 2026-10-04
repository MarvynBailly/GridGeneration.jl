"""
    SmoothBlocks(blocks; solver=:ellipticSS, params)

Smooth grid blocks using the elliptic PDE solver. `params` is either a single
`EllipticParams` applied to every block or a vector with one `EllipticParams` per block.

Returns `(smoothBlocks, finalErrors, finalIterations)`.
"""
function SmoothBlocks(blocks; solver=:ellipticSS, params)
    solver == :ellipticSS || throw(ArgumentError("Unknown smoothing method :$solver (expected :ellipticSS)"))

    if params isa EllipticParams
        params = fill(params, length(blocks))
    end
    length(params) == length(blocks) ||
        throw(ArgumentError("got $(length(params)) EllipticParams for $(length(blocks)) blocks"))

    smoothBlocks = Vector{Array{Float64,3}}(undef, length(blocks))
    finalErrors = Vector{Float64}(undef, length(blocks))
    finalIterations = Vector{Int}(undef, length(blocks))

    if solver == :ellipticSS

        for i in eachindex(blocks)
            if params[i].skipBlock
                @info "Skipping block $i as per parameters."
                smoothBlocks[i] = blocks[i]
                finalErrors[i] = 0.0
                finalIterations[i] = 0
                continue
            end
            xr, yr, finalError, finalIter = GridGeneration.EllipticSolver(blocks[i][1, :, :], blocks[i][2, :, :], params = params[i]) 

            smoothBlocks[i] = permutedims(cat(xr, yr, dims=3), (3,1,2))
            finalErrors[i] = finalError
            finalIterations[i] = finalIter
        end
    end
     
    @info "Smoothing convergence: $finalErrors"
    @info "Smoothing iterations: $finalIterations"
    
    return smoothBlocks, finalErrors, finalIterations
end