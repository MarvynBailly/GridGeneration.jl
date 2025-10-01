function SmoothBlocks(blocks; solver=:ellipticSS, params)
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
    
    return smoothBlocks, finalErrors, finalIterations
end






######################## ADD A THREADED OPTION? ########################
    # @threads for i in eachindex(blocks)
    #     println("Solving for block $i on thread $(threadid())...")
    #     xr, yr, finalError, finalIter = GridGeneration.EllipticSolver(x = blocks[i][1, :, :], y = blocks[i][2, :, :],
    #                             max_iter = params.max_iter, tol = params.tol, ω = params.ω,
    #                             s_left = params.s_left, a_decay_left = params.a_decay_left, b_decay_left = params.b_decay_left,
    #                             s_right = params.s_right, a_decay_right = params.a_decay_right, b_decay_right = params.b_decay_right,
    #                             s_top = params.s_top, a_decay_top = params.a_decay_top, b_decay_top = params.b_decay_top,
    #                             s_bottom = params.s_bottom, a_decay_bottom = params.a_decay_bottom, b_decay_bottom = params.b_decay_bottom,
    #                             use_top_wall = params.use_top_wall, use_bottom_wall = params.use_bottom_wall, use_left_wall = params.use_left_wall, use_right_wall = params.use_right_wall,
    #                             verbose = params.EllipticSolver_verbose
    #                                 ) 

    #     smoothBlocks[i] = permutedims(cat(xr, yr, dims=3), (3,1,2))
    #     finalErrors[i] = finalError
    #     finalIterations[i] = finalIter
    # end