function GenerateGrid(initialGrid, bndInfo, interInfo, M; params=params)
    @assert ndims(initialGrid) == 3 "initialGrid must be a 3D array"

    if params.useSplitting
        @info "Splitting blocks..."

        blocks, bndInfo, interInfo = GridGeneration.SplitBlock(initialGrid, params.splitLocations, bndInfo, interInfo)

        if params.useEdgeSolver === false && params.useSmoothing === false
            return blocks, bndInfo, interInfo
        end
    end

    if params.useEdgeSolver
        @info "Solving boundary edges..."

        if params.useSplitting === false
            blocks = [initialGrid]
        end


        blocksNew, bndInfoNew, interInfoNew = GridGeneration.SolveAllBlocks(M, deepcopy(blocks), deepcopy(bndInfo), deepcopy(interInfo); solver=params.boundarySolver)
        return blocksNew, bndInfoNew, interInfoNew
    end

    if params.useSmoothing
        @info "Smoothing blocks..."
        smoothBlocks, finalErrors, finalIterations = GridGeneration.SmoothBlocks(blocksNew; solver=params.smoothMethod, params=params.elliptic)
        return smoothBlocks, blocksNew, bndInfoNew, interInfoNew, finalErrors, finalIterations
    end
end