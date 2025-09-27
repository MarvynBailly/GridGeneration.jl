function GenerateGrid(initialGrid, bndInfo, interInfo, M; params=params)
    if params.useSplitting
        blocks, bndInfo, interInfo = GridGeneration.SplitBlock(initialGrid, params.splitLocations, bndInfo, interInfo)
    end
    if params.useEdgeSolver
        blocks, bndInfo, interInfo = GridGeneration.SolveAllBlocks(M, blocks, bndInfo, interInfo; solver=params.boundarySolver)
    end

    if params.useSmoothing
        smoothBlocks, finalErrors, finalIterations = GridGeneration.SmoothBlocks(blocks; solver=params.smoothMethod, params=params.elliptic)
    end

    return smoothBlocks, blocks, bndInfo, interInfo, finalErrors, finalIterations
end