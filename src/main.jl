"""
    GenerateGrid(initialGrid, bndInfo, interInfo, M; params=SimParams())

Run the grid generation pipeline on a single-block `initialGrid` (a `[2, Ni, Nj]` array):

1. split the block at `params.splitLocations` (if `params.useSplitting`),
2. redistribute the block edges according to the metric `M` and rebuild each block with TFI
   (if `params.useEdgeSolver`),
3. smooth each block with the elliptic solver (if `params.useSmoothing`).

Every stage is optional; a disabled stage passes its input through unchanged.

Returns a `NamedTuple` `(smoothBlocks, blocks, bndInfo, interInfo, finalErrors, finalIterations)`,
which can also be destructured positionally. `blocks` is the grid before smoothing; when smoothing
is disabled `smoothBlocks === blocks` and `finalErrors`/`finalIterations` are empty.
"""
function GenerateGrid(initialGrid, bndInfo, interInfo, M; params::SimParams=SimParams())
    @assert ndims(initialGrid) == 3 "initialGrid must be a 3D array"

    blocks = [initialGrid]

    if params.useSplitting
        @info "Splitting blocks..."
        blocks, bndInfo, interInfo = SplitBlock(initialGrid, params.splitLocations, bndInfo, interInfo)
    end

    if params.useEdgeSolver
        @info "Solving boundary edges..."
        blocks, bndInfo, interInfo = SolveAllBlocks(M, deepcopy(blocks), deepcopy(bndInfo), deepcopy(interInfo); solver=params.boundarySolver)
    end

    smoothBlocks = blocks
    finalErrors = Float64[]
    finalIterations = Int[]

    if params.useSmoothing
        @info "Smoothing blocks..."
        smoothBlocks, finalErrors, finalIterations = SmoothBlocks(blocks; solver=params.smoothMethod, params=params.elliptic)
    end

    return (; smoothBlocks, blocks, bndInfo, interInfo, finalErrors, finalIterations)
end
