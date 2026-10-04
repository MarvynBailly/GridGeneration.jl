"""
    GenerateGrid(initialGrid, bndInfo, interInfo, M; params=SimParams())
    GenerateGrid(blocks::Vector, bndInfo, interInfo, M; params=SimParams(), splitRequests=[])

Run the grid generation pipeline on a single-block `initialGrid` (a `[2, Ni, Nj]` array) or on a
vector of blocks:

1. split the block(s) (if `params.useSplitting`): a single block at `params.splitLocations` with
   `SplitBlock`; multiple blocks with [`SplitMultiBlock`](@ref) using `splitRequests`
   (`[(blockId, [[i_splits...], [j_splits...]]), ...]`), which propagates splits across interfaces,
2. redistribute the block edges according to the metric `M` and rebuild each block with TFI
   (if `params.useEdgeSolver`),
3. smooth each block with the elliptic solver (if `params.useSmoothing`).

Every stage is optional; a disabled stage passes its input through unchanged.

Returns a `NamedTuple` `(smoothBlocks, blocks, bndInfo, interInfo, finalErrors, finalIterations)`,
which can also be destructured positionally. `blocks` is the grid before smoothing; when smoothing
is disabled `smoothBlocks === blocks` and `finalErrors`/`finalIterations` are empty.
"""
function GenerateGrid(initialGrid::AbstractArray{<:Real,3}, bndInfo, interInfo, M; params::SimParams=SimParams())
    blocks = [initialGrid]

    if params.useSplitting
        @info "Splitting blocks..."
        blocks, bndInfo, interInfo = SplitBlock(initialGrid, params.splitLocations, bndInfo, interInfo)
    end

    return run_pipeline(blocks, bndInfo, interInfo, M, params)
end

function GenerateGrid(blocks::AbstractVector, bndInfo, interInfo, M; params::SimParams=SimParams(), splitRequests=[])
    all(b -> b isa AbstractArray{<:Real,3}, blocks) || throw(ArgumentError("blocks must be [2, Ni, Nj] arrays"))
    isempty(params.splitLocations) ||
        @warn "params.splitLocations is ignored for multi-block input; pass splitRequests instead"

    blocks = collect(blocks)
    if params.useSplitting && !isempty(splitRequests)
        @info "Splitting blocks..."
        blocks, bndInfo, interInfo = SplitMultiBlock(blocks, splitRequests, bndInfo, interInfo)
    end

    return run_pipeline(blocks, bndInfo, interInfo, M, params)
end

function run_pipeline(blocks, bndInfo, interInfo, M, params)
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
