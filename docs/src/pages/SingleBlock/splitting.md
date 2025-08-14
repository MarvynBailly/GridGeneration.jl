# Single Block with Splitting
Next, let's build `GridGeneration` functions that allow the user to add splits into the blocks. The flow I'm thinking will be
- Input single Tortuga block and desired split locations with the block
- break the block into sub-blocks at the split locations
  - Create new interfaces and add to `interInfo`
  - Update `bndInfo` for the blocks
  - Save block in blocks
- For (blockId, block) in enumerate(blocks)
  - for each dir in unvistedDirs[blockId]
    - get neighboring block ids
    - add this block and neighboring blocks to computeBlocks
    - for each edge in computeBlocks, solve for optimal $N$
    - take optimal $N$ for this direction to be the max of all opt $N$s
    - solve ODE with optimal $N$ along dir for each computeBlock
    - mark dir as visited for each computeBlock

Let's make this into two main `GridGeneration` functions: `GridGeneration.SplitBlock(block, splitLocations, bndInfo, interInfo)` and `GridGeneration.SolveAllBlocks(metric, blocks, bndInfo, interInfo).` And thus all we have to do is
- Input single Tortuga block and desired split locations with the block
- Split block into multi-block grid and update bndInfo and interInfo
- Solve for the optimal distribution in each block



## Split Block Algorithm
We want something along the lines:
```
    splitLocations = [
        [ horizontal split indices... ],   # split along the x axis
        [ vertical split indices... ]      # split along the y axis
    ]

    blocks, bndInfo, interInfo = GridGeneration.SplitBlock(block, splitLocations, bndInfo, interInfo)
```

## Solve All Blocks Algorithm

