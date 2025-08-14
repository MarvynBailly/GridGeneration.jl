# Multi-Block with Splitting

To allow the user to input a multi-block `Tortuga` grid, all that we have to do is add compatibility to the `GridGeneration.SplitBlock(block, splitLocations, bndInfo, interInfo)` to allow for multiple blocks rather than the single block. Let's call this function `GridGeneration.SplitBlocks(blocks, splitLocations, bndInfo, interInfo).`

## Split Blocks Algorithm
