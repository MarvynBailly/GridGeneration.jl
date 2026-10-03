# Multi-Block with Splitting

A multi-block input grid (for example one read from a Tortuga file with [`ImportTurtleGrid`](@ref))
is split with [`SplitMultiBlock`](@ref):

```julia
newBlocks, newBndInfo, newInterInfo = SplitMultiBlock(blocks, splitRequests, bndInfo, interInfo)
```

- `blocks` is a vector of `[2, Ni, Nj]` arrays.
- `splitRequests` lists, for each block to split, the `i` indices and the `j` indices to split at.
  It is either a vector of `(blockId, [[i_splits...], [j_splits...]])` tuples or a
  `Dict(blockId => [[i_splits...], [j_splits...]])`.
- `bndInfo` and `interInfo` use the package's 1-based indexing.

## Split propagation

Point-to-point matching across block interfaces must survive the split. A split line that ends on
an interface is therefore propagated into the neighbouring block at the matching index, and from
there onward until no new splits appear. For example, in a 2×2 layout, a `j`-split of the
lower-left block is carried across the vertical interface into the lower-right block. The upper
blocks are left alone.

```
Before                 After splitting block 1 at j = s
┌─────┬─────┐          ┌─────┬─────┐
│  3  │  4  │          │  5  │  6  │
├─────┼─────┤          ├─────┼─────┤
│  1  │  2  │          │  2  │  4  │
└─────┴─────┘          ├─────┼─────┤
                       │  1  │  3  │
                       └─────┴─────┘
```

Blocks are renumbered globally after splitting. Boundary faces are reassigned to the sub-blocks
that touch them, and interfaces are created between the new sub-blocks and remapped onto the
neighbours.

## Single-block input

For a single block, [`GenerateGrid`](@ref) uses `SplitBlock(block, splitLocations, bndInfo, interInfo)`,
where `splitLocations = [[i_splits...], [j_splits...]]`; see
[Single Block with Splitting](../SingleBlock/splitting.md).
