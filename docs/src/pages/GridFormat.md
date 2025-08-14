# Grid Format

To save multi-block grids, we will use "Tortuga Format" used by Dr. Larsson's [Computational Turbulence Laboratory](https://larsson.umd.edu/research/). 

The format saves three pieces of information which can easily be dumped to binary and loaded into his in-house solver. All grids inputted into GridGeneration.jl are expected to be in this format.

## Boundary Information
An array of dicts of the form

```julia
Dict("faces" => faces, "name" => name)
```

where faces is an array of dicts

```julia
faces = Any[ 
    Dict("start" => [startni, startnj, startnk], "end" => [endni, endnj, endnk], "blockId" = blockId),
    Dict("start" => [startni, startnj, startnk], "end" => [endni, endnj, endnk], "blockId" = blockId),
]
```


## Interface Information
An array of dicts of the forms

```julia
Dict(
    "blockA" => blockAId, "start_blkA" => [niStartA,njStartA,nkStartA], "end_blkA" => [niEndA,njEndA,nkEndA],
    "blockB" => blockBId, "start_blkB" => [niStartB,njStartB,nkStartB], "end_blkB" => [niEndB,njEndB,nkEndB],
    "offset" => [0.0, 0.0, 0.0], "angle" => 0.0
)
```

## Block Information
An array of arrays 

```julia
blocks = [block1, block2, ...]
```
where each block is of size `(3, ni, nj, nk)` containing the $(i,j,k)$ coordinates for each grid point within this block.