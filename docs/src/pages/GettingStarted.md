# Getting Started

## Installation

GridGeneration.jl is registered in the General registry:

```julia
using Pkg
Pkg.add("GridGeneration")
```

To work on the package itself, clone the repository and run `Pkg.develop(path="path/to/GridGeneration.jl")`.

## Grid and metric conventions

- A **block** is a `[2, Ni, Nj]` array: `block[1, i, j]` and `block[2, i, j]` are the x and y coordinates of node `(i, j)`.
- An initial block can be built from its four edges with [`TFI`](@ref). The edges are given as
  `[top, right, bottom, left]`, each an `N×2` array, all running left→right or bottom→top.
- A **metric** is any function `M(x, y) -> (M11, M22)` returning the diagonal of the metric tensor.
  Larger values ask for finer spacing: with a constant metric `m`, the target spacing is `1/sqrt(m)`.
  [`make_getMetric`](@ref) builds distance-based metrics around a polyline and/or a hotspot.
- Boundary conditions (`bndInfo`) and block-to-block interfaces (`interInfo`) are vectors of
  dictionaries; see [Grid Format](./GridFormat.md). For a single block with no special boundaries,
  empty vectors are fine.

## A first grid

Build a 4×2 rectangle, ask for fine spacing near the point `(1, 0)`, split it into four blocks,
redistribute the points, and smooth:

```@example quickstart
using GridGeneration

N = 41
top    = [range(0, 4, length=N) fill(2.0, N)]
right  = [fill(4.0, N) range(0, 2, length=N)]
bottom = [range(0, 4, length=N) zeros(N)]
left   = [zeros(N) range(0, 2, length=N)]
initialGrid = TFI([top, right, bottom, left])

# background spacing ≈ 0.1 everywhere, refined around (1, 0)
hotspot = make_getMetric(nothing; A_origin = 2000.0, ℓ_origin = 0.3, p_origin = 2,
                         origin_center = (1.0, 0.0), floor = 100.0)

params = SimParams(
    splitLocations = [[11], [11]],          # split at i = 11 and j = 11
    boundarySolver = :analytic,
    smoothMethod   = :ellipticSS,
    elliptic       = EllipticParams(max_iter = 2000, useBottomWall = true),
)

result = GenerateGrid(initialGrid, Any[], Any[], hotspot; params = params)
size.(result.smoothBlocks)
```

`GenerateGrid` returns a `NamedTuple`:

```@example quickstart
keys(result)
```

- `blocks` is the grid after edge redistribution, `smoothBlocks` after elliptic smoothing.
- `bndInfo` and `interInfo` are updated to the new block sizes.
- `finalErrors` and `finalIterations` report the smoother's convergence for each block.

```@example quickstart
result.finalIterations
```

Each stage can be switched off through [`SimParams`](@ref) (`useSplitting`, `useEdgeSolver`,
`useSmoothing`), or run directly with [`SplitMultiBlock`](@ref), [`SolveAllBlocks`](@ref) and
[`SmoothBlocks`](@ref).

## Saving to a Tortuga grid file

Tortuga grids are 3D, so the 2D blocks are extruded first:

```julia
mesh3D, bnd3D, itf3D = convert_2D_to_3D(result.smoothBlocks, result.bndInfo, result.interInfo,
                                        0.1, 20)   # extrusion length, number of k-layers
write_turtle_grid(mesh3D, itf3D, bnd3D, "mygrid.grid")
```

An existing grid is read back with [`ImportTurtleGrid`](@ref), and a metric field file with
[`readTurtleFields`](@ref).

## More examples

The repository's `examples/` folder has a complete airfoil case and a template for regenerating
an existing Tortuga grid; see `examples/README.md`.
