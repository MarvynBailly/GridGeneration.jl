# GridGeneration.jl - AI Coding Agent Instructions

## Project Overview
GridGeneration.jl is a Julia package for generating structured computational grids using metric-based mesh adaptation. The core approach solves a second-order ODE derived from a metric tensor field to achieve optimal point spacing, with support for multi-block grids and elliptic smoothing.

## Architecture

### Core Pipeline (see `src/main.jl::GenerateGrid`)
1. **Block Splitting** (`blocksplitting/`) - Subdivide initial grid at specified indices
2. **Edge Solving** (`numerics/`) - Solve ODE on block boundaries for optimal spacing  
3. **Interior Interpolation** (`interpolators/TFI.jl`) - Use transfinite interpolation to fill block interiors
4. **Smoothing** (`smoothing/`) - Apply elliptic PDE solver for final grid quality

### Data Flow: 2D Grid → 1D Projection → ODE Solution → 2D Reconstruction
- **Projection** (`projections/1D2DFunctions.jl`): Maps 2D boundary curves to 1D parametric space
- **Metric Computation** (`projections/Get1DMetric.jl`): Evaluates metric tensor along projected curves
- **ODE Solving** (`numerics/NumSolver.jl`): Computes optimal point distribution via `SolveODE()` or `SolveODEFixedN()`
- **Reconstruction** (`projections/1D2DFunctions.jl`): Maps 1D solution back to 2D physical space

## Key Conventions

### Grid Representation
```julia
# Grids are 3D arrays: [coordinate, i-direction, j-direction]
block = zeros(2, Ni, Nj)  # block[1,:,:] = x-coords, block[2,:,:] = y-coords
```

### Boundary Order Convention
Boundaries are **always** specified as `[top, right, bottom, left]`:
```julia
boundary = [topEdge, rightEdge, bottomEdge, leftEdge]  # Each is Nx2 array
initialGrid = GridGeneration.TFI(boundary)
```

### Solver Selection
- `:analytic` - Semi-analytical method (default, stable, uses numerical integration)
- `:numeric` - Iterative central differencing (unstable, avoid unless necessary)

### Block Metadata
- **bndInfo**: Boundary condition dictionaries grouped by name
- **interInfo**: Interface connectivity between blocks (stored as `Dict` with blockA/blockB keys)
- Both updated automatically by `UpdateBndInfo!()` and `UpdateInterInfo()` after grid changes

## Critical Workflows

### Running Examples
```julia
# From project root
using Pkg; Pkg.activate(".")
include("examples/generalExample.jl")
```
Examples use `include()` pattern, not module imports. See `examples/airfoil/` and `examples/rectangle/`.

### Building Documentation
```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```
Docs auto-deploy via GitHub Actions on push to `main` or version tags.

### Parameter Configuration
Use `SimParams` and `EllipticParams` structs (see `src/types.jl`):
```julia
params = GridGeneration.SimParams(
    useSplitting = true,
    splitLocations = [[300, 400], [30]],  # [x-splits, y-splits] as grid indices
    boundarySolver = :analytic,
    useSmoothing = true,
    smoothMethod = :ellipticSS,
    elliptic = GridGeneration.EllipticParams(...)
)
```

## Module Structure Peculiarities

### No Explicit Exports
Functions are accessed via `GridGeneration.FunctionName()` - the module does **not** export symbols. Always use fully qualified names.

### Include-Based Organization
`src/GridGeneration.jl` uses `include()` to compose submodules. Dependencies must be manually managed between files.

### Metric Tensor Interface
Metric functions must return `[M11, M22]` (diagonal components) for a given `(x, y)`:
```julia
M = GetAirfoilMetric(problem; scale=0.05)
metricValues = M(x, y)  # Returns [M11, M22]
```
See `examples/airfoil/metric/GetMetric.jl` for loading from `.mat` files using k-d tree lookup.

## Common Patterns

### Multi-Block Neighbor Resolution
`GetNeighbors(blockId, interInfo, dir)` returns all blocks sharing an edge in direction `dir` (1=horizontal, 2=vertical). Used to synchronize optimal point counts across interfaces.

### Optimal Point Calculation
```julia
optN = GetOptNEdgePair(leftEdge, rightEdge, metricFunc)  # Computes max of both edges
block = SolveBlockFixedN(block, bndInfo, interInfo, metricFunc, [Ni, Nj])
```

### Logging
Use `@info` for pipeline stages (splitting, solving, smoothing). Most debug `println()` calls are commented out - search for them when debugging.

## Testing
**No formal test suite exists.** Validation is done via examples. When modifying core numerics, run `examples/generalExample.jl` with both `:airfoil` and `:rectangle` cases.

## Version & Dependencies
- Julia 1.10.4 (pinned)
- Key deps: `NearestNeighbors.jl` (metric lookup), `Documenter.jl`, `LinearAlgebra`
- Versioning follows semver (currently v2.2.0)
