# GridGeneration.jl - AI Coding Agent Instructions

## Project Overview
GridGeneration.jl is a Julia package for generating structured computational grids using metric-based adaptive refinement. It solves 2D grid generation problems by reformulating them as 1D ODEs along grid lines, supporting both single and multi-block domains with elliptic smoothing.

**Core workflow**: Domain boundary → Initial TFI grid → Optional block splitting → Metric-adaptive edge solving → Optional elliptic smoothing

## Architecture & Key Components

### Three-Tier Grid Generation Pipeline
1. **Initialization** (`src/interpolators/TFI.jl`): Transfinite interpolation creates initial grids from boundary curves
2. **Block Processing** (`src/blocksplitting/`): Optional splitting into multi-block domains with automatic interface propagation
3. **Refinement** (`src/numerics/`, `src/smoothing/`): Metric-adaptive point distribution + elliptic smoothing

### Grid Data Format
- Grids are **3D arrays**: `[2, Ni, Nj]` where dim 1 = [x, y], dim 2 = ξ-direction, dim 3 = η-direction
- Boundaries follow **[top, right, bottom, left]** ordering as N×2 arrays
- Block splitting uses **index-based locations** on the initial grid: `[[i_splits...], [j_splits...]]`

### Metric Tensor System
Metric fields `M(x,y)` define desired point spacing. Two approaches:
- **Custom functions** (`src/metric/CustomMetric.jl`): `make_getMetric()` for analytical distance-based metrics (e.g., airfoil clustering)
- **Discrete data** (`src/metric/Metric.jl`): KDTree-based nearest-neighbor lookup for pre-computed metric fields

The metric field is queried during ODE solving to adapt grid point distribution.

### Solver Hierarchy
- **Boundary solvers** (`:analytic` or `:numeric`): Solve 1D ODEs along grid edges using semi-analytical integration or Thomas algorithm
- **Elliptic smoothers** (`:ellipticSS`): SOR iteration with optional wall forcing for orthogonality/spacing control
- Solver selection via `SimParams.boundarySolver` and `SimParams.smoothMethod`

## Critical Workflows

### Running Tests
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
Tests are modular (`test_solvers.jl`, `test_integration.jl`, `test_multiblock_splitting.jl`, `test_tortuga.jl`, ...),
orchestrated by `test/runtests.jl`, and include Aqua.jl quality checks.

### Building Documentation
```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'   # once
julia --project=docs docs/make.jl
```
Docs use Documenter.jl with sources in `docs/src/pages/`; output goes to `docs/build/` (ignored).
`docs/src/pages/Theory/` is generated from `math_notes/` by `make.jl` (edit the notes, not the copies).
`checkdocs = :exports`, so every exported symbol needs a docstring listed in `docs/src/pages/api.md`.

### Example Execution Pattern
See `examples/generalExample.jl` for canonical workflow:
1. Define domain boundaries (use helper functions in `examples/airfoil/data/` or `examples/rectangle/data/`)
2. Create initial grid via `TFI(boundary)`
3. Define metric field `M(x, y) -> (M11, M22)` (custom or from file)
4. Configure `SimParams` and `EllipticParams`
5. Call `GenerateGrid(initialGrid, bndInfo, interInfo, M; params=params)`

Returns a NamedTuple `(smoothBlocks, blocks, bndInfo, interInfo, finalErrors, finalIterations)`.
Run examples with `julia --project=examples examples/generalExample.jl` (see `examples/README.md`).

### Interactive GUI
A GLMakie desktop GUI lives in `gui/` (`gui/GridGenerationGUI.jl`, with its own `gui/Project.toml`).
Launch with `julia --project=gui gui/GridGenerationGUI.jl`. It loads a Tortuga grid and metric field
from `gui/` (these large data files are not tracked in git) and saves grids to `gui/output/`.
Provides interactive controls for split selection, solver configuration, and visualization.

## Project-Specific Patterns

### Multi-Block Splitting with Propagation
**Critical**: Splits automatically propagate across interfaces to maintain one-to-one grid correspondence
- Horizontal splits (j-direction) propagate through **vertical** interfaces (side-by-side blocks)
- Vertical splits (i-direction) propagate through **horizontal** interfaces (stacked blocks)
- Implementation: `src/blocksplitting/MultiBlockSplitting.jl::MapSplitAcrossInterface()`
- Tests verify propagation logic in `test/test_multiblock_splitting.jl` with detailed 2×2 grid scenarios

### Interface/Boundary Info Dictionaries
**bndInfo**: Array of boundary conditions with `"name"` and `"faces"` (each face has `"start"` and `"end"` indices)
**interInfo**: Array of interfaces with `"blockA"`, `"blockB"`, `"start_blkA"`, `"end_blkA"`, etc.

Example from `examples/airfoil/data/GetBoundary.jl` shows structure.

### Solver Selection Logic
- **Analytic solver** (`src/numerics/AnalyticSolver.jl`): Semi-analytic integration via cumulative trapezoidal rule, suitable for smooth metrics
- **Numeric solver** (`src/numerics/SecondOrderSolver.jl`): fixed-point iteration with a Thomas tridiagonal solve and under-relaxation (ω = 0.5)
- Choice affects `SolveAllBlocks()` and `GetOptNEdgePair()` behavior

### Elliptic Forcing Parameters
`EllipticParams` wall decay parameters (`a_decay_*`, `b_decay_*`) control boundary orthogonality enforcement:
- Typical values: 0.4-0.5 for moderate forcing
- Set `use*Wall=true` only for walls requiring orthogonality
- Verbose mode shows convergence per block

## Module Exports & API
**Main entry point**: `GenerateGrid()`
**Public types**: `SimParams`, `EllipticParams`
**Pipeline stages**: `SplitMultiBlock()`, `SolveAllBlocks()`, `SmoothBlocks()`
**Utilities**: `TFI()`, `make_getMetric()`, `setup_metric_tree()`, `find_nearest_kd()`
**Tortuga I/O**: `ImportTurtleGrid()`, `readTurtleFields()`, `convert_2D_to_3D()`, `write_turtle_grid()`

Private functions use `GridGeneration.` prefix. See `src/GridGeneration.jl` for complete export list.

## Development Notes
- **CI**: `.github/workflows/CI.yml` runs `Pkg.test()` on Julia 1.10 and 1; `Docs.yml` builds/deploys docs; CompatHelper and TagBot are configured
- **Documenter assets**: GIFs in `docs/src/assets/gifs/`, images in `docs/src/assets/images/` (prefer PNG over SVG for dense meshes)
- **Julia version**: `julia = "1.10"` compat in `Project.toml`
- **Dependencies**: Minimal (LinearAlgebra, NearestNeighbors); tests also use Aqua, Statistics, Test

### Common Pitfalls
1. Grid array indexing: Remember `[coord, i, j]` not `[i, j, coord]`
2. Boundary order: Always [top, right, bottom, left], and every edge runs left→right / bottom→top, or TFI will produce a twisted grid
3. Split locations: Indices refer to the **initial grid** before any splitting
4. Interface orientation: Vertical interfaces have varying i-index, horizontal have varying j-index

## File Organization Logic
- `src/`: Core module with subdirectories by functionality (numerics, interpolators, blocksplitting, metric, projections, smoothing)
- `examples/`: Complete use cases with domain setup and metric definition
- `test/`: One test file per major component, orchestrated by `runtests.jl`
- `gui/`: GLMakie GUI (separate environment)
- `math_notes/`: Theory notes, published into the docs' Theory section
- `docs/`: Documenter.jl setup with markdown sources in `src/pages/`

When adding features, follow existing subdirectory patterns and update relevant test files.
