# GridGeneration.jl

[![Docs (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://MarvynBailly.github.io/GridGeneration.jl/dev/)
[![CI](https://github.com/MarvynBailly/GridGeneration.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/MarvynBailly/GridGeneration.jl/actions/workflows/CI.yml)

Metric-driven structured multi-block grid generation in Julia. Starting from an initial grid,
GridGeneration.jl splits it into blocks, redistributes each block's edge points so the spacing
follows a user-supplied metric field (choosing the number of points automatically), rebuilds the
blocks by transfinite interpolation, and smooths them with an elliptic solver. Grids can be read
from and written to the Tortuga (`.grid`) format.

## Installation

```julia
using Pkg
Pkg.add("GridGeneration")
```

## Quick start

```julia
using GridGeneration

N = 41
top    = [range(0, 4, length=N) fill(2.0, N)]
right  = [fill(4.0, N) range(0, 2, length=N)]
bottom = [range(0, 4, length=N) zeros(N)]
left   = [zeros(N) range(0, 2, length=N)]
initialGrid = TFI([top, right, bottom, left])           # [2, Ni, Nj] array

# metric M(x, y) -> (M11, M22): spacing ≈ 0.1, refined near (1, 0)
M = make_getMetric(nothing; A_origin = 2000.0, ℓ_origin = 0.3, p_origin = 2,
                   origin_center = (1.0, 0.0), floor = 100.0)

params = SimParams(splitLocations = [[11], [11]])
result = GenerateGrid(initialGrid, Any[], Any[], M; params = params)
result.smoothBlocks   # the final blocks
```

See [Getting Started](https://MarvynBailly.github.io/GridGeneration.jl/dev/pages/GettingStarted/)
for the conventions and the full pipeline.

## Examples and GUI

- [`examples/`](examples/README.md): a complete airfoil case and a template for regenerating an existing Tortuga grid.
- [`gui/`](gui/README.md): an interactive GLMakie GUI.

## Development

```bash
julia --project=. -e 'using Pkg; Pkg.test()'      # tests (includes Aqua checks)
julia --project=docs docs/make.jl                 # docs (after Pkg.develop(path=".") in docs/)
```

The theory behind the elliptic smoother lives in [`math_notes/`](math_notes/) and is published in the docs.
