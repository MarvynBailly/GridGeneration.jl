# Examples

## Setup (once)

From the repository root, create the examples environment and point it at your local checkout:

```bash
julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
```

## Scripts

| Script | What it does |
|---|---|
| `generalExample.jl` | Full pipeline on a built-in case: builds an initial grid (`case = :airfoil` or `:rectangle`), splits it, redistributes the block edges according to a metric, smooths the blocks, and saves a before/after plot to `examples/output/`. Takes about a minute for the airfoil case. |
| `generalExample_blank.jl` | Template for regenerating an existing Tortuga (`.grid` + `.metric`) grid. Fill in `metricFile` and `gridFile` at the top; the result is written to `examples/output/regenerated.grid`. |

Run either with:

```bash
julia --project=examples examples/generalExample.jl
```

## Case helpers

- `airfoil/`: C-grid setup around the A-airfoil (`GetAirfoilSetup`) and metric choices (`GetAirfoilMetric(problem)`, problems 1 to 6; problem 6 uses the bundled `A-airfoil_grid_data.mat`).
- `rectangle/`: rectangular domain (`GetRectangleDomain`) and hotspot metrics (`GetRectangleMetric(problem)`, problems 1 to 5).

A metric is any function `M(x, y) -> (M11, M22)` giving the diagonal metric tensor at a point;
`GridGeneration.make_getMetric` builds distance-based metrics around a polyline and/or a hotspot.

`archive/` holds old prototypes that no longer run against the current API.
