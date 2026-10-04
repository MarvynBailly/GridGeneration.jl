# GridGeneration.jl release notes

## v4.0.0

### Breaking changes
- **Smoothed grids change.** The left/right (ξ-wall) forcing in the elliptic smoother had its
  α and γ coefficients swapped; this is fixed, so grids smoothed with left/right wall forcing
  differ from earlier versions. Decay parameters tuned against the old behaviour may need adjusting.
- **`:numeric` boundary solver results change.** It previously stopped after 10 heavily
  under-relaxed iterations, well short of convergence; it now converges and agrees with `:analytic`.
- **`SimParams` defaults** are now `boundarySolver = :analytic` and `smoothMethod = :ellipticSS`
  (previously the invalid `:none`); unknown solver names throw an `ArgumentError`.
- **`EllipticParams`** has a new `sweep` field, so its positional constructor changed (the keyword
  constructor is unchanged).
- `ComputeOptimalNumberofPoints` integrates over every interval and never returns fewer than 3
  points, so optimal point counts can differ slightly.

### New features
- `GenerateGrid` accepts a vector of blocks with `splitRequests` (multi-block input via
  `SplitMultiBlock`) and returns a `NamedTuple`
  `(smoothBlocks, blocks, bndInfo, interInfo, finalErrors, finalIterations)`, which still
  destructures positionally.
- `load_turtle_grid` and `setup_turtle_grid_domain` load Tortuga grids and metric fields in the
  package's conventions; the metric is looked up on the grid it was computed on (`metricGridFile`).
- `EllipticParams(sweep = :line)`: alternating line Gauss–Seidel option for the elliptic smoother.
- `ComputeAngleDeviation` for grid orthogonality.
- Newly exported: `SolveAllBlocks`, `SmoothBlocks`, `ImportTurtleGrid`, `readTurtleFields`,
  `convert_2D_to_3D`, `write_turtle_grid`.

### Fixes
- `GenerateGrid` could not reach the smoothing stage and errored without an explicit `params`.
- `SmoothBlocks` failed with a single `EllipticParams`.
- `write_turtle_grid` failed on output of `convert_2D_to_3D` (`KeyError: "faceInfo"`).
- Uneven-interface alignment split blocks one node too far; interfaces whose second block was
  split were not remapped.
- `SplitMultiBlock` could map interfaces onto the wrong edges when the first block of an
  interface was the upper/right one, or when the two sides ran in opposite directions.
- The elliptic convergence check now includes y, errors on divergence, and warns at `max_iter`.
- Tortuga readers raise clear errors on unsupported files instead of returning `nothing`.

### Performance
- The elliptic smoother is 8–13× faster and allocates almost nothing per iteration.

### Other
- `Pkg.test()` works (`test/runtests.jl`); the test suite covers the solvers, splitting,
  Tortuga I/O and the full pipeline, and includes Aqua.jl checks.
- Documenter is no longer a package dependency.
- New documentation: Getting Started, API reference, and a Theory section on elliptic smoothing.
