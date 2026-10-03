# Airfoil Examples

These figures come from `examples/generalExample.jl` with `case = :airfoil`. The example builds a
C-grid around the A-airfoil, splits it into blocks, redistributes the block edges according to a
metric, and smooths the result. Each panel shows the input grid, the blocks after splitting, the
metric, the blocks after the edge solve, and the final grid.

The metric choices are `GetAirfoilMetric(problem)` in `examples/airfoil/metric/GetMetric.jl`:

```julia
initialGrid, bndInfo, interInfo = GetAirfoilSetup(radius = 3, type = :cgrid)
M = GetAirfoilMetric(problem; scale = 4000)   # problem = 1, ..., 6
params = SimParams(splitLocations = [[300, 400], [30]], boundarySolver = :analytic,
                   smoothMethod = :ellipticSS)
result = GenerateGrid(initialGrid, bndInfo, interInfo, M; params = params)
```

See `examples/README.md` for how to set up and run the example.

## Uniform (problem 1)
![uniform-metric](../../assets/images/Examples/airfoil/airfoil_all_1.png)

## Leading Edge (problem 2)
![leading-metric](../../assets/images/Examples/airfoil/airfoil_all_2.png)

## Trailing Edge (problem 3)
![trailing-metric](../../assets/images/Examples/airfoil/airfoil_all_3.png)

## Leading and Trailing Edge (problem 4)
![leadingandtrailing-metric](../../assets/images/Examples/airfoil/airfoil_all_4.png)

## Custom Metric (problem 5)
![custom-metric](../../assets/images/Examples/airfoil/airfoil_all_5.png)
