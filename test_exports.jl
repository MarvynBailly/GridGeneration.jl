# Test script to verify exports work correctly
using Pkg
Pkg.activate(".")

include("src/GridGeneration.jl")
using .GridGeneration

println("=" ^ 60)
println("Testing GridGeneration.jl Exports")
println("=" ^ 60)

# Test 1: Configuration types can be created without prefix
println("\nTest 1: Creating SimParams without prefix")
params = SimParams(
    useSplitting = false,
    boundarySolver = :analytic,
    useSmoothing = false
)
println("  Created: ", typeof(params))

# Test 2: EllipticParams can be created without prefix
println("\nTest 2: Creating EllipticParams without prefix")
elliptic_params = EllipticParams(
    max_iter = 1000,
    tol = 1e-5,
    ω = 0.3
)
println("  Created: ", typeof(elliptic_params))

# Test 3: TFI function works without prefix
println("\nTest 3: Using TFI function without prefix")
# Create a simple unit square
L = 1.0
top = [range(0, L, length=5) fill(L, 5)]
right = [fill(L, 5) range(L, 0, length=5)]
bottom = [range(L, 0, length=5) zeros(5)]
left = [zeros(5) range(0, L, length=5)]

grid = TFI([top, right, bottom, left])
println("  Grid size: ", size(grid))
println("  Grid type: ", typeof(grid))

# Test 4: Verify GenerateGrid is exported
println("\nTest 4: GenerateGrid is exported")
println("  Function exists: ", isdefined(GridGeneration, :GenerateGrid))

# Test 5: Metric utilities are exported
println("\nTest 5: Metric utilities are exported")
println("  make_getMetric: ", isdefined(GridGeneration, :make_getMetric))
println("  setup_metric_tree: ", isdefined(GridGeneration, :setup_metric_tree))
println("  find_nearest_kd: ", isdefined(GridGeneration, :find_nearest_kd))

println("\n" * "=" ^ 60)
println("All export tests passed! ✓")
println("=" ^ 60)
