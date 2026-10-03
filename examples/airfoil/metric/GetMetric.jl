using DelimitedFiles
using MAT: matread

const AIRFOIL_DATA_PATH = joinpath(@__DIR__, "..", "data", "A-airfoil.txt")
const AIRFOIL_METRIC_PATH = joinpath(@__DIR__, "A-airfoil_grid_data.mat")

"""
    GetAirfoilMetric(problem = 1; scale = 4000, airfoilPath = AIRFOIL_DATA_PATH)

Return a metric function `M(x, y) -> (M11, M22)` for the airfoil example.

- problem 1: constant metric
- problem 2: leading edge
- problem 3: trailing edge
- problem 4: leading edge and trailing edge
- problem 5: for fun
- problem 6: real metric (nearest-neighbour lookup in `A-airfoil_grid_data.mat`)
"""
function GetAirfoilMetric(problem = 1; scale = 4000, airfoilPath = AIRFOIL_DATA_PATH)
    # airfoil surface as an N×2 polyline (the data file has x, y, z columns)
    airfoil = Float64.(readdlm(airfoilPath, '\t'; skipstart=1)[:, 1:2])

    hotspot(center; A = scale) = GridGeneration.make_getMetric(airfoil;
        A_airfoil = 100.0, ℓ_airfoil = 0.5, p_airfoil = 2,
        A_origin = A, ℓ_origin = 0.1, p_origin = 10,
        floor = 1e-4, origin_center = center,
        profile = :rational
    )
    sum_metrics(fs...) = (x, y) -> reduce((a, b) -> a .+ b, (f(x, y) for f in fs))

    if problem == 1
        return (x, y) -> (scale, scale)
    elseif problem == 2
        return hotspot((0, 0))
    elseif problem == 3
        return hotspot((1, 0))
    elseif problem == 4
        return sum_metrics(hotspot((0, 0)), hotspot((1, 0)))
    elseif problem == 5
        return sum_metrics(hotspot((0, 0); A = 1000.0), hotspot((1, 0); A = 1000.0),
                           hotspot((0.30, 0.13); A = 1000.0), hotspot((0.50, -0.09); A = 1000.0))
    elseif problem == 6
        metricData = matread(AIRFOIL_METRIC_PATH)
        tree, refs = GridGeneration.setup_metric_tree(metricData)
        return (x, y) -> scale * GridGeneration.find_nearest_kd(metricData, tree, refs, x, y)
    else
        error("Unknown problem type: $problem")
    end
end
