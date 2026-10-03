"""
    GetRectangleMetric(problem = 1; scale = 4000)

Return a metric function `M(x, y) -> (M11, M22)` for the rectangle example.

- problem 1: constant metric
- problem 2: hotspot at (0, 0)
- problem 3: hotspot at (1, 0)
- problem 4: hotspots at (0, 0) and (1, 0)
- problem 5: four hotspots, for fun
"""
function GetRectangleMetric(problem = 1; scale = 4000)
    hotspot(center; A = scale) = GridGeneration.make_getMetric(nothing;
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
    else
        error("Unknown problem type: $problem")
    end
end
