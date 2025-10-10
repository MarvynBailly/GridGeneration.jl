
"""
- problem 1: Constant metric
- problem 2: leading edge
- problem 3: trailing edge
- problem 4: leading edge and trailing edge
- problem 5: for fun
"""
function GetRectangleMetric(problem = 1; scale = 4000)
    if problem == 1
        M11 = scale
        M22 = scale
        return [M11, M22]
    elseif problem == 2
        metricFunction =  make_getMetric(airfoil;
                A_airfoil = 100,  ℓ_airfoil = 0.5, p_airfoil = 2,   
                A_origin  = scale,  ℓ_origin  = 0.1, p_origin  = 10,   
                floor     = 1e-4,  origin_center=(0, 0),
                profile   = :rational
            ) 
        return metricFunction(x,y)
    elseif problem == 3
        metricFunction =  make_getMetric(airfoil;
                A_airfoil = 100.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
                A_origin  = scale,  ℓ_origin  = 0.1, p_origin  = 10,   
                floor     = 1e-4,  origin_center=(1, 0),
                profile   = :rational
            ) 
        return metricFunction(x,y)
    elseif problem == 4
        metricFunction1 =  make_getMetric(airfoil;
                A_airfoil = 100.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
                A_origin  = scale,  ℓ_origin  = 0.1, p_origin  = 10,   
                floor     = 1e-4,  origin_center=(0, 0),
                profile   = :rational
            ) 
        metricFunction2 =  make_getMetric(airfoil;
                A_airfoil = 100.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
                A_origin  = scale,  ℓ_origin  = 0.1, p_origin  = 10,   
                floor     = 1e-4,  origin_center=(1, 0),
                profile   = :rational
            )

        return metricFunction1(x,y) .+ metricFunction2(x,y)
    elseif problem == 5
        metricFunction1 =  make_getMetric(airfoil;
                A_airfoil = 100.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
                A_origin  = 1000.0,  ℓ_origin  = 0.1, p_origin  = 10,   
                floor     = 1e-4,  origin_center=(0, 0),
                profile   = :rational
            ) 
        metricFunction2 =  make_getMetric(airfoil;
                A_airfoil = 100.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
                A_origin  = 1000.0,  ℓ_origin  = 0.1, p_origin  = 10,   
                floor     = 1e-4,  origin_center=(1, 0),
                profile   = :rational
            )
        metricFunction3 =  make_getMetric(airfoil;
                A_airfoil = 100.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
                A_origin  = 1000.0,  ℓ_origin  = 0.1, p_origin  = 10,   
                floor     = 1e-4,  origin_center=(0.30, 0.13),
                profile   = :rational
            ) 
        metricFunction4 =  make_getMetric(airfoil;
                A_airfoil = 100.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
                A_origin  = 1000.0,  ℓ_origin  = 0.1, p_origin  = 10,   
                floor     = 1e-4,  origin_center=(0.50, -0.09),
                profile   = :rational
            )

        return metricFunction1(x,y) .+ metricFunction2(x,y) .+ metricFunction3(x,y) .+ metricFunction4(x,y)
    else
        error("Unknown problem type: $problem")
    end
end