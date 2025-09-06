function Metric(x,y, scale, problem)
    if problem == 1
        M11 = scale * (1 .+ 15 * (x)).^(-2)
        M22 = scale * (1 .+ 15 * (x)).^(-2)
        return [M11, M22]
    elseif problem == 2
        M11 = scale * (1 .+ 15 * (1 .- x)).^(-2)
        M22 = 0 #scale * (1 .+ 15 * (1 .- y)).^(-2)
        return [M11, M22]
    elseif problem == 3
        M11 = scale
        M22 = scale
        return [M11, M22]
    end
end

function MetricDerivative(x, y, scale, problem)
    if problem == 1
        M11 = -2 * scale * (1 .+ 15 * (x)).^(-3) * (15)
        M22 =  -2 * scale * (1 .+ 15 * (y)).^(-3) * (15)
        return [M11, M22]
    elseif problem == 2
        M11 = -2 * scale * (1 .+ 15 * (1 .- x)).^(-3) * (-15)
        M22 = 0 # -2 * scale * (1 .+ 15 * (1 .- y)).^(-3) * (-15)
        return [M11, M22]
    elseif problem == 3
        M11 = 0
        M22 = 0 # 0
        return [M11, M22]
    end
end