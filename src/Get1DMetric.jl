""" 
Get1DMetric.jl
- Description: Compute metric between points depending on method 
- Input: 2xn array of points 
- Output: 1xn array of metric 
"""
function Get1DMetric(points, getMetric; method = "local")

    function norm(v)
        return sqrt(v[1]^2 + v[2]^2)
    end


    n = size(points, 2)
    m_vals = zeros(Float64, n)
    diff = zeros(Float64, 2, n)

    diff[:, 2:n-1] = (points[:, 3:n] - points[:, 1:n-2])     
    diff[:, n] = (points[:, n] - points[:, n-1]) 
    diff[:,1] = points[:, 2] - points[:, 1] 
    
    for i in 1:n
        # get metric value for the points M
        metricValues = getMetric(points[1, i], points[2, i])
        M = zeros(Float64, 2, 2)
        M[1, 1] = metricValues[1]
        M[2, 2] = metricValues[2]
        
        
        localDiff = diff[:, i]
        if method == "local"
            # normalize the localDiff vector
            normLocalDiff = norm(localDiff)
            m_vals[i] = localDiff' * M * localDiff / normLocalDiff^2
        elseif method == "nuclear"
            m_vals[i] = M[1, 1] + M[2, 2]
        end
    end


    return m_vals
end