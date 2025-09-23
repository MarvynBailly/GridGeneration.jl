"""
    ComputeOptimalSpacing(x, M, s)
Compute the optimal grid spacing based on the metric values `M` at points `x` and spacing `s`.
Integration is done using the trapezoidal rule.
return Ceil(Int, 1 / sigma_opt) as the optimal number of points.
"""
function ComputeOptimalNumberofPoints(x, M)
    Nn = length(x)
    # set up computational coordinate
    s = range(0, 1, length=Nn)

    # compute x_s
    x_s = zeros(Nn)
    @inbounds for i in 2:Nn-1
        Δs = s[i+1] - s[i-1]
        x_s[i] = (x[i+1] - x[i-1]) / Δs
    end
    x_s[1] = (x[2] - x[1]) / (s[2] - s[1])               # forward difference
    x_s[end] = (x[end] - x[end-1]) / (s[end] - s[end-1]) # backward difference

    # compute integral with trapezoidal rule
    numer = 0.0
    denom = 0.0
    @inbounds for i in 2:Nn-1
        Δs   = s[i+1] - s[i]
        p    = 0.5*(M(x[i-1])*x_s[i-1]^2 + M(x[i])*x_s[i]^2)
        numer += p * Δs
        denom += (p^2) * Δs
    end
    # println("numer: $numer, denom: $denom")

    sigma_opt = sqrt(numer / denom)

    N_opt = floor(Int, 1/(sigma_opt))
    return N_opt
end