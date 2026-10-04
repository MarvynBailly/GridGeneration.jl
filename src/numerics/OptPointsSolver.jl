"""
Compute optimal number of grid points based on metric field.
Uses trapezoidal integration to evaluate σ_opt = sqrt(∫p ds / ∫p² ds).
Returns floor(1/σ_opt).
"""
function ComputeOptimalNumberofPoints(x, M)
    Nn = length(x)
    s = range(0, 1, length=Nn)

    # compute derivatives x_s using central differences
    x_s = zeros(Nn)
    @inbounds for i in 2:Nn-1
        Δs = s[i+1] - s[i-1]
        x_s[i] = (x[i+1] - x[i-1]) / Δs
    end
    x_s[1] = (x[2] - x[1]) / (s[2] - s[1])
    x_s[end] = (x[end] - x[end-1]) / (s[end] - s[end-1])

    # trapezoidal integration
    numer = 0.0
    denom = 0.0
    p_prev = M(x[1])*x_s[1]^2
    @inbounds for i in 2:Nn
        Δs   = s[i] - s[i-1]
        p    = M(x[i])*x_s[i]^2
        numer += 0.5*(p_prev + p) * Δs
        denom += 0.5*(p_prev^2 + p^2) * Δs
        p_prev = p
    end

    (numer > 0 && denom > 0) || throw(ArgumentError("cannot compute the optimal number of points: metric integrals must be positive (got ∫p = $numer, ∫p² = $denom)"))

    sigma_opt = sqrt(numer / denom)
    # at least 3 points so one-sided second-order stencils remain valid
    N_opt = max(3, floor(Int, 1/(sigma_opt)))
    return N_opt
end