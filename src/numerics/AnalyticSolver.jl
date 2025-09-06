function SolveAnalytic(xs, M, N; method="linear")
    x_sol = equidistribute_linear_inverse(xs, M, N)
    sol = zeros(2, N)
    sol[1, :] = x_sol
    return sol
end


"""
    semi_analytic_solution_from_data(x_vals, M_vals, N)
Compute a semi-analytic solution for the second-order BVP using the provided metric values `M_vals` at points `x_vals`. 
Uses numerical integration and interpolation to derive the solution.
The function returns a tuple of `s_vals` (the normalized parameter values) and `x_sol` (the corresponding x values).
"""
function equidistribute_linear_inverse(x, m, N)
    @assert length(x) == length(m) "x and m must have same length"
    @assert isapprox(first(x), 0.0; atol=1e-12) && isapprox(last(x), 1.0; atol=1e-12)
    @assert issorted(x) "x must be strictly increasing"
    @assert minimum(m) >= 0 "m must be positive"

    n = length(x) - 1
    Δx = x[2:end] .- x[1:end-1]
    w  = sqrt.(m)  # √m

    # cumulative trapezoid
    I = similar(x, Float64)
    I[1] = 0.0
    for i in 1:n
        I[i+1] = I[i] + 0.5*(w[i] + w[i+1]) * Δx[i]
    end
    Itot = I[end]

    # uniform s-grid and targets
    s = range(0.0, 1.0; length=N)
    t = Itot .* s

    # invert by linear interpolation on (I, x)
    x_nodes = similar(s, Float64)
    for (j, tj) in enumerate(t)
        # find i with I[i] <= tj <= I[i+1]
        i = searchsortedlast(I, tj)
        if i == length(I)         # tj == I[end]
            x_nodes[j] = x[end]
        elseif I[i+1] == I[i]     # flat segment (m ~ 0)
            x_nodes[j] = x[i]     # or x[i] + θ*Δx[i] if you want uniform spread
        else
            θ = (tj - I[i]) / (I[i+1] - I[i])
            x_nodes[j] = x[i] + θ * (x[i+1] - x[i])
        end
    end
    return x_nodes
end