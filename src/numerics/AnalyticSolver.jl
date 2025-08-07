function make_linear_inverse_interpolator(x_vals, y_vals)
    @assert length(x_vals) == length(y_vals)
    @assert all(diff(y_vals) .> 0) "y_vals must be strictly increasing for invertibility"

    function Iinv(y_query)
        @assert y_query ≥ y_vals[1] && y_query ≤ y_vals[end] "Query out of bounds"

        # Binary search to find the interval [i, i+1] such that y_vals[i] ≤ y_query ≤ y_vals[i+1]
        low, high = 1, length(y_vals) - 1
        while low ≤ high
            mid = div(low + high, 2)
            if y_vals[mid] ≤ y_query ≤ y_vals[mid+1]
                # linear interpolation
                y1, y2 = y_vals[mid], y_vals[mid+1]
                x1, x2 = x_vals[mid], x_vals[mid+1]
                t = (y_query - y1) / (y2 - y1)
                return x1 + t * (x2 - x1)
            elseif y_query < y_vals[mid]
                high = mid - 1
            else
                low = mid + 1
            end
        end

        error("Value not found in data range.")
    end

    return Iinv
end


"""
    semi_analytic_solution_from_data(x_vals, M_vals, N)
Compute a semi-analytic solution for the second-order BVP using the provided metric values `M_vals` at points `x_vals`. 
Uses numerical integration and interpolation to derive the solution.
The function returns a tuple of `s_vals` (the normalized parameter values) and `x_sol` (the corresponding x values).
"""
function SolveAnalytic(x_vals, M_vals, N)
    @assert length(x_vals) == length(M_vals) "x_vals and M_vals must be the same length"
    @assert all(diff(x_vals) .>= 0) "x_vals must be increasing"

    # Step 1: Compute I(x) via trapezoidal integration of sqrt(M)
    sqrtM = sqrt.(M_vals)
    I_vals = zeros(length(x_vals))
    for i in 2:length(x_vals)
        Δx = x_vals[i] - x_vals[i-1]
        I_vals[i] = I_vals[i-1] + 0.5 * Δx * (sqrtM[i] + sqrtM[i-1])
    end

    I_total = I_vals[end]

    # Step 2: Interpolation of I(x) to get inverse: I(x) → x
    # Itp = LinearInterpolation(I_vals, x_vals, extrapolation_bc=Throw())
    Iinv = make_linear_inverse_interpolator(x_vals, I_vals)

    # Step 3: Compute x(s) from I(x(s)) = s * I_total
    s_vals = range(0, 1; length=N)
    x_of_s = [Iinv(s * I_total) for s in s_vals]

    return s_vals, x_of_s
end
