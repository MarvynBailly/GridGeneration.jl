# Semi Analytic Solution

## Set up 

As shown in [math work](../ODE/MathematicalWork.md), we can reduce the ODE to

We have an initial value boundary value problem (IVBP ODE), so let's prescribe boundary conditions at $x(s=0) = 0$ and $x(s=1) = 1$ such that our computational domain is $s \in [0,1]$ and the physical domain is $x \in [0,1] \subset \R$.

$\int_{x(0)}^{x(s)} \sqrt{M(\xi)} d\xi = C_1 s$

Finally let's enforce the boundary condition at $x(1) = 1$ to solve for $C_1$ and find that

$C_1 = \int_{0}^{1} \sqrt{M(\xi)} d \xi = I.$

Therefore the final solution becomes

$\int_{0}^{x(s)} \sqrt{M(\xi)} d\xi = I s.$

If we let $I(x) = \int_0^x \sqrt{M(\xi)} d \xi$, we can further clean up the express as

$I(x(s)) = s I(x(s=1))$

## Algorithm
Since the real function of $M(x)$ is unknown, let's use trapezoid rule to integrate. Then the algorithm will follow the steps

- Compute $I(x)$ via trapezoid rule.
- Compute $l = I(1)$.
- Solve $I(x) = s * I(1)$
  - We can do this via interpolation $x(s) = I^{-1}(s \cdot I(1))$.


```julia
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

function semi_analytic_solution_from_data(x_vals, M_vals, N)
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
```

## Results