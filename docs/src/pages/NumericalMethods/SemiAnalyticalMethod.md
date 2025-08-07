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
function SemiAnalyticSolution(x, M, N)
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

    

```