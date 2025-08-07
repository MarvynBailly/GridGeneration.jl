# Second Order Nonlinear ODE Boundary Value Solver
we have the ODE

$x_{ss} + \frac{ M_x x_s^2}{2M} = 0,$

with boundary values $x(0) = 0$ and $x(1) = 1$ from [math work](../ODE/MathematicalWork.md). Let $f(x) = \frac{M_x}{2 M}$ to get

$x_{ss} + f(x) x_{s}^2 = 0$.

## Numerical Method
Let's begin by discretizing $x_{ss}$ and $x_s$ via a second order central difference:

$x_{s} = \frac{x_{i+1} - x_{i-1}}{2h}, \quad x_{ss} = \frac{x_{i-1} - 2x_i + x_{i+1}}{h^2},$ 

where $h = \Delta x$. Plugging this into the ODE yields the ith equation to be

$\frac{1}{2h} \left(x_{i-1} - 2x_i + x_{i+1}\right) + \frac{f_i}{4 h^2} \left(x_{i+1} - x_{i-1} \right)^2 = 0$.

Now let's use an iterative method to solve for the nonlinearity. Now on the $k+1$th iteration, we treat the nonlinear terms as known from the $k$th iteration to get

$\frac{1}{2h} \left(x_{i-1}^{(k+1)} - 2x_i^{(k+1)} + x_{i+1}^{(k+1)}\right) = - \frac{f_i^{(k)}}{4 h^2} \left(x_{i+1}^{(k)} - x_{i-1}^{(k)} \right)^2$.

This is a tri-diagonal system which we can solve with [Thomas's algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm) on each iteration.

## Algorithms
### Thomas Algorithm
A very straightforward guy

```julia
function thomas_algorithm(a, b, c, d)
    n = length(d)
    cp = similar(c)
    dp = similar(d)

    # forward sweep
    cp[1] = c[1] / b[1]
    dp[1] = d[1] / b[1]

    for i in 2:n
        denom = b[i] - a[i] * cp[i-1]
        cp[i] = c[i] / denom
        dp[i] = (d[i] - a[i] * dp[i-1]) / denom
    end

    # backward substitution
    x = zeros(n)
    x[end] = dp[end]
    for i = n-1:-1:1
        x[i] = dp[i] - cp[i] * x[i+1]
    end

    return x
end
```

### Solver

With this in hand, we set up the solver

```julia
function solve_bvp_fixed_point(f, x0, x1; N=100, max_iter=100, tol=1e-8)
    h = 1.0 / (N + 1)
    x = range(0, 1, length=N+2)
    u = [x0; collect(x[2:end-1]); x1]  # initial guess: linear

    u_new = similar(u)

    for iter in 1:max_iter
        rhs = zeros(N)

        for i in 2:N+1
            fi = f(u[i])
            du = u[i+1] - u[i-1]
            rhs[i-1] = - (fi / (4h^2)) * du^2
        end

        # Build tridiagonal matrix
        a = fill(1/h^2, N)    # subdiagonal
        b = fill(-2/h^2, N)   # diagonal
        c = fill(1/h^2, N)    # superdiagonal

        # Solve linear system
        δu = thomas_algorithm(a, b, c, rhs)

        # Update interior points
        u_new[1] = x0
        u_new[end] = x1
        u_new[2:N+1] = δu

        # Check convergence
        if norm(u_new - u, Inf) < tol
            println("Converged in $iter iterations")
            return x, u_new
        end

        u .= u_new
    end

    error("Did not converge after $max_iter iterations.")
end
```

Note here that $f$ is expecting a function. We can use our piecewise linear interpolator to take the discrete points given by $M_x(x_i)/(2M(x_i))$ where $M_x$ is approximated through central differences.

## Results