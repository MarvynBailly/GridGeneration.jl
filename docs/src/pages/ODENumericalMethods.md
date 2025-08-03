# ODE Numerical Methods
We wish to solve the ODE

$\boxed{8 \sigma^4  M^2 x_s^2 x_{ss}  + 4 \sigma^4  M M_x x_s^4 + 4\sigma^2 m^2 M x_{ss} + 2 \sigma^2 m M_x x_s^2  = 0.}$

## Numerical Solver 
### DifferentialEquations.jl
Let's use [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) to try and solve the ODE with low effort.

Let $u_1 = x$ and $u_2 = x_s$ to find

$\begin{cases} 
u_1' = u_2 = x',\\
u'_2 = x'' = - \frac{1}{2} \frac{M_x}{M} x_s^2.
\end{cases}$

with dirichlet boundary conditions:

$u_1(1) = x_1, u_1(0) = x_0.$

We can code this in Julia with the following:

```julia
# Spacing ODE
function SpacingODE!(du, u, p, s)
    M_u1_func, M_func, sigma, _, _ = p
    u1, u2 = u

    M_u1 = M_u1_func(u1)
    M = M_func(u1)

    du[1] = u2
    du[2] = - (M_u1 * u2^2) / (2 * M)    
end


# Set the boundary conditions
function BoundaryConditions!(residual, u, p, s)
    _, _, _, x0, x1 = p
    residual[1] = u[1][1] - x0
    residual[2] = u[end][1] - x1
end
```
Now as the user we can define the metric function and derivative (analytically for now) and solve the spacing

```@example
using GridGeneration

# define metric and derivative
scale = 8000
M(x) = scale * (1 + 15 * (x))^(-2)
M_u1(x) = -2 * scale * (1 + 15 * (x))^(-3) * (15)

# boundary values
x0 = 0.0
x1 = 1.0

# define number of grid points
N = 100

# pass to the numerical solver
sol = GridGeneration.SolveODE(M, M_u1, N, x0, x1);
``` 

## Semi-Analytical Solution


## Results
Let's compare the distribution of points for four different metrics
- Uniform
- Clustering at $x=0$
- Clustering at $x=1$
- Clustering at $x=0.5$
  
Here are the results for the numerical approach

### Uniform

![Uniform](../assets/images/ODENumericalMethods/Uniform_N50_numeric.svg)

### Clustering Near $x=0$

![x=0 clustering](../assets/images/ODENumericalMethods/x=0_N50_numeric.svg)

### Clustering Near $x=1$

![x=1 clustering](../assets/images/ODENumericalMethods/x=1_N50_numeric.svg)

