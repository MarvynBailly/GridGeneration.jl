# ODE Numerical Methods
We wish to solve the ODE

$\boxed{8 \sigma^4  M^2 x_s^2 x_{ss}  + 4 \sigma^4  M M_x x_s^4 + 4\sigma^2 m^2 M x_{ss} + 2 \sigma^2 m M_x x_s^2  = 0.}$

## Numerical Solver 
### DifferentialEquations.jl
Let's use [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) to try and solve the ODE with low effort.

Let $u_1 = x$ and $u_2 = x_s$ to find

$\begin{cases} 
u_1' = u_2 = x',\\
u'_2 = x'' = - \frac{1}{2} \frac{M_x}{M}.
\end{cases}$

with dirichlet boundary conditions:

$u_1(1) = x_1, u_1(0) = x_0.$

    



## Semi-Analytical Solution
