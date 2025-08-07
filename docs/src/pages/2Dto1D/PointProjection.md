# Projecting Points
## 2D to 1D Projection

Suppose we wish to solve the grid spacing along a discrete boundary $\Gamma$ given by the points $\gamma_i \in \R^2$ for $i=1,2,\dots, n$ where $n$ is the total number of points along the boundary. In-between each point is a linear interpolation $\Gamma_i$ for $i=1,2,\dots,n-1$ which defines the piecewise-continuous boundary $\Gamma$. That is

$\Gamma(x) = \begin{cases} \Gamma_i(\eta_i(x)) &\text{if } x \in [\gamma_i, \gamma_{i+1}], \\ 0 &\text{otherwise}, \end{cases}$

where

$\Gamma_i(\eta) = \eta \gamma_{i+1} + (1 - \eta) \gamma_i, \quad \eta \in [0,1].$

Here $\eta_i(x)$ would be a parameterization between $\gamma_i$ and $\gamma_{i+1}$. That is 

$\eta_i(x) = \frac{|x - \gamma_i|}{|\gamma_{i+1} - \gamma_i|}.$ 

Okay writing this out mathematically is a pain, let's just explain the algorithm.


## Algorithm

Given $\gamma_i \in \R^2$, find a 1D representation $s \in \R$. We want to do this in such a way that we can find an invertible mapping $f: \R^2 \to \R$. Now with $f(\Gamma_i) \in \R$ and $m \in \R$, we can solve the spacing ODE to get the optimal distribution $\xi_i$. Finally we can project the points back onto the $\Gamma_i$ to get the optimal distribution on the boundary using $f^{-1}$.

### Basic Method in Words
#### Step 1 - 2D to 1D
Given points $\{ \gamma_i \}_{i=1}^n$ we can $\{ f(\gamma_i) \}_{i=1}^n$ by 

- Computing the local spacing in array `diff`,
- Accumulate the spacing `sp = cumsum(diff)`,
- Add zero back in `xi = [0, sp]`,
- Normalize `xs = xi / xi[end]`.

Now $\{ \text{xs} \}_{i=1}^n$ is within $[0,1]$ with normalized spacing according to the original boundary spacing. Algorithmically:

```julia
function Boundary2Dto1D(boundary)
    # boundarySection is 2×N (rows: x,y; columns: points in order)
    x = boundary[1, :]
    y = boundary[2, :]

    # segment lengths (N-1)
    Δx = diff(x)
    Δy = diff(y)
    Δs = sqrt.(Δx.^2 .+ Δy.^2)

    xs = [0.0; cumsum(Δs)]   # length N, xs[1]=0, xs[end]=arclength

    # normalize 
    xs = xs ./ xs[end]  # now xs is in [0, 1]
    return xs
end
```

#### Step 2 - Spacing
Now we get the optimal spacing

- Compute non-optimal solution `sol = GridGeneration.SolveODE(M, Mx, N, xs[1], xs[end])`
- Compute optimal spacing `N_opt = ceil(Int, 1 / GridGeneration.ComputeOptimalSpacing(sol[1, :], M, xs))`
- Compute optimal solution `sol_opt = GridGeneration.SolveODE(M, Mx, N_opt, xs[1], xs[end])`  

```julia
function GetOptimalSolution(m, mx, N, xs; method = "system of odes")
    m_func = GridGeneration.build_interps_linear(xs, m)
    mx_func = GridGeneration.build_interps_linear(xs, mx)

    if method == "system of odes"
        sol = GridGeneration.SolveODE(m_func, mx_func, N, xs[1], xs[end])
        N_opt = GridGeneration.ComputeOptimalNumberofPoints(sol[1, :], m, xs)
        @info("Optimal number of points: ", N_opt)
        sol_opt = GridGeneration.SolveODE(m_func, mx_func, N_opt, xs[1], xs[end])
    end

    return sol_opt, sol
end
```


#### Step 3 - 1D to 2D
And finally let's project the points back onto the boundary

- Extract solution `x_sol = sol_opt[1,:]`
- For each `x` in `x_sol`
  - determine `i` such that $x \in [\text{xs}_i, \text{xs}_{i+1}]$
  - compute $\eta = \frac{|x - \text{xs}_i|}{|\text{xs}_{i+1} - \text{xs}_i|}$
  - Interpolate `x` onto $\Gamma_i$ according to $\eta$

The algorithm for this can be pretty straightforward if we assert that the incoming `xs` is increasing. We have already taken care of that. Let's also add clamping to the edges to make sure that the start and ending points don't move.

```julia
function ProjectBoundary1Dto2D(boundary, points, xs)
    projectedPoints = zeros(2, length(xs))
    for (i, pnt) in enumrate(points)
        intervalIndex = FindContainingIntervalIndex(pnt, xs)
        dist = pnt - xs[intervalIndex]
        normalDist = dist / (xs[intervalIndex + 1] - xs[intervalIndex])
        projectPoint = ProjectPointOntoBoundary(normalDist, intervalIndex, boundary)
        projectedPoints[i] = projectPoint
    end
    return projectedPoints
end
```

with the two helper functions

```julia
function FindContainingIntervalIndex(pnt, dist)
    # assert that dist is increasing
    index = - 1
    
    # clamp
    if pnt < dist[1]
        index = 1
    end

    if pnt > dist[end]
        index = length(dist) - 1
    end

    for (i,d) in enumerate(dist)
        if d > pnt
            index = i - 1
        end
    end

    return index
end

function ProjectPointOntoBoundary(s, ind, boundary)
    # interpolate boundary[ind] and boundary[ind + 1] with s
    projectedPoint = s * boundary[ind] + (1 - s) * boundary[ind+1]
    return projectPoint
end

```

