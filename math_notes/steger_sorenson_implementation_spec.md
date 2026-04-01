# Implementation Specification: Steger & Sorenson Elliptic Grid Smoother

A complete, self-contained specification for implementing the Steger & Sorenson elliptic
grid smoothing method in Julia. This document contains every equation, algorithm, data
structure, and function signature needed to produce a working implementation for a
single-block structured grid.

---

## 1. Scope and Interface

### 1.1 What This Solver Does

Takes a structured 2D grid (from TFI, edge solving, or any other source) and smooths the
interior points by solving an elliptic PDE system. Boundary points are held fixed.
Optionally enforces grid orthogonality at selected boundaries.

### 1.2 Function Signature

```julia
function StegerSorensonSolver(
    x::Matrix{Float64},       # Ni x Nj matrix of x-coordinates
    y::Matrix{Float64};       # Ni x Nj matrix of y-coordinates
    params::SSParams           # Parameters struct (see Section 2)
) -> (x, y, error, iterations)
```

**Input grid format**:
- `x[i, j]` and `y[i, j]` are $N_i \times N_j$ matrices
- $i$ indexes the $\xi$-direction ($i = 1$ is the left boundary, $i = N_i$ is the right)
- $j$ indexes the $\eta$-direction ($j = 1$ is the bottom boundary, $j = N_j$ is the top)
- Boundary points ($i = 1$, $i = N_i$, $j = 1$, $j = N_j$) are **fixed** throughout the solve

**Output**:
- `x, y`: the smoothed coordinate matrices (same size, boundaries unchanged)
- `error`: final convergence error (Float64)
- `iterations`: number of iterations performed (Int)

### 1.3 Integration with GridGeneration.jl

The solver replaces the inner call in `SmoothBlocks` (`src/smoothing/SmoothBlocks.jl`).
The 3D block format `blocks[k][2, Ni, Nj]` is unpacked to 2D matrices before calling
this solver:

```julia
x = blocks[k][1, :, :]    # Ni x Nj
y = blocks[k][2, :, :]    # Ni x Nj
x, y, err, iters = StegerSorensonSolver(x, y; params=params[k])
blocks[k] = permutedims(cat(x, y, dims=3), (3,1,2))  # back to [2, Ni, Nj]
```

---

## 2. Parameters Struct

```julia
struct SSParams
    # Iteration control
    max_iter::Int              # Maximum iterations (default: 5000)
    tol::Float64               # Convergence tolerance (default: 1e-6)
    omega::Float64             # SOR relaxation factor (default: 1.2)

    # Wall orthogonality forcing
    useBottomWall::Bool        # Enforce orthogonality at j = 1 (default: false)
    useTopWall::Bool           # Enforce orthogonality at j = Nj (default: false)
    useLeftWall::Bool          # Enforce orthogonality at i = 1 (default: false)
    useRightWall::Bool         # Enforce orthogonality at i = Ni (default: false)

    # Forcing decay rates (one per wall)
    decay_bottom::Float64      # Exponential decay rate for bottom wall (default: 0.4)
    decay_top::Float64         # Exponential decay rate for top wall (default: 0.4)
    decay_left::Float64        # Exponential decay rate for left wall (default: 0.4)
    decay_right::Float64       # Exponential decay rate for right wall (default: 0.4)

    # Diagnostics
    verbose::Bool              # Print convergence info (default: false)
    print_interval::Int        # Iterations between status prints (default: 500)

    # Skip flag
    skipBlock::Bool            # If true, return the input grid unchanged (default: false)
end
```

**Default constructor**:

```julia
function SSParams(;
    max_iter=5000, tol=1e-6, omega=1.2,
    useBottomWall=false, useTopWall=false,
    useLeftWall=false, useRightWall=false,
    decay_bottom=0.4, decay_top=0.4,
    decay_left=0.4, decay_right=0.4,
    verbose=false, print_interval=500,
    skipBlock=false
)
    return SSParams(
        max_iter, tol, omega,
        useBottomWall, useTopWall, useLeftWall, useRightWall,
        decay_bottom, decay_top, decay_left, decay_right,
        verbose, print_interval, skipBlock
    )
end
```

### 2.1 Parameter Guidance

| Parameter | Conservative | Moderate | Aggressive |
|-----------|-------------|----------|------------|
| $\omega$ | $0.8 - 1.0$ | $1.0 - 1.4$ | $1.4 - 1.8$ |
| decay | $0.1$ | $0.4$ | $1.0$ |
| tol | $10^{-8}$ | $10^{-6}$ | $10^{-4}$ |

- Start with $\omega = 1.0$ (pure Gauss-Seidel). Increase toward 1.4 once convergence
  is confirmed.
- If diverging, reduce $\omega$ toward 0.8.
- $\text{decay} \sim 0.1$ forces orthogonality deep into the domain;
  $\text{decay} \sim 1.0$ only near the wall.

---

## 3. Metric Computation

### 3.1 Function: `compute_metrics`

Computes first derivatives and metric coefficients at every grid point.

```julia
function compute_metrics(x::Matrix, y::Matrix)
    Ni, Nj = size(x)

    x_xi  = zeros(Ni, Nj);  y_xi  = zeros(Ni, Nj)
    x_eta = zeros(Ni, Nj);  y_eta = zeros(Ni, Nj)

    # --- xi-derivatives ---
    for j in 1:Nj
        # Left boundary: forward difference
        x_xi[1, j] = x[2, j] - x[1, j]
        y_xi[1, j] = y[2, j] - y[1, j]
        # Right boundary: backward difference
        x_xi[Ni, j] = x[Ni, j] - x[Ni-1, j]
        y_xi[Ni, j] = y[Ni, j] - y[Ni-1, j]
        # Interior: central difference
        for i in 2:Ni-1
            x_xi[i, j] = (x[i+1, j] - x[i-1, j]) / 2.0
            y_xi[i, j] = (y[i+1, j] - y[i-1, j]) / 2.0
        end
    end

    # --- eta-derivatives ---
    for i in 1:Ni
        # Bottom boundary: forward difference
        x_eta[i, 1] = x[i, 2] - x[i, 1]
        y_eta[i, 1] = y[i, 2] - y[i, 1]
        # Top boundary: backward difference
        x_eta[i, Nj] = x[i, Nj] - x[i, Nj-1]
        y_eta[i, Nj] = y[i, Nj] - y[i, Nj-1]
        # Interior: central difference
        for j in 2:Nj-1
            x_eta[i, j] = (x[i, j+1] - x[i, j-1]) / 2.0
            y_eta[i, j] = (y[i, j+1] - y[i, j-1]) / 2.0
        end
    end

    # --- Metric coefficients ---
    alpha = x_eta.^2 + y_eta.^2
    beta  = x_xi .* x_eta + y_xi .* y_eta
    gamma = x_xi.^2 + y_xi.^2
    J     = x_xi .* y_eta - x_eta .* y_xi

    return alpha, beta, gamma, J, x_xi, y_xi, x_eta, y_eta
end
```

The metric coefficients are:

$$\alpha_{i,j} = (x_\eta)_{i,j}^2 + (y_\eta)_{i,j}^2$$

$$\beta_{i,j} = (x_\xi)_{i,j}\,(x_\eta)_{i,j} + (y_\xi)_{i,j}\,(y_\eta)_{i,j}$$

$$\gamma_{i,j} = (x_\xi)_{i,j}^2 + (y_\xi)_{i,j}^2$$

$$J_{i,j} = (x_\xi)_{i,j}\,(y_\eta)_{i,j} - (x_\eta)_{i,j}\,(y_\xi)_{i,j}$$

**Note**: This function returns the first-derivative arrays (`x_xi`, `y_xi`, `x_eta`,
`y_eta`) in addition to the metric coefficients. The wall-forcing computation needs
these derivatives, so returning them avoids redundant computation.

---

## 4. Control Function Computation (Wall Forcing)

This section implements the core Steger & Sorenson innovation: automatic computation of
control functions $P$ and $Q$ from boundary orthogonality conditions.

### 4.1 Function: `compute_control_functions`

Computes $P_{i,j}$ and $Q_{i,j}$ arrays over the entire domain by evaluating boundary
conditions at each enabled wall and decaying into the interior.

```julia
function compute_control_functions(
    x::Matrix, y::Matrix,
    alpha::Matrix, beta::Matrix, gamma::Matrix, J::Matrix,
    x_xi::Matrix, y_xi::Matrix, x_eta::Matrix, y_eta::Matrix;
    params::SSParams
)
    Ni, Nj = size(x)
    P = zeros(Ni, Nj)
    Q = zeros(Ni, Nj)

    if params.useBottomWall
        Pb, Qb = compute_PQ_eta_wall(x, y, x_xi, y_xi, x_eta, y_eta; wall=1)
        for j in 1:Nj
            decay = exp(-params.decay_bottom * abs(j - 1))
            for i in 1:Ni
                P[i,j] += Pb[i] * decay
                Q[i,j] += Qb[i] * decay
            end
        end
    end

    if params.useTopWall
        Pt, Qt = compute_PQ_eta_wall(x, y, x_xi, y_xi, x_eta, y_eta; wall=Nj)
        for j in 1:Nj
            decay = exp(-params.decay_top * abs(j - Nj))
            for i in 1:Ni
                P[i,j] += Pt[i] * decay
                Q[i,j] += Qt[i] * decay
            end
        end
    end

    if params.useLeftWall
        Pl, Ql = compute_PQ_xi_wall(x, y, x_xi, y_xi, x_eta, y_eta; wall=1)
        for i in 1:Ni
            decay = exp(-params.decay_left * abs(i - 1))
            for j in 1:Nj
                P[i,j] += Pl[j] * decay
                Q[i,j] += Ql[j] * decay
            end
        end
    end

    if params.useRightWall
        Pr, Qr = compute_PQ_xi_wall(x, y, x_xi, y_xi, x_eta, y_eta; wall=Ni)
        for i in 1:Ni
            decay = exp(-params.decay_right * abs(i - Ni))
            for j in 1:Nj
                P[i,j] += Pr[j] * decay
                Q[i,j] += Qr[j] * decay
            end
        end
    end

    return P, Q
end
```

The decay propagation for each wall is:

$$P_{i,j} = \sum_{\text{walls}} P_{\text{wall}} \cdot e^{-c\,|\,\text{dist from wall}\,|}$$

### 4.2 Function: `compute_PQ_eta_wall`

Computes $P$ and $Q$ at a bottom ($j{=}1$) or top ($j{=}N_j$) wall.

**Algorithm at each interior wall point $i = 2, \ldots, N_i{-}1$:**

**Step 1** -- Wall tangent derivatives (central differences along wall):

$$(x_\xi)_w = \frac{x_{i+1,w} - x_{i-1,w}}{2}, \qquad (x_{\xi\xi})_w = x_{i+1,w} - 2\,x_{i,w} + x_{i-1,w}$$

**Step 2** -- Wall-normal spacing (linear interpolation of corner spacings):

$$s_1 = \left\|\mathbf{r}_{1,w} - \mathbf{r}_{1,w+\text{dir}}\right\|, \qquad s_{N_i} = \left\|\mathbf{r}_{N_i,w} - \mathbf{r}_{N_i,w+\text{dir}}\right\|$$

$$s = s_1 + (s_{N_i} - s_1)\,\frac{i - 1}{N_i - 1}$$

**Step 3** -- Desired $\eta$-derivatives (orthogonal + prescribed spacing):

$$(x_\eta)_d = \frac{-s\,(y_\xi)_w}{\sqrt{(x_\xi)_w^2 + (y_\xi)_w^2}}, \qquad (y_\eta)_d = \frac{s\,(x_\xi)_w}{\sqrt{(x_\xi)_w^2 + (y_\xi)_w^2}}$$

**Step 4** -- Second normal derivatives (one-sided, using desired first derivative):

$$x_{\eta\eta} = \frac{1}{2}\!\left(-7\,x_{i,w} + 8\,x_{i,w+\text{dir}} - x_{i,w+2\text{dir}}\right) - \text{dir}\cdot 3\,(x_\eta)_d$$

$$y_{\eta\eta} = \frac{1}{2}\!\left(-7\,y_{i,w} + 8\,y_{i,w+\text{dir}} - y_{i,w+2\text{dir}}\right) - \text{dir}\cdot 3\,(y_\eta)_d$$

**Step 5** -- Boundary metrics (with desired derivatives):

$$\alpha_b = (x_\eta)_d^2 + (y_\eta)_d^2 = s^2, \qquad \gamma_b = (x_\xi)_w^2 + (y_\xi)_w^2$$

$$J_b = (x_\xi)_w\,(y_\eta)_d - (x_\eta)_d\,(y_\xi)_w = s\,\|(x_\xi, y_\xi)_w\|$$

**Step 6** -- Solve $2 \times 2$ system for $P$, $Q$:

$$\begin{bmatrix} (x_\xi)_w & (x_\eta)_d \\ (y_\xi)_w & (y_\eta)_d \end{bmatrix} \begin{bmatrix} P \\ Q \end{bmatrix} = \frac{-1}{J_b^2} \begin{bmatrix} \alpha_b\,(x_{\xi\xi})_w + \gamma_b\,x_{\eta\eta} \\ \alpha_b\,(y_{\xi\xi})_w + \gamma_b\,y_{\eta\eta} \end{bmatrix}$$

Solve via Cramer's rule with $\det = J_b$:

$$P = \frac{1}{J_b}\left(\frac{R_x}{J_b^2}\,(y_\eta)_d - \frac{R_y}{J_b^2}\,(x_\eta)_d\right), \qquad Q = \frac{1}{J_b}\left(\frac{R_y}{J_b^2}\,(x_\xi)_w - \frac{R_x}{J_b^2}\,(y_\xi)_w\right)$$

where $R_x = -\!\left(\alpha_b\,(x_{\xi\xi})_w + \gamma_b\,x_{\eta\eta}\right)$ and
$R_y = -\!\left(\alpha_b\,(y_{\xi\xi})_w + \gamma_b\,y_{\eta\eta}\right)$.

```julia
function compute_PQ_eta_wall(
    x::Matrix, y::Matrix,
    x_xi::Matrix, y_xi::Matrix,
    x_eta::Matrix, y_eta::Matrix;
    wall::Int
)
    Ni, Nj = size(x)
    P_wall = zeros(Ni)
    Q_wall = zeros(Ni)

    # Direction: +1 for bottom (j=1), -1 for top (j=Nj)
    dir = wall == 1 ? 1 : -1
    w = wall

    for i in 2:Ni-1
        # Step 1: Wall tangent derivatives
        xxi  = (x[i+1, w] - x[i-1, w]) / 2.0
        yxi  = (y[i+1, w] - y[i-1, w]) / 2.0
        xxixi = x[i+1, w] - 2.0*x[i, w] + x[i-1, w]
        yxixi = y[i+1, w] - 2.0*y[i, w] + y[i-1, w]

        # Step 2: Wall-normal spacing
        s1 = sqrt((x[1, w] - x[1, w + dir])^2 + (y[1, w] - y[1, w + dir])^2)
        s2 = sqrt((x[Ni, w] - x[Ni, w + dir])^2 + (y[Ni, w] - y[Ni, w + dir])^2)
        s = s1 + (s2 - s1) * (i - 1) / (Ni - 1)

        # Step 3: Desired eta-derivatives
        tangent_mag = sqrt(xxi^2 + yxi^2)
        xeta_d = -s * yxi / tangent_mag
        yeta_d =  s * xxi / tangent_mag

        # Step 4: Second normal derivatives
        xetaeta = 0.5*(-7.0*x[i,w] + 8.0*x[i,w+dir] - x[i,w+2*dir]) - dir*3.0*xeta_d
        yetaeta = 0.5*(-7.0*y[i,w] + 8.0*y[i,w+dir] - y[i,w+2*dir]) - dir*3.0*yeta_d

        # Step 5: Boundary metrics
        alpha_b = xeta_d^2 + yeta_d^2
        gamma_b = xxi^2 + yxi^2
        J_b     = xxi * yeta_d - xeta_d * yxi

        # Step 6: Solve 2x2 system
        rhs_x = -(alpha_b * xxixi + gamma_b * xetaeta)
        rhs_y = -(alpha_b * yxixi + gamma_b * yetaeta)

        J_b_sq = J_b^2
        if abs(J_b_sq) < 1e-30
            continue  # Degenerate point, skip
        end

        det_M = J_b
        P_wall[i] = (rhs_x/J_b_sq * yeta_d - rhs_y/J_b_sq * xeta_d) / det_M
        Q_wall[i] = (rhs_y/J_b_sq * xxi    - rhs_x/J_b_sq * yxi   ) / det_M
    end

    return P_wall, Q_wall
end
```

### 4.3 Function: `compute_PQ_xi_wall`

Computes $P$ and $Q$ at a left ($i{=}1$) or right ($i{=}N_i$) wall. Same logic as the
$\eta$-wall version but with $\xi$ and $\eta$ roles swapped.

**Desired $\xi$-derivatives** (orthogonal to $\eta$-tangent at wall):

$$(x_\xi)_d = \frac{-s\,(y_\eta)_w}{\sqrt{(x_\eta)_w^2 + (y_\eta)_w^2}}, \qquad (y_\xi)_d = \frac{s\,(x_\eta)_w}{\sqrt{(x_\eta)_w^2 + (y_\eta)_w^2}}$$

**Second normal derivatives** (note sign difference from $\eta$-wall):

$$x_{\xi\xi} = \frac{1}{2}\!\left(-7\,x_{w,j} + 8\,x_{w+\text{dir},j} - x_{w+2\text{dir},j}\right) + \text{dir}\cdot 3\,(x_\xi)_d$$

$$y_{\xi\xi} = \frac{1}{2}\!\left(-7\,y_{w,j} + 8\,y_{w+\text{dir},j} - y_{w+2\text{dir},j}\right) + \text{dir}\cdot 3\,(y_\xi)_d$$

```julia
function compute_PQ_xi_wall(
    x::Matrix, y::Matrix,
    x_xi::Matrix, y_xi::Matrix,
    x_eta::Matrix, y_eta::Matrix;
    wall::Int
)
    Ni, Nj = size(x)
    P_wall = zeros(Nj)
    Q_wall = zeros(Nj)

    dir = wall == 1 ? 1 : -1
    w = wall

    for j in 2:Nj-1
        # Step 1: Wall tangent derivatives
        xeta    = (x[w, j+1] - x[w, j-1]) / 2.0
        yeta    = (y[w, j+1] - y[w, j-1]) / 2.0
        xetaeta = x[w, j+1] - 2.0*x[w, j] + x[w, j-1]
        yetaeta = y[w, j+1] - 2.0*y[w, j] + y[w, j-1]

        # Step 2: Wall-normal spacing
        s1 = sqrt((x[w, 1] - x[w+dir, 1])^2 + (y[w, 1] - y[w+dir, 1])^2)
        s2 = sqrt((x[w, Nj] - x[w+dir, Nj])^2 + (y[w, Nj] - y[w+dir, Nj])^2)
        s = s1 + (s2 - s1) * (j - 1) / (Nj - 1)

        # Step 3: Desired xi-derivatives
        tangent_mag = sqrt(xeta^2 + yeta^2)
        xxi_d = -s * yeta / tangent_mag
        yxi_d =  s * xeta / tangent_mag

        # Step 4: Second normal derivatives (NOTE: +dir*3 for xi-walls)
        xxixi = 0.5*(-7.0*x[w,j] + 8.0*x[w+dir,j] - x[w+2*dir,j]) + dir*3.0*xxi_d
        yxixi = 0.5*(-7.0*y[w,j] + 8.0*y[w+dir,j] - y[w+2*dir,j]) + dir*3.0*yxi_d

        # Step 5: Boundary metrics
        alpha_b = xeta^2 + yeta^2
        gamma_b = xxi_d^2 + yxi_d^2
        J_b     = xxi_d * yeta - xeta * yxi_d

        # Step 6: Solve 2x2 system
        rhs_x = -(alpha_b * xxixi + gamma_b * xetaeta)
        rhs_y = -(alpha_b * yxixi + gamma_b * yetaeta)

        J_b_sq = J_b^2
        if abs(J_b_sq) < 1e-30
            continue
        end

        det_M = J_b
        P_wall[j] = (rhs_x/J_b_sq * yeta  - rhs_y/J_b_sq * xeta ) / det_M
        Q_wall[j] = (rhs_y/J_b_sq * xxi_d - rhs_x/J_b_sq * yxi_d) / det_M
    end

    return P_wall, Q_wall
end
```

### 4.4 Sign Convention for the One-Sided Second Derivative

The sign difference between the $\eta$-wall and $\xi$-wall formulas warrants explanation.

**$\eta$-wall (bottom, $\text{dir}{=}{+}1$)**: The one-sided first-derivative formula
pointing into the domain is:

$$f'_{\text{forward}} = \frac{-3\,f_{j=1} + 4\,f_{j=2} - f_{j=3}}{2}$$

The second derivative consistent with the desired $f'_d$ is:

$$f'' = \frac{1}{2}\!\left(-7\,f_1 + 8\,f_2 - f_3\right) - 3\,f'_d$$

**$\xi$-wall (left, $\text{dir}{=}{+}1$)**: Same structure but:

$$f'' = \frac{1}{2}\!\left(-7\,f_1 + 8\,f_2 - f_3\right) + 3\,f'_d$$

The sign rule: use $-\text{dir}\cdot 3\,f'_d$ for **$\eta$-walls** and
$+\text{dir}\cdot 3\,f'_d$ for **$\xi$-walls**. This matches the existing code in
`EllipticSolver.jl` and arises from the orientation convention of the outward normal
relative to the coordinate direction.

---

## 5. Thomas Algorithm (Tridiagonal Solver)

### 5.1 Function: `thomas_solve!`

Solves the tridiagonal system:

$$a_i\,u_{i-1} + b_i\,u_i + c_i\,u_{i+1} = d_i \qquad \text{for } i = 1, \ldots, n$$

in $O(n)$ operations. Arrays `a`, `b`, `c`, `d` are **modified in place**.

```julia
function thomas_solve!(a::Vector, b::Vector, c::Vector, d::Vector, u::Vector)
    n = length(b)

    # Forward elimination
    for i in 2:n
        w = a[i] / b[i-1]
        b[i] -= w * c[i-1]
        d[i] -= w * d[i-1]
    end

    # Back substitution
    u[n] = d[n] / b[n]
    for i in n-1:-1:1
        u[i] = (d[i] - c[i] * u[i+1]) / b[i]
    end
end
```

### 5.2 Numerical Stability

The Thomas algorithm is stable when the system is diagonally dominant:

$$|b_i| \geq |a_i| + |c_i| \qquad \forall\; i$$

For our elliptic system: $|b_i| = 2(\alpha + \gamma)$ and
$|a_i| + |c_i| = 2\alpha + |J^2 P|$. This is diagonally dominant as long as
$|J^2 P| \leq 2\gamma$, which is satisfied for reasonable grids and moderate $P$ values.

---

## 6. Line Gauss-Seidel Sweeps

### 6.1 Xi-Line Sweep (sweep along constant-$j$ lines)

For each $j$-level from $j{=}2$ to $j{=}N_j{-}1$, form and solve a tridiagonal system
for $x_{2:N_i-1,\,j}$ (and separately for $y$).

**Tridiagonal coefficients** at each interior point $i$:

$$a_i = \alpha_{i,j} - \frac{J_{i,j}^2\,P_{i,j}}{2}, \qquad b_i = -2\!\left(\alpha_{i,j} + \gamma_{i,j}\right), \qquad c_i = \alpha_{i,j} + \frac{J_{i,j}^2\,P_{i,j}}{2}$$

$$d_i = -\gamma_{i,j}\!\left(x_{i,j+1} + x_{i,j-1}\right) - \frac{J_{i,j}^2\,Q_{i,j}}{2}\!\left(x_{i,j+1} - x_{i,j-1}\right) + \frac{\beta_{i,j}}{2}\!\left(x_{i+1,j+1} - x_{i-1,j+1} - x_{i+1,j-1} + x_{i-1,j-1}\right)$$

Boundary terms are moved to the RHS: $d_1 \mathrel{{-}{=}} a_1 \cdot x_{1,j}$ and
$d_n \mathrel{{-}{=}} c_n \cdot x_{N_i,j}$.

After solving with Thomas algorithm, apply SOR:

$$x_{i,j} \leftarrow (1 - \omega)\,x_{i,j} + \omega\,u_i^{\text{Thomas}}$$

```julia
function xi_line_sweep!(
    x::Matrix, y::Matrix,
    alpha::Matrix, beta::Matrix, gamma::Matrix, J::Matrix,
    P::Matrix, Q::Matrix,
    omega::Float64
)
    Ni, Nj = size(x)
    n = Ni - 2

    # Pre-allocate work arrays
    a = zeros(n);  b = zeros(n);  c = zeros(n);  d = zeros(n);  u = zeros(n)

    for j in 2:Nj-1
        # --- Solve for x-coordinates ---
        for k in 1:n
            i = k + 1
            a[k] = alpha[i,j] - J[i,j]^2 * P[i,j] / 2.0
            b[k] = -2.0 * (alpha[i,j] + gamma[i,j])
            c[k] = alpha[i,j] + J[i,j]^2 * P[i,j] / 2.0
            d[k] = -gamma[i,j] * (x[i,j+1] + x[i,j-1]) -
                    J[i,j]^2 * Q[i,j] / 2.0 * (x[i,j+1] - x[i,j-1]) +
                    beta[i,j] / 2.0 * (x[i+1,j+1] - x[i-1,j+1] - x[i+1,j-1] + x[i-1,j-1])
        end
        d[1] -= a[1] * x[1, j];   a[1] = 0.0
        d[n] -= c[n] * x[Ni, j];  c[n] = 0.0

        thomas_solve!(a, b, c, d, u)
        for k in 1:n
            i = k + 1
            x[i, j] = (1.0 - omega) * x[i, j] + omega * u[k]
        end

        # --- Solve for y-coordinates (identical structure) ---
        for k in 1:n
            i = k + 1
            a[k] = alpha[i,j] - J[i,j]^2 * P[i,j] / 2.0
            b[k] = -2.0 * (alpha[i,j] + gamma[i,j])
            c[k] = alpha[i,j] + J[i,j]^2 * P[i,j] / 2.0
            d[k] = -gamma[i,j] * (y[i,j+1] + y[i,j-1]) -
                    J[i,j]^2 * Q[i,j] / 2.0 * (y[i,j+1] - y[i,j-1]) +
                    beta[i,j] / 2.0 * (y[i+1,j+1] - y[i-1,j+1] - y[i+1,j-1] + y[i-1,j-1])
        end
        d[1] -= a[1] * y[1, j];   a[1] = 0.0
        d[n] -= c[n] * y[Ni, j];  c[n] = 0.0

        thomas_solve!(a, b, c, d, u)
        for k in 1:n
            i = k + 1
            y[i, j] = (1.0 - omega) * y[i, j] + omega * u[k]
        end
    end
end
```

### 6.2 Eta-Line Sweep (sweep along constant-$i$ lines)

For each $i$-level from $i{=}2$ to $i{=}N_i{-}1$, form and solve a tridiagonal system
for $x_{i,\,2:N_j-1}$ (and separately for $y$).

**Tridiagonal coefficients** at each interior point $j$:

$$a_j = \gamma_{i,j} - \frac{J_{i,j}^2\,Q_{i,j}}{2}, \qquad b_j = -2\!\left(\alpha_{i,j} + \gamma_{i,j}\right), \qquad c_j = \gamma_{i,j} + \frac{J_{i,j}^2\,Q_{i,j}}{2}$$

$$d_j = -\alpha_{i,j}\!\left(x_{i+1,j} + x_{i-1,j}\right) - \frac{J_{i,j}^2\,P_{i,j}}{2}\!\left(x_{i+1,j} - x_{i-1,j}\right) + \frac{\beta_{i,j}}{2}\!\left(x_{i+1,j+1} - x_{i-1,j+1} - x_{i+1,j-1} + x_{i-1,j-1}\right)$$

```julia
function eta_line_sweep!(
    x::Matrix, y::Matrix,
    alpha::Matrix, beta::Matrix, gamma::Matrix, J::Matrix,
    P::Matrix, Q::Matrix,
    omega::Float64
)
    Ni, Nj = size(x)
    n = Nj - 2

    a = zeros(n);  b = zeros(n);  c = zeros(n);  d = zeros(n);  u = zeros(n)

    for i in 2:Ni-1
        # --- Solve for x-coordinates ---
        for k in 1:n
            j = k + 1
            a[k] = gamma[i,j] - J[i,j]^2 * Q[i,j] / 2.0
            b[k] = -2.0 * (alpha[i,j] + gamma[i,j])
            c[k] = gamma[i,j] + J[i,j]^2 * Q[i,j] / 2.0
            d[k] = -alpha[i,j] * (x[i+1,j] + x[i-1,j]) -
                    J[i,j]^2 * P[i,j] / 2.0 * (x[i+1,j] - x[i-1,j]) +
                    beta[i,j] / 2.0 * (x[i+1,j+1] - x[i-1,j+1] - x[i+1,j-1] + x[i-1,j-1])
        end
        d[1] -= a[1] * x[i, 1];   a[1] = 0.0
        d[n] -= c[n] * x[i, Nj];  c[n] = 0.0

        thomas_solve!(a, b, c, d, u)
        for k in 1:n
            j = k + 1
            x[i, j] = (1.0 - omega) * x[i, j] + omega * u[k]
        end

        # --- Solve for y-coordinates ---
        for k in 1:n
            j = k + 1
            a[k] = gamma[i,j] - J[i,j]^2 * Q[i,j] / 2.0
            b[k] = -2.0 * (alpha[i,j] + gamma[i,j])
            c[k] = gamma[i,j] + J[i,j]^2 * Q[i,j] / 2.0
            d[k] = -alpha[i,j] * (y[i+1,j] + y[i-1,j]) -
                    J[i,j]^2 * P[i,j] / 2.0 * (y[i+1,j] - y[i-1,j]) +
                    beta[i,j] / 2.0 * (y[i+1,j+1] - y[i-1,j+1] - y[i+1,j-1] + y[i-1,j-1])
        end
        d[1] -= a[1] * y[i, 1];   a[1] = 0.0
        d[n] -= c[n] * y[i, Nj];  c[n] = 0.0

        thomas_solve!(a, b, c, d, u)
        for k in 1:n
            j = k + 1
            y[i, j] = (1.0 - omega) * y[i, j] + omega * u[k]
        end
    end
end
```

### 6.3 Derivation of the Tridiagonal Coefficients

The discrete governing equation at interior point $(i,j)$ is (equation 7 from the
theory document):

$$\alpha\!\left(x_{i+1,j} - 2x_{i,j} + x_{i-1,j}\right) - \frac{\beta}{2}\!\left(x_{i+1,j+1} - x_{i-1,j+1} - x_{i+1,j-1} + x_{i-1,j-1}\right) + \gamma\!\left(x_{i,j+1} - 2x_{i,j} + x_{i,j-1}\right)$$

$$+ \;\frac{J^2 P}{2}\!\left(x_{i+1,j} - x_{i-1,j}\right) + \frac{J^2 Q}{2}\!\left(x_{i,j+1} - x_{i,j-1}\right) = 0$$

**For the $\xi$-line sweep** (unknowns at $j$-level: $x_{i-1,j}$, $x_{i,j}$, $x_{i+1,j}$):

Collecting terms by unknown:

| Unknown | Coefficient |
|---------|------------|
| $x_{i-1,j}$ | $\alpha - J^2 P/2 \;=\; a_i$ |
| $x_{i,j}$ | $-2(\alpha + \gamma) \;=\; b_i$ |
| $x_{i+1,j}$ | $\alpha + J^2 P/2 \;=\; c_i$ |

Everything else goes to $d_i$ (the RHS).

**For the $\eta$-line sweep** (unknowns at $i$-level: $x_{i,j-1}$, $x_{i,j}$, $x_{i,j+1}$):

| Unknown | Coefficient |
|---------|------------|
| $x_{i,j-1}$ | $\gamma - J^2 Q/2 \;=\; a_j$ |
| $x_{i,j}$ | $-2(\alpha + \gamma) \;=\; b_j$ |
| $x_{i,j+1}$ | $\gamma + J^2 Q/2 \;=\; c_j$ |

### 6.4 Treatment of Cross-Derivative Terms

The cross-derivative terms (involving $\beta$) reference **corner neighbors** like
$x_{i+1,j+1}$. These always use the **most recently computed values** (Gauss-Seidel
ordering):
- In the $\xi$-line sweep at $j$-level $j$: $j{-}1$ values are already updated
  (sweeping upward), $j{+}1$ values are from the previous iteration
- In the $\eta$-line sweep at $i$-level $i$: $i{-}1$ values are already updated
  (sweeping rightward), $i{+}1$ values are from the previous $\xi$-sweep

---

## 7. Main Solver Loop

### 7.1 Function: `StegerSorensonSolver`

```julia
using LinearAlgebra

function StegerSorensonSolver(x::Matrix{Float64}, y::Matrix{Float64}; params::SSParams)
    if params.skipBlock
        return x, y, 0.0, 0
    end

    Ni, Nj = size(x)
    @assert size(y) == (Ni, Nj) "x and y must have the same dimensions"
    @assert Ni >= 4 && Nj >= 4 "Grid must be at least 4x4 for wall forcing stencils"

    x = copy(x)   # Don't modify input
    y = copy(y)

    final_error = 0.0
    final_iter = 0

    for iter in 1:params.max_iter
        final_iter = iter

        # Save old grid for convergence check
        x_old = copy(x)
        y_old = copy(y)

        # --- Step 1: Compute metrics from current grid ---
        alpha, beta, gamma, J, x_xi, y_xi, x_eta, y_eta = compute_metrics(x, y)

        # --- Step 2: Compute control functions P, Q ---
        P, Q = compute_control_functions(
            x, y, alpha, beta, gamma, J,
            x_xi, y_xi, x_eta, y_eta;
            params=params
        )

        # --- Step 3: Xi-line sweep ---
        xi_line_sweep!(x, y, alpha, beta, gamma, J, P, Q, params.omega)

        # --- Step 4: Eta-line sweep (recompute metrics after xi-sweep) ---
        alpha, beta, gamma, J, _, _, _, _ = compute_metrics(x, y)
        eta_line_sweep!(x, y, alpha, beta, gamma, J, P, Q, params.omega)

        # --- Step 5: Convergence check (L-infinity norm) ---
        err_x = maximum(abs.(x .- x_old))
        err_y = maximum(abs.(y .- y_old))
        final_error = max(err_x, err_y)

        if params.verbose && (iter % params.print_interval == 0)
            println("  Iter $iter: error = $final_error")
        end

        if final_error < params.tol
            if params.verbose
                println("  Converged at iteration $iter, error = $final_error")
            end
            break
        end

        if iter == params.max_iter && params.verbose
            println("  WARNING: Max iterations ($iter) reached, error = $final_error")
        end
    end

    return x, y, final_error, final_iter
end
```

### 7.2 Design Notes

**Metric recomputation between sweeps** (Step 4): The metrics are recomputed after the
$\xi$-line sweep and before the $\eta$-line sweep. The $P$, $Q$ values are NOT
recomputed between sweeps (they change slowly and recomputing them is expensive).

**Minimum grid size**: The wall-forcing stencil requires three points into the domain
from each wall. Grids smaller than $4 \times 4$ cannot use wall forcing.

---

## 8. Wrapper for SmoothBlocks Integration

To integrate with the existing `SmoothBlocks` function in `src/smoothing/SmoothBlocks.jl`,
add a new solver branch:

```julia
function SmoothBlocks(blocks; solver=:stegerSorenson, params)
    smoothBlocks = Vector{Array{Float64,3}}(undef, length(blocks))
    finalErrors = Vector{Float64}(undef, length(blocks))
    finalIterations = Vector{Int}(undef, length(blocks))

    for i in eachindex(blocks)
        if params[i].skipBlock
            @info "Skipping block $i"
            smoothBlocks[i] = blocks[i]
            finalErrors[i] = 0.0
            finalIterations[i] = 0
            continue
        end

        xr, yr, err, iters = StegerSorensonSolver(
            blocks[i][1, :, :],
            blocks[i][2, :, :];
            params=params[i]
        )

        smoothBlocks[i] = permutedims(cat(xr, yr, dims=3), (3, 1, 2))
        finalErrors[i] = err
        finalIterations[i] = iters
    end

    @info "Smoothing convergence: $finalErrors"
    @info "Smoothing iterations: $finalIterations"
    return smoothBlocks, finalErrors, finalIterations
end
```

---

## 9. Complete File Structure

```
src/
  numerics/
    StegerSorensonSolver.jl    # Contains ALL of the following:
                                #   - compute_metrics()      (Section 3)
                                #   - compute_PQ_eta_wall()  (Section 4.2)
                                #   - compute_PQ_xi_wall()   (Section 4.3)
                                #   - compute_control_functions()  (Section 4.1)
                                #   - thomas_solve!()        (Section 5)
                                #   - xi_line_sweep!()       (Section 6.1)
                                #   - eta_line_sweep!()      (Section 6.2)
                                #   - StegerSorensonSolver() (Section 7)
  types.jl                     # Add SSParams struct (Section 2)
  smoothing/
    SmoothBlocks.jl            # Add :stegerSorenson branch (Section 8)
```

---

## 10. Testing Strategy

### 10.1 Unit Grid (Trivial Test)

A uniform rectangular grid should be unchanged by the smoother (it already satisfies
the Laplace equation exactly):

```julia
Ni, Nj = 20, 10
x = repeat(LinRange(0, 1, Ni), 1, Nj)
y = repeat(LinRange(0, 1, Nj)', Ni, 1)
params = SSParams(max_iter=100, tol=1e-10, omega=1.0, verbose=true)
x_out, y_out, err, iters = StegerSorensonSolver(x, y; params=params)
@assert err < 1e-10 "Uniform grid should converge in 1 iteration"
@assert maximum(abs.(x_out - x)) < 1e-12 "Grid should not move"
```

### 10.2 Perturbed Grid

Start with a uniform grid, add random perturbation to interior points, and verify the
smoother recovers something close to uniform:

```julia
x_pert = copy(x)
x_pert[2:Ni-1, 2:Nj-1] .+= 0.01 * randn(Ni-2, Nj-2)
params = SSParams(max_iter=5000, tol=1e-8, omega=1.2, verbose=true)
x_out, y_out, err, iters = StegerSorensonSolver(x_pert, y; params=params)
# Should converge and produce a smooth grid
```

### 10.3 Curved Domain

Test with a curved boundary (e.g., quarter annulus) generated by TFI:

```julia
params = SSParams(
    omega=1.2, tol=1e-6,
    useBottomWall=true, useTopWall=true,
    useLeftWall=true, useRightWall=true,
    decay_bottom=0.4, decay_top=0.4,
    decay_left=0.4, decay_right=0.4,
    verbose=true
)
x_out, y_out, err, iters = StegerSorensonSolver(x_tfi, y_tfi; params=params)
# Verify orthogonality at walls: beta should be near zero at boundaries
```

### 10.4 Orthogonality Verification

After smoothing with wall forcing, measure the angle between grid lines at each
boundary. The dot product of tangent and normal grid-line directions should be near zero:

```julia
# At bottom wall (j=1), check orthogonality:
for i in 2:Ni-1
    t = [x_out[i+1,1] - x_out[i-1,1], y_out[i+1,1] - y_out[i-1,1]]  # tangent
    n = [x_out[i,2] - x_out[i,1], y_out[i,2] - y_out[i,1]]            # normal
    angle = acosd(abs(dot(t, n) / (norm(t) * norm(n))))
    @assert abs(90.0 - angle) < 2.0 "Orthogonality error > 2 degrees at i=$i"
end
```

### 10.5 Convergence Rate Test

Compare convergence rates at different $\omega$ values to find the optimal relaxation:

```julia
for omega in [0.8, 1.0, 1.2, 1.4, 1.6]
    params = SSParams(omega=omega, max_iter=10000, tol=1e-8, verbose=false)
    _, _, err, iters = StegerSorensonSolver(x_test, y_test; params=params)
    println("omega=$omega: $iters iterations, error=$err")
end
```

---

## 11. Equation Quick Reference

All equations in one place for fast implementation reference.

**Metric coefficients** (compute at all grid points):

$$\alpha = x_\eta^2 + y_\eta^2, \qquad \beta = x_\xi\,x_\eta + y_\xi\,y_\eta, \qquad \gamma = x_\xi^2 + y_\xi^2, \qquad J = x_\xi\,y_\eta - x_\eta\,y_\xi$$

**Desired derivatives for orthogonality** (bottom/top wall):

$$(x_\eta)_d = \frac{-s\,y_\xi}{\sqrt{x_\xi^2 + y_\xi^2}}, \qquad (y_\eta)_d = \frac{s\,x_\xi}{\sqrt{x_\xi^2 + y_\xi^2}}$$

**Desired derivatives for orthogonality** (left/right wall):

$$(x_\xi)_d = \frac{-s\,y_\eta}{\sqrt{x_\eta^2 + y_\eta^2}}, \qquad (y_\xi)_d = \frac{s\,x_\eta}{\sqrt{x_\eta^2 + y_\eta^2}}$$

**One-sided second derivative with desired first derivative**:

$$f'' = \frac{1}{2}\!\left(-7\,f_{\text{wall}} + 8\,f_{\text{wall}+\text{dir}} - f_{\text{wall}+2\text{dir}}\right) + \sigma \cdot 3\,f'_d$$

where $\sigma = -\text{dir}$ for $\eta$-walls, $\sigma = +\text{dir}$ for $\xi$-walls.

**$2 \times 2$ system for $P$, $Q$ at boundary**:

$$J_b^2\!\left(P\,t_x + Q\,n_x\right) = -\!\left(\alpha_b\,f_{tt} + \gamma_b\,f_{nn}\right)$$

$$J_b^2\!\left(P\,t_y + Q\,n_y\right) = -\!\left(\alpha_b\,g_{tt} + \gamma_b\,g_{nn}\right)$$

**Exponential decay**:

$$P_{i,j} = \sum_{\text{walls}} P_{\text{wall}} \cdot e^{-c\,|\text{dist}|}, \qquad Q_{i,j} = \sum_{\text{walls}} Q_{\text{wall}} \cdot e^{-c\,|\text{dist}|}$$

**$\xi$-line tridiagonal coefficients**:

$$a_i = \alpha - \frac{J^2 P}{2}, \qquad b_i = -2(\alpha + \gamma), \qquad c_i = \alpha + \frac{J^2 P}{2}$$

$$d_i = -\gamma\!\left(x_{i,j+1}+x_{i,j-1}\right) - \frac{J^2 Q}{2}\!\left(x_{i,j+1}-x_{i,j-1}\right) + \frac{\beta}{2}\!\left(x_{i+1,j+1}-x_{i-1,j+1}-x_{i+1,j-1}+x_{i-1,j-1}\right)$$

**$\eta$-line tridiagonal coefficients**:

$$a_j = \gamma - \frac{J^2 Q}{2}, \qquad b_j = -2(\alpha + \gamma), \qquad c_j = \gamma + \frac{J^2 Q}{2}$$

$$d_j = -\alpha\!\left(x_{i+1,j}+x_{i-1,j}\right) - \frac{J^2 P}{2}\!\left(x_{i+1,j}-x_{i-1,j}\right) + \frac{\beta}{2}\!\left(x_{i+1,j+1}-x_{i-1,j+1}-x_{i+1,j-1}+x_{i-1,j-1}\right)$$

**SOR relaxation**:

$$x_{i,j} \leftarrow (1 - \omega)\,x_{i,j}^{\text{old}} + \omega\,x_{i,j}^{\text{Thomas}}$$

**Convergence**:

$$\varepsilon = \max\!\left(\max_{i,j}|x_{i,j}^{\text{new}} - x_{i,j}^{\text{old}}|,\;\max_{i,j}|y_{i,j}^{\text{new}} - y_{i,j}^{\text{old}}|\right) < \text{tol}$$
