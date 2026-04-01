# Mathematical Deep Dive: Elliptic Grid Smoothing (`ellipticSS`)

This document provides a rigorous mathematical derivation and analysis of the elliptic
smoothing method implemented in `GridGeneration.jl`. We trace the theory from the
continuous PDE formulation through discretization, the SOR iteration scheme, and the
wall-forcing orthogonality constraints.

---

## 1. The Fundamental Idea

Structured grid generation seeks a smooth mapping from a rectangular computational domain
to a curved physical domain:

```
(xi, eta) in [0,1] x [0,1]  -->  (x, y) in physical space
```

The initial grid (from Transfinite Interpolation) satisfies the boundary constraints but
may have poor interior quality -- skewed cells, non-smooth spacing, or non-orthogonal
grid lines. Elliptic smoothing corrects this by requiring the mapping to satisfy an
elliptic PDE system, which inherits the **maximum principle**: interior values are
governed by boundary values, guaranteeing smooth, fold-free grids.

---

## 2. Continuous Formulation

### 2.1 The Winslow Equations (Inverse Formulation)

The classical approach begins with the Laplace equations in physical space:

```
nabla^2(xi)  = P(xi, eta)
nabla^2(eta) = Q(xi, eta)
```

where P and Q are **control functions** (source terms). When P = Q = 0, these are the
Winslow equations. However, solving these in physical space is impractical because the
physical domain has a complex shape. We instead solve in computational space.

### 2.2 Transformation to Computational Space

Applying the chain rule and inverting the mapping, we obtain the **Poisson system in
computational coordinates**:

```
alpha * x_{xi,xi} - 2*beta * x_{xi,eta} + gamma * x_{eta,eta} + J^2 * (P * x_xi + Q * x_eta) = 0
alpha * y_{xi,xi} - 2*beta * y_{xi,eta} + gamma * y_{eta,eta} + J^2 * (P * y_xi + Q * y_eta) = 0
```

where the **metric coefficients** are:

```
alpha = x_eta^2 + y_eta^2       (measures eta-line stretching)
beta  = x_xi * x_eta + y_xi * y_eta   (measures non-orthogonality)
gamma = x_xi^2 + y_xi^2         (measures xi-line stretching)
J     = x_xi * y_eta - x_eta * y_xi   (Jacobian determinant)
```

**Key insight**: alpha, beta, gamma depend on the solution itself, making this a
**nonlinear** system even when P = Q = 0.

### 2.3 The Homogeneous Case (P = Q = 0)

With zero control functions, the system implemented in this code reduces to:

```
alpha * x_{xi,xi} - 2*beta * x_{xi,eta} + gamma * x_{eta,eta} = -S_x(xi, eta)
alpha * y_{xi,xi} - 2*beta * y_{xi,eta} + gamma * y_{eta,eta} = -S_y(xi, eta)
```

where S_x, S_y are the wall-forcing source terms (Section 5). Without wall forcing,
the right-hand side is zero and this is the **Laplace system** -- the simplest elliptic
smoother.

### 2.4 Why Elliptic?

The coefficient matrix of the second-order terms is:

```
A = [ alpha   -beta ]
    [ -beta   gamma ]
```

The eigenvalues of A are always positive because:

```
det(A) = alpha * gamma - beta^2
       = (x_eta^2 + y_eta^2)(x_xi^2 + y_xi^2) - (x_xi * x_eta + y_xi * y_eta)^2
       = (x_xi * y_eta - x_eta * y_xi)^2
       = J^2
```

By the Cauchy-Schwarz inequality, `alpha * gamma >= beta^2`, with equality only when the
Jacobian vanishes (degenerate mapping). As long as J != 0, the system is **strictly
elliptic**, ensuring smooth solutions with no interior extrema.

---

## 3. Finite Difference Discretization

### 3.1 Grid and Index Conventions

The computational grid has Ni points in the xi-direction and Nj points in the
eta-direction, with uniform spacing Delta_xi = Delta_eta = 1 (absorbed into coefficients):

```
x[i, j], y[i, j]     for i in 1:Ni, j in 1:Nj
```

Boundaries are fixed:
- Bottom: j = 1
- Top: j = Nj
- Left: i = 1
- Right: i = Ni

### 3.2 Derivative Approximations

**Interior points** (central differences, second-order accurate):

```
x_xi[i,j]   = (x[i+1,j] - x[i-1,j]) / 2
x_eta[i,j]  = (x[i,j+1] - x[i,j-1]) / 2
```

**Boundary points** (one-sided differences, first-order accurate):

```
x_xi[1,j]   = x[2,j] - x[1,j]          (forward)
x_xi[Ni,j]  = x[Ni,j] - x[Ni-1,j]      (backward)
```

(Analogous for eta-derivatives at j = 1 and j = Nj.)

**Second derivatives** (central differences):

```
x_{xi,xi}[i,j]   = x[i+1,j] - 2*x[i,j] + x[i-1,j]
x_{eta,eta}[i,j]  = x[i,j+1] - 2*x[i,j] + x[i,j-1]
```

**Cross derivative** (central differences on the 2D stencil):

```
x_{xi,eta}[i,j] = (x[i+1,j+1] - x[i-1,j+1] - x[i+1,j-1] + x[i-1,j-1]) / 4
```

### 3.3 Discrete Elliptic Equation

Substituting into the PDE for the x-coordinate at interior point (i,j):

```
alpha[i,j] * (x[i+1,j] - 2*x[i,j] + x[i-1,j])
 - 2*beta[i,j] * (x[i+1,j+1] - x[i-1,j+1] - x[i+1,j-1] + x[i-1,j-1]) / 4
 + gamma[i,j] * (x[i,j+1] - 2*x[i,j] + x[i,j-1])
 = RHS_x[i,j]
```

Solving for x[i,j]:

```
x[i,j] = [ alpha[i,j] * (x[i+1,j] + x[i-1,j])
          + gamma[i,j] * (x[i,j+1] + x[i,j-1])
          - beta[i,j]/2 * (x[i+1,j+1] - x[i-1,j+1] - x[i+1,j-1] + x[i-1,j-1])
          + RHS_x[i,j] ]
        / [ 2*(alpha[i,j] + gamma[i,j]) ]
```

This is exactly what appears in the code at lines 217-226 of `EllipticSolver.jl`.

### 3.4 The 9-Point Stencil

The discretization uses a **9-point computational stencil**:

```
        (i-1,j+1)   (i,j+1)   (i+1,j+1)
             \          |          /
              \         |         /
        (i-1,j) --- (i,j) --- (i+1,j)
              /         |         \
             /          |          \
        (i-1,j-1)   (i,j-1)   (i+1,j-1)
```

The weights on each neighbor are:

| Neighbor       | Weight in x-update                            |
|----------------|-----------------------------------------------|
| (i+1, j)       | +alpha[i,j]                                   |
| (i-1, j)       | +alpha[i,j]                                   |
| (i, j+1)       | +gamma[i,j]                                   |
| (i, j-1)       | +gamma[i,j]                                   |
| (i+1, j+1)     | -beta[i,j] / 2                                |
| (i-1, j+1)     | +beta[i,j] / 2                                |
| (i+1, j-1)     | +beta[i,j] / 2                                |
| (i-1, j-1)     | -beta[i,j] / 2                                |

The **diagonal weight** (coefficient of x[i,j]) is `2*(alpha + gamma)`.

When beta = 0 (orthogonal grid), the stencil reduces to the standard **5-point Laplacian**
with only cardinal neighbors.

---

## 4. SOR Iteration Scheme

### 4.1 Point Relaxation

The elliptic system is solved iteratively. At each iteration, every interior point is
updated using the latest available values (Gauss-Seidel ordering), then blended with the
old value using a relaxation parameter omega:

```
x*[i,j] = f(neighbors, RHS)     <-- Gauss-Seidel update
x^{n+1}[i,j] = (1 - omega) * x^n[i,j] + omega * x*[i,j]
```

### 4.2 Relaxation Parameter omega

The parameter omega controls convergence behavior:

| omega Range | Name               | Behavior                              |
|-------------|-------------------|---------------------------------------|
| 0 < omega < 1 | Under-relaxation  | Damps oscillations, very stable       |
| omega = 1     | Gauss-Seidel     | Standard iterative method             |
| 1 < omega < 2 | Over-relaxation   | Accelerates convergence (classic SOR) |

**Default value: omega = 0.2** (under-relaxation).

This is notably conservative. The reason: because the metric coefficients alpha, beta,
gamma are **recomputed from the solution at every iteration**, the system is nonlinear.
Large relaxation parameters can cause the nonlinear coupling to oscillate or diverge.
Under-relaxation stabilizes this feedback loop at the cost of slower convergence.

### 4.3 Why Not Over-Relax?

For a linear system A*x = b with known optimal omega, SOR theory (Young, 1954) gives:

```
omega_opt = 2 / (1 + sqrt(1 - rho(B)^2))
```

where rho(B) is the spectral radius of the Gauss-Seidel iteration matrix. For the
5-point Laplacian on an N x N grid, omega_opt ~ 2 - O(pi/N), which is close to 2.

However, this theory **does not apply here** because:
1. The coefficients alpha, beta, gamma change every iteration (nonlinear)
2. The forcing terms RHS are also recomputed every iteration
3. The mixed-derivative (beta) term breaks the standard convergence analysis

Under-relaxation (omega << 1) effectively linearizes the problem locally by taking small
steps, ensuring the nonlinear metric feedback remains stable.

### 4.4 Convergence Criterion

Convergence is measured by the L2 norm of the displacement in x-coordinates:

```
error = || x^{n+1} - x^n ||_2 = sqrt( sum_{i,j} (x^{n+1}[i,j] - x^n[i,j])^2 )
```

The iteration terminates when `error < tol` (default: 1e-6) or after `max_iter`
iterations (default: 5000).

**Note**: Only x-displacement is monitored, not y. This is a simplification; in practice,
both should converge together since they are driven by similar metric coefficients.

### 4.5 Iteration Complexity

Each SOR iteration requires:
- O(Ni * Nj) work to compute metrics (first derivatives at all points)
- O(Ni * Nj) work for the wall forcing (if enabled)
- O((Ni-2) * (Nj-2)) work for the interior sweep
- Total: **O(Ni * Nj) per iteration**

With under-relaxation (omega = 0.2), expect O(N^2) to O(N^3) iterations for an N x N
grid, giving total complexity of O(N^4) to O(N^5). This is why the solver can be slow
for large grids -- a multigrid approach would reduce this to O(N^2 log N).

---

## 5. Wall Forcing for Orthogonality

The most mathematically involved part of the implementation is the wall-forcing mechanism,
which computes source terms that drive grid lines to meet boundaries at right angles.

### 5.1 The Orthogonality Condition

At a boundary, we want the grid lines arriving at the wall to be **perpendicular** to
the wall. Mathematically, at a bottom wall (j = 1):

```
(x_xi, y_xi) . (x_eta, y_eta) = 0   at j = 1
```

where (x_xi, y_xi) is tangent to the wall and (x_eta, y_eta) is the direction of grid
lines leaving the wall.

### 5.2 Desired eta-Derivatives (Bottom/Top Walls)

Given the wall tangent vector (x_xi, y_xi) at each point along the boundary, the
orthogonal outward direction is:

```
n = (-y_xi, x_xi) / ||(x_xi, y_xi)||
```

We want the first grid line off the wall to be at a prescribed spacing s. So the
**desired** eta-derivatives at the wall are:

```
x_eta_desired = -s * y_xi / sqrt(x_xi^2 + y_xi^2)
y_eta_desired =  s * x_xi / sqrt(x_xi^2 + y_xi^2)
```

This appears in the code at lines 81-83 of `EllipticSolver.jl`.

The spacing s is linearly interpolated between the actual spacings at the two corners:

```
s1 = ||(x[1,1] - x[1,2], y[1,1] - y[1,2])||      (left corner spacing)
s2 = ||(x[Ni,1] - x[Ni,2], y[Ni,1] - y[Ni,2])||   (right corner spacing)
s(i) = s1 + (s2 - s1) * (i - 1) / (Ni - 1)
```

### 5.3 Computing Second Derivatives at the Wall

With the desired first derivatives known, we need the **second derivatives** at the wall
to compute forcing terms. The code uses a one-sided, second-order finite difference
formula. For the bottom wall (j = 1, dir = +1):

The standard one-sided second-order approximation for f'(0) given f(0), f(h), f(2h) is:

```
f'(0) = (-3*f(0) + 4*f(h) - f(2h)) / (2*h)
```

Rearranging and solving for f''(0), starting from the Taylor expansion:

```
f(h)  = f(0) + h*f'(0) + (h^2/2)*f''(0) + ...
f(2h) = f(0) + 2h*f'(0) + 2h^2*f''(0) + ...
```

Solving for f''(0):

```
f''(0) = f(2h) - 2*f(h) + f(0)   (standard second difference)
```

But we override f'(0) with the **desired** derivative rather than the actual one. The
code computes a modified second derivative that is consistent with the desired first
derivative:

```
x_etaeta = 0.5*(-7*x[i, wall] + 8*x[i, wall+1] - x[i, wall+2]) - 3*x_eta_desired
```

This formula arises from combining the one-sided second-derivative stencil with the
constraint that the first derivative takes the desired value. Specifically, if we define:

```
x_eta_actual ~ (-3*x[i,1] + 4*x[i,2] - x[i,3]) / 2
```

and the actual second derivative is:

```
x_etaeta_actual = x[i,3] - 2*x[i,2] + x[i,1]
```

then the modified second derivative that accounts for replacing x_eta_actual with
x_eta_desired involves a correction proportional to the difference. The factor of
-7, 8, -1 and the 3x_eta_desired term encode this correction.

### 5.4 Boundary Metric Coefficients

At the wall, the metric coefficients use the **desired** derivatives:

```
alpha_b = x_eta_desired^2 + y_eta_desired^2 = s^2
gamma_b = x_xi^2 + y_xi^2
```

Note that beta_b = 0 by construction (orthogonality enforced), so the mixed-derivative
term vanishes at the boundary.

### 5.5 The Forcing Source Terms

The RHS forcing at the boundary is:

```
RHS_x = -(alpha_b * x_{xi,xi} + gamma_b * x_{eta,eta})
RHS_y = -(alpha_b * y_{xi,xi} + gamma_b * y_{eta,eta})
```

This is the **residual** of the elliptic equation at the boundary with the desired
derivatives imposed. By feeding this residual as a source term, we drive the interior
solution toward a state that satisfies both the elliptic equation and the orthogonality
condition at the wall.

### 5.6 Exponential Decay into the Domain

The forcing is not applied only at the wall -- it is propagated into the domain with
exponential decay:

```
RHS_x_full[i, j] = RHS_x[i] * exp(-a_decay * |j - j_wall|)
RHS_y_full[i, j] = RHS_y[i] * exp(-b_decay * |j - j_wall|)
```

The decay parameters `a_decay` and `b_decay` control how far the wall influence extends:

| a_decay | Influence Region  | Effect                                     |
|---------|------------------|--------------------------------------------|
| ~0.1    | ~10 grid lines   | Deep penetration, strong orthogonality far from wall |
| ~0.4    | ~3 grid lines    | Moderate influence (default)               |
| ~0.9    | ~1 grid line     | Shallow forcing, orthogonal only at wall   |
| >> 1    | ~0 grid lines    | Essentially no forcing                     |

The 1/e decay distance is `1/a_decay` grid lines. For a_decay = 0.4, the forcing
drops to ~37% at 2.5 grid lines from the wall.

### 5.7 Left/Right Wall Forcing (xi-direction)

The same logic applies for left/right walls, but with the roles of xi and eta swapped.
At a left wall (i = 1), the desired xi-derivatives enforce orthogonality to eta-lines:

```
x_xi_desired = -s * y_eta / sqrt(x_eta^2 + y_eta^2)
y_xi_desired =  s * x_eta / sqrt(x_eta^2 + y_eta^2)
```

And the forcing decays in the i-direction:

```
RHS_x_full[i, j] = RHS_x_boundary[j] * exp(-a_decay * |i - i_wall|)
```

### 5.8 Superposition of Wall Forcing

When multiple walls have forcing enabled, the total source term is the **sum** of all
individual wall contributions:

```
RHS_total = RHS_bottom + RHS_top + RHS_left + RHS_right
```

This linear superposition means that in corners (where two walls meet), the forcing
from both walls combines. The exponential decay ensures that far from the walls, the
forcing is negligible and the interior behaves like a standard Laplace solve.

---

## 6. Complete Algorithm Summary

```
INPUT:  x[Ni, Nj], y[Ni, Nj]  (initial grid from TFI/edge solver)
        EllipticParams (omega, tol, max_iter, wall flags, decay rates)

FOR iter = 1 to max_iter:

    1. COMPUTE WALL FORCING (if enabled):
       For each enabled wall:
         a. Compute wall tangent vectors (x_xi or x_eta along wall)
         b. Compute perpendicular direction with desired spacing
         c. Compute modified second derivatives at wall
         d. Evaluate boundary metric coefficients
         e. Compute RHS residual at wall
         f. Propagate into domain with exponential decay
       Sum all wall contributions into RHS_x_full, RHS_y_full

    2. COMPUTE METRICS from current grid:
       For all (i,j):
         alpha[i,j] = x_eta^2 + y_eta^2
         beta[i,j]  = x_xi * x_eta + y_xi * y_eta
         gamma[i,j] = x_xi^2 + y_xi^2

    3. SOR SWEEP over interior (i=2:Ni-1, j=2:Nj-1):
       For each interior point:
         x* = [alpha*(x[i+1,j]+x[i-1,j]) + gamma*(x[i,j+1]+x[i,j-1])
               - beta/2*(x[i+1,j+1]-x[i-1,j+1]-x[i+1,j-1]+x[i-1,j-1])
               + RHS_x] / [2*(alpha + gamma)]
         x[i,j] = (1 - omega)*x[i,j] + omega*x*
         (same for y)

    4. CONVERGENCE CHECK:
       error = ||x^{n+1} - x^n||_2
       If error < tol: STOP

OUTPUT: smoothed x[Ni, Nj], y[Ni, Nj], final error, iteration count
```

---

## 7. Theoretical Properties

### 7.1 Existence and Uniqueness

For the continuous Laplace system with Dirichlet boundary conditions, existence and
uniqueness of a smooth solution follow from standard elliptic PDE theory (Gilbarg &
Trudinger). The discrete system inherits these properties: the coefficient matrix is
an M-matrix (diagonally dominant with positive diagonal, non-positive off-diagonal for
the 5-point part), guaranteeing a unique solution and convergence of iterative methods.

### 7.2 Grid Folding Prevention

The **maximum principle** for elliptic PDEs guarantees that interior grid coordinates
cannot exceed their boundary values. This means:
- Grid lines cannot cross (no folding)
- Cell volumes remain positive
- The Jacobian J stays nonzero (same sign as initial grid)

However, this guarantee is for the continuous problem. The discrete 9-point stencil
(with the beta cross-derivative term) can occasionally violate this if the grid is
extremely skewed. The under-relaxation parameter helps prevent such issues.

### 7.3 Smoothness

Elliptic solutions are infinitely differentiable in the interior (assuming smooth
boundaries). This translates to:
- Smooth variation of cell sizes
- Gradual changes in grid line angles
- No abrupt transitions in point spacing

These properties are precisely what makes elliptic smoothing desirable for computational
grids.

### 7.4 Relationship to Variational Principles

The Laplace equations for grid generation are the Euler-Lagrange equations of the
**Dirichlet energy functional**:

```
E[x, y] = integral integral [ (x_xi^2 + y_xi^2 + x_eta^2 + y_eta^2) ] d_xi d_eta
```

Minimizing E distributes grid points to minimize total stretching energy. This
variational interpretation explains why the smoothed grid tends toward uniform spacing
where boundaries allow it.

### 7.5 Effect of Control Functions

The wall forcing terms effectively add control functions P and Q that modify the
variational problem:

```
E_modified = E + integral integral [ P*x + Q*y ] d_xi d_eta
```

This allows biasing the minimization toward grids that also satisfy boundary
orthogonality, at the cost of slightly increased interior distortion (controlled by
the decay rate).

---

## 8. Parameter Sensitivity Guide

### 8.1 Relaxation Factor (omega)

```
omega = 0.05-0.1  : Very conservative. Use for extremely distorted initial grids.
omega = 0.15-0.3  : Standard range. Good balance of stability and speed. (default: 0.2)
omega = 0.4-0.6   : Aggressive. May work for mild smoothing tasks.
omega > 0.6       : Risk of divergence due to nonlinear metric coupling.
```

### 8.2 Decay Parameters (a_decay, b_decay)

The two decay parameters per wall (a and b) control decay of the x and y forcing
components independently. In practice they are often set equal.

```
decay ~ 0.1  : Strong influence, orthogonality enforced deep into domain
decay ~ 0.4  : Moderate influence (default), good general-purpose setting
decay ~ 1.0  : Weak influence, orthogonality only near boundary
decay ~ 2.0+ : Effectively no wall forcing
```

### 8.3 Convergence Tolerance (tol)

```
tol = 1e-4  : Rough smoothing, fast convergence
tol = 1e-6  : Standard (default), adequate for most applications
tol = 1e-8  : High precision, many more iterations required
tol = 1e-10 : Machine-precision level, rarely needed
```

### 8.4 Practical Tuning Strategy

1. Start with defaults (omega=0.2, decay=0.4, tol=1e-6)
2. If convergence is too slow: increase omega cautiously (try 0.3, then 0.4)
3. If diverging: decrease omega (try 0.1, then 0.05)
4. For boundary layer grids: use small decay (~0.1) on viscous wall, larger elsewhere
5. For inviscid grids: moderate decay (~0.4) on all walls
6. Monitor convergence history: should decrease monotonically after initial transient

---

## 9. Connection to the Broader Grid Generation Pipeline

```
Boundary Curves
      |
      v
  TFI (Transfinite Interpolation)  -- algebraic initial grid
      |
      v
  Block Splitting (optional)  -- multi-block decomposition
      |
      v
  Metric-Adaptive Edge Solver  -- redistribute points on edges using metric M
      |
      v
  Elliptic Smoothing (this document)  -- smooth interior using elliptic PDE
      |
      v
  Final Grid
```

The elliptic smoother receives a grid where boundaries are already well-resolved by
the edge solver. Its job is to:
1. Propagate the boundary point distribution smoothly into the interior
2. Remove any kinks or discontinuities from the TFI or block-splitting steps
3. Optionally enforce orthogonality at selected walls

**Note**: The current implementation does NOT incorporate the problem metric tensor M
(from `src/metric/`) into the elliptic smoothing. The smoothing uses only the geometric
metrics derived from the grid itself. Incorporating M would require additional source
terms in the elliptic equations to cluster points in regions of high metric variation.

---

## 10. References

1. **Thompson, Warsi, Mastin** (1985). *Numerical Grid Generation: Foundations and
   Applications*. The foundational text on elliptic grid generation.

2. **Winslow** (1967). "Numerical solution of the quasilinear Poisson equation in a
   nonuniform triangle mesh." *J. Comp. Phys.* 1(2), 149-172.

3. **Thomas & Middlecoff** (1980). "Direct control of the grid point distribution in
   meshes generated by elliptic equations." *AIAA Journal* 18(6), 652-656.
   (Control functions for boundary orthogonality.)

4. **Steger & Sorenson** (1979). "Automatic mesh-point clustering near a boundary in
   grid generation with elliptic partial differential equations." *J. Comp. Phys.* 33,
   405-410. (Wall forcing with exponential decay -- the approach used here.)

5. **Young** (1954). "Iterative methods for solving partial difference equations of
   elliptic type." *Trans. Amer. Math. Soc.* 76, 92-111. (SOR convergence theory.)

6. **Gilbarg & Trudinger** (2001). *Elliptic Partial Differential Equations of Second
   Order*. Springer. (Maximum principle and existence theory.)
