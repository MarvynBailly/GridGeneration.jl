# The Steger & Sorenson Elliptic Grid Smoothing Method: Mathematical Theory

A self-contained mathematical treatment of the elliptic grid generation method of
Steger & Sorenson (1979), as implemented in the GRAPE code (Sorenson, NASA TM-81198,
1980). This document derives every equation from first principles.

---

## 1. Problem Statement

We seek a smooth, invertible mapping from a rectangular computational domain to a
curved physical domain:

```
    Computational (xi, eta)              Physical (x, y)
    +-------------------+        F       +-----------------+
    |                   |     ------>   /                   \
    |   uniform grid    |              |    curved domain    |
    |                   |               \                   /
    +-------------------+                +-----------------+
```

Given: boundary curves defining the physical domain.
Find: interior point positions $x(\xi, \eta)$ and $y(\xi, \eta)$ such that the grid is
smooth, non-folding, and (optionally) orthogonal at boundaries.

---

## 2. Why Elliptic PDEs?

The key idea is to require the mapping to satisfy an elliptic PDE. Elliptic equations
have three critical properties for grid generation:

1. **Maximum principle**: Interior values are bounded by boundary values. This prevents
   grid lines from crossing (no folding) and keeps the Jacobian nonzero.

2. **Smoothness**: Solutions are infinitely differentiable in the interior (given smooth
   boundaries). This produces smoothly varying cell sizes and angles.

3. **Boundary control**: The solution throughout the domain is completely determined by
   boundary conditions. We fix boundary point positions (Dirichlet) and optionally
   prescribe how grid lines meet boundaries (Neumann for orthogonality).

---

## 3. The Governing Equations

### 3.1 Starting Point: Laplace Equations in Physical Space

The simplest elliptic system requires the computational coordinates to be harmonic
functions of the physical coordinates:

$$\frac{\partial^2 \xi}{\partial x^2} + \frac{\partial^2 \xi}{\partial y^2} = P(\xi, \eta) \tag{1a}$$

$$\frac{\partial^2 \eta}{\partial x^2} + \frac{\partial^2 \eta}{\partial y^2} = Q(\xi, \eta) \tag{1b}$$

where $P$ and $Q$ are control (source) functions. When $P = Q = 0$, these are the
Winslow equations (Winslow, 1967). The control functions allow clustering and
orthogonality control.

### 3.2 Transformation to Computational Coordinates

Equations (1) are defined on the physical domain, which has a complex shape. To solve
them on the uniform rectangular computational domain, we apply the chain rule to
transform all derivatives.

The mapping between coordinate systems is characterized by the Jacobian matrix:

$$\mathbf{J} = \begin{bmatrix} x_\xi & x_\eta \\ y_\xi & y_\eta \end{bmatrix}, \qquad J = \det(\mathbf{J}) = x_\xi \, y_\eta - x_\eta \, y_\xi$$

where subscripts denote partial derivatives (e.g., $x_\xi = \partial x / \partial \xi$).

The inverse transformation gives:

$$\xi_x = \frac{y_\eta}{J}, \quad \xi_y = \frac{-x_\eta}{J}, \quad \eta_x = \frac{-y_\xi}{J}, \quad \eta_y = \frac{x_\xi}{J}$$

Substituting into (1a) and expanding:

$$\xi_{xx} + \xi_{yy} = \frac{1}{J^2}\left[\alpha\,\xi_{\xi\xi} - 2\beta\,\xi_{\xi\eta} + \gamma\,\xi_{\eta\eta}\right] = P$$

But $\xi_{\xi\xi} = 0$, $\xi_{\xi\eta} = 0$, $\xi_{\eta\eta} = 0$ (since $\xi$ is
itself a coordinate). So this approach directly gives us nothing useful for $\xi$.
Instead, we use the **inverse formulation**: we solve for $x(\xi,\eta)$ and $y(\xi,\eta)$.

Starting from the identity $\nabla^2\xi = P$ in physical space and transforming entirely
to computational space, after lengthy algebra (see Thompson, Warsi & Mastin, Ch. 5), we
obtain:

$$\alpha\,x_{\xi\xi} - 2\beta\,x_{\xi\eta} + \gamma\,x_{\eta\eta} + J^2\left(P\,x_\xi + Q\,x_\eta\right) = 0 \tag{2a}$$

$$\alpha\,y_{\xi\xi} - 2\beta\,y_{\xi\eta} + \gamma\,y_{\eta\eta} + J^2\left(P\,y_\xi + Q\,y_\eta\right) = 0 \tag{2b}$$

### 3.3 Metric Coefficients

The coefficients in (2) are the components of the **covariant metric tensor** of the
mapping:

$$\alpha = x_\eta^2 + y_\eta^2 \qquad \text{(magnitude squared of the } \eta\text{-tangent vector)}$$

$$\beta = x_\xi\,x_\eta + y_\xi\,y_\eta \qquad \text{(dot product of } \xi\text{- and } \eta\text{-tangent vectors)}$$

$$\gamma = x_\xi^2 + y_\xi^2 \qquad \text{(magnitude squared of the } \xi\text{-tangent vector)}$$

$$J = x_\xi\,y_\eta - x_\eta\,y_\xi \qquad \text{(Jacobian determinant)}$$

**Geometric interpretation**:
- $\sqrt{\alpha}$ = length of the $\eta$-direction tangent vector (local $\eta$-line stretching)
- $\sqrt{\gamma}$ = length of the $\xi$-direction tangent vector (local $\xi$-line stretching)
- $\beta = 0$ when grid lines are locally orthogonal (tangent vectors perpendicular)
- $J > 0$ when the mapping preserves orientation (right-handed grid)

### 3.4 Ellipticity Proof

The coefficient matrix of second-order terms is:

$$\mathbf{A} = \begin{bmatrix} \alpha & -\beta \\ -\beta & \gamma \end{bmatrix}$$

Its determinant is:

$$\det(\mathbf{A}) = \alpha\gamma - \beta^2$$

Expanding:

$$\alpha\gamma = \left(x_\eta^2 + y_\eta^2\right)\left(x_\xi^2 + y_\xi^2\right)$$

$$\beta^2 = \left(x_\xi\,x_\eta + y_\xi\,y_\eta\right)^2$$

By the Cauchy-Schwarz inequality applied to vectors $(x_\xi, y_\xi)$ and $(x_\eta, y_\eta)$:

$$\left(x_\xi\,x_\eta + y_\xi\,y_\eta\right)^2 \leq \left(x_\xi^2 + y_\xi^2\right)\left(x_\eta^2 + y_\eta^2\right)$$

with equality iff the vectors are parallel. Therefore $\det(\mathbf{A}) \geq 0$.

In fact, the **Lagrange identity** gives us exactly:

$$\alpha\gamma - \beta^2 = \left(x_\xi\,y_\eta - x_\eta\,y_\xi\right)^2 = J^2$$

So $\det(\mathbf{A}) = J^2 > 0$ whenever the mapping is non-degenerate, proving the
system is **strictly elliptic**.

### 3.5 Nonlinearity

Equations (2) are **nonlinear** because $\alpha$, $\beta$, $\gamma$, and $J$ are
functions of the first derivatives of $x$ and $y$ -- the unknowns. This means:
- No closed-form solution exists in general
- Iterative methods are required
- Convergence is not guaranteed without care (relaxation parameters)

---

## 4. Control Functions for Boundary Orthogonality

### 4.1 The Steger & Sorenson Innovation

Prior to Steger & Sorenson, the control functions $P$ and $Q$ had to be specified by
trial and error. Their key contribution was an **automatic algorithm** to compute $P$
and $Q$ from boundary conditions, specifically to enforce:
1. Grid orthogonality at boundaries
2. Prescribed wall-normal spacing

### 4.2 Orthogonality Condition

At a boundary, orthogonality means the grid line arriving at the wall is perpendicular
to the wall tangent. At the bottom wall ($j = 1$), the wall tangent is
$(x_\xi,\, y_\xi)$ and the departing grid line direction is $(x_\eta,\, y_\eta)$.
Perpendicularity requires:

$$x_\xi\,x_\eta + y_\xi\,y_\eta = 0 \qquad (\text{i.e., } \beta = 0 \text{ at the wall}) \tag{3}$$

### 4.3 Desired Normal Derivatives

The unit outward normal to the wall at the bottom boundary is:

$$\hat{\mathbf{n}} = \frac{(-y_\xi,\; x_\xi)}{\sqrt{x_\xi^2 + y_\xi^2}}$$

(choosing the sign so $\hat{\mathbf{n}}$ points into the domain, i.e., in the $+\eta$
direction).

For a prescribed wall-normal spacing $s$ (the distance from the wall to the first
interior grid line), the desired $\eta$-derivatives are:

$$(x_\eta)_d = -\frac{s\,y_\xi}{\sqrt{x_\xi^2 + y_\xi^2}} \tag{4a}$$

$$(y_\eta)_d = \frac{s\,x_\xi}{\sqrt{x_\xi^2 + y_\xi^2}} \tag{4b}$$

**Verification of orthogonality** (condition 3):

$$x_\xi\,(x_\eta)_d + y_\xi\,(y_\eta)_d = x_\xi\!\left(\frac{-s\,y_\xi}{\|\mathbf{t}\|}\right) + y_\xi\!\left(\frac{s\,x_\xi}{\|\mathbf{t}\|}\right) = \frac{s}{\|\mathbf{t}\|}\left(-x_\xi\,y_\xi + y_\xi\,x_\xi\right) = 0 \;\checkmark$$

**Verification of spacing**:

$$\sqrt{(x_\eta)_d^2 + (y_\eta)_d^2} = \frac{s\sqrt{y_\xi^2 + x_\xi^2}}{\sqrt{x_\xi^2 + y_\xi^2}} = s \;\checkmark$$

### 4.4 Computing P and Q at the Boundary

At the bottom wall ($j = 1$), we know:
- The point positions $x_{i,1}$, $y_{i,1}$ (Dirichlet boundary condition)
- The desired first derivatives $(x_\eta)_d$, $(y_\eta)_d$ (from orthogonality)
- The tangential derivatives $x_\xi$, $y_\xi$, $x_{\xi\xi}$, $y_{\xi\xi}$ (from boundary geometry)
- The second normal derivatives $x_{\eta\eta}$, $y_{\eta\eta}$ (computed from a one-sided stencil, see Section 5.5)

The metric coefficients evaluated with the desired derivatives are:

$$\alpha_b = (x_\eta)_d^2 + (y_\eta)_d^2 = s^2$$

$$\beta_b = 0 \qquad \text{(by construction)}$$

$$\gamma_b = x_\xi^2 + y_\xi^2$$

$$J_b = x_\xi\,(y_\eta)_d - (x_\eta)_d\,y_\xi = \frac{s\left(x_\xi^2 + y_\xi^2\right)}{\|\mathbf{t}\|} = s\,\|\mathbf{t}\|$$

Substituting into equations (2a) and (2b) at the boundary:

$$\alpha_b\,x_{\xi\xi} + \gamma_b\,x_{\eta\eta} + J_b^2\!\left(P\,x_\xi + Q\,(x_\eta)_d\right) = 0$$

$$\alpha_b\,y_{\xi\xi} + \gamma_b\,y_{\eta\eta} + J_b^2\!\left(P\,y_\xi + Q\,(y_\eta)_d\right) = 0$$

This is a $2 \times 2$ linear system for $P$ and $Q$ at each boundary point:

$$\begin{bmatrix} x_\xi & (x_\eta)_d \\ y_\xi & (y_\eta)_d \end{bmatrix} \begin{bmatrix} P \\ Q \end{bmatrix} = \frac{-1}{J_b^2} \begin{bmatrix} \alpha_b\,x_{\xi\xi} + \gamma_b\,x_{\eta\eta} \\ \alpha_b\,y_{\xi\xi} + \gamma_b\,y_{\eta\eta} \end{bmatrix}$$

Since the coefficient matrix has determinant $= J_b \neq 0$ for a valid grid, this
system always has a unique solution.

### 4.5 Propagation into the Interior

The control functions $P$ and $Q$ are computed at the boundary by the above procedure.
They must be extended into the interior. Steger & Sorenson use **exponential decay**:

$$P(\xi,\,j) = P_{\text{boundary}}(\xi)\;\exp\!\left(-c\,|j - j_{\text{wall}}|\right) \tag{5a}$$

$$Q(\xi,\,j) = Q_{\text{boundary}}(\xi)\;\exp\!\left(-c\,|j - j_{\text{wall}}|\right) \tag{5b}$$

where $c$ is a decay rate parameter (larger $c$ = faster decay = less interior influence).

The $1/e$ influence depth is $1/c$ grid lines. Typical values:
- $c \sim 0.1$: influence extends $\sim 10$ grid lines from wall
- $c \sim 0.4$: influence extends $\sim 2{-}3$ grid lines (common default)
- $c \sim 1.0$: influence extends $\sim 1$ grid line
- $c \gg 1$: effectively no wall forcing

When multiple walls have forcing, the total $P$ and $Q$ are the **sum** of contributions
from all enabled walls. Since the exponential decays are localized near each wall, the
contributions are approximately disjoint except near corners.

### 4.6 Left/Right Walls

For left/right walls ($i = \text{const}$ boundaries), the roles of $\xi$ and $\eta$ swap.
The wall tangent direction is $(x_\eta,\,y_\eta)$, and the desired derivatives are in $\xi$:

$$(x_\xi)_d = \frac{-s\,y_\eta}{\sqrt{x_\eta^2 + y_\eta^2}} \tag{6a}$$

$$(y_\xi)_d = \frac{s\,x_\eta}{\sqrt{x_\eta^2 + y_\eta^2}} \tag{6b}$$

The $P$, $Q$ computation at the wall follows the same logic, and the decay is in the
$i$-direction:

$$P(i,\,\eta) = P_{\text{boundary}}(\eta)\;\exp\!\left(-c\,|i - i_{\text{wall}}|\right)$$

$$Q(i,\,\eta) = Q_{\text{boundary}}(\eta)\;\exp\!\left(-c\,|i - i_{\text{wall}}|\right)$$

---

## 5. Discretization

### 5.1 Grid Conventions

The computational domain is discretized with $N_i$ points in the $\xi$-direction and
$N_j$ points in the $\eta$-direction, with unit spacing
($\Delta\xi = \Delta\eta = 1$):

$$x_{i,j},\; y_{i,j} \qquad i \in \{1, \ldots, N_i\},\;\; j \in \{1, \ldots, N_j\}$$

Boundary identification:
- Bottom: $j = 1$
- Top: $j = N_j$
- Left: $i = 1$
- Right: $i = N_i$

### 5.2 First Derivatives (Central Differences, $O(h^2)$)

Interior ($2 \leq i \leq N_i{-}1$ or $2 \leq j \leq N_j{-}1$):

$$(x_\xi)_{i,j} = \frac{x_{i+1,j} - x_{i-1,j}}{2}, \qquad (x_\eta)_{i,j} = \frac{x_{i,j+1} - x_{i,j-1}}{2}$$

Boundaries (one-sided, $O(h)$):

$$(x_\xi)_{1,j} = x_{2,j} - x_{1,j}, \qquad (x_\xi)_{N_i,j} = x_{N_i,j} - x_{N_i-1,j}$$

$$(x_\eta)_{i,1} = x_{i,2} - x_{i,1}, \qquad (x_\eta)_{i,N_j} = x_{i,N_j} - x_{i,N_j-1}$$

### 5.3 Second Derivatives (Central Differences, $O(h^2)$)

$$(x_{\xi\xi})_{i,j} = x_{i+1,j} - 2\,x_{i,j} + x_{i-1,j}$$

$$(x_{\eta\eta})_{i,j} = x_{i,j+1} - 2\,x_{i,j} + x_{i,j-1}$$

$$(x_{\xi\eta})_{i,j} = \frac{x_{i+1,j+1} - x_{i-1,j+1} - x_{i+1,j-1} + x_{i-1,j-1}}{4}$$

### 5.4 Discrete Governing Equation

Substituting into (2a) at an interior point $(i, j)$:

$$\alpha_{i,j}\!\left(x_{i+1,j} - 2\,x_{i,j} + x_{i-1,j}\right) - \frac{\beta_{i,j}}{2}\!\left(x_{i+1,j+1} - x_{i-1,j+1} - x_{i+1,j-1} + x_{i-1,j-1}\right)$$

$$+ \;\gamma_{i,j}\!\left(x_{i,j+1} - 2\,x_{i,j} + x_{i,j-1}\right) + J_{i,j}^2\!\left(P_{i,j}\,\frac{x_{i+1,j} - x_{i-1,j}}{2} + Q_{i,j}\,\frac{x_{i,j+1} - x_{i,j-1}}{2}\right) = 0 \tag{7}$$

### 5.5 One-Sided Second Derivatives at the Wall

To compute $P$ and $Q$ at the boundary (Section 4.4), we need $x_{\eta\eta}$ at $j = 1$
where only forward points are available. We use the known desired first derivative to
improve accuracy.

Taylor expansions with unit spacing ($h = 1$) at the bottom wall ($j = 1$):

$$x_{i,2} = x_{i,1} + x_\eta + \tfrac{1}{2}\,x_{\eta\eta} + \tfrac{1}{6}\,x_{\eta\eta\eta} + \cdots$$

$$x_{i,3} = x_{i,1} + 2\,x_\eta + 2\,x_{\eta\eta} + \tfrac{4}{3}\,x_{\eta\eta\eta} + \cdots$$

We have two unknowns ($x_{\eta\eta}$ and $x_{\eta\eta\eta}$) and we **replace** $x_\eta$
with the desired value $(x_\eta)_d$.

**Derivation**: Let $f_0 = x_{i,1}$, $f_1 = x_{i,2}$, $f_2 = x_{i,3}$,
$f' = (x_\eta)_d$, $f'' =$ unknown.

From the Taylor series:

$$f_1 = f_0 + f' + \frac{f''}{2} + \frac{f'''}{6} \tag{A}$$

$$f_2 = f_0 + 2f' + 2f'' + \frac{4f'''}{3} \tag{B}$$

From (A): $\;f''' = 6\!\left(f_1 - f_0 - f' - f''/2\right)$

Substituting into (B):

$$f_2 = f_0 + 2f' + 2f'' + 8\!\left(f_1 - f_0 - f' - f''/2\right) = -7f_0 + 8f_1 - 6f' - 2f''$$

Solving for $f''$:

$$\boxed{f'' = \frac{1}{2}\!\left(-7f_0 + 8f_1 - f_2\right) - 3f'} \tag{8}$$

This uses three grid points ($j{=}1, j{=}2, j{=}3$) plus the desired first derivative
to compute a second-order accurate second derivative at the wall.

For the **top wall** ($j = N_j$, $\text{dir} = -1$):

$$x_{\eta\eta} = \frac{1}{2}\!\left(-7\,x_{i,N_j} + 8\,x_{i,N_j-1} - x_{i,N_j-2}\right) + 3\,(x_\eta)_d$$

The sign of the $3(x_\eta)_d$ term flips because the one-sided stencil points in the
opposite direction.

For **left/right walls**, the same formulas apply with $\xi$ replacing $\eta$ and $i$
replacing $j$.

---

## 6. Line Gauss-Seidel with Thomas Algorithm

This is the core solver and the key algorithmic feature of the GRAPE implementation.
Instead of updating each grid point independently (point SOR), the line method solves
for all points along a grid line simultaneously.

### 6.1 Forming the Tridiagonal System

We sweep through $j$-levels from $j = 2$ to $j = N_j{-}1$. At each $j$-level, we
rearrange equation (7) so that terms involving $x_{i-1,j}$, $x_{i,j}$, and $x_{i+1,j}$
(all at the same $j$-level) form the tridiagonal unknowns, and everything else is
treated as known.

Collecting the $\xi$-direction implicit terms from the second-derivative and
first-derivative discretizations:

$$a_i = \alpha_{i,j} - \frac{J_{i,j}^2\,P_{i,j}}{2} \qquad \text{(sub-diagonal, coefficient of } x_{i-1,j}\text{)}$$

$$b_i = -2\!\left(\alpha_{i,j} + \gamma_{i,j}\right) \qquad \text{(diagonal, coefficient of } x_{i,j}\text{)}$$

$$c_i = \alpha_{i,j} + \frac{J_{i,j}^2\,P_{i,j}}{2} \qquad \text{(super-diagonal, coefficient of } x_{i+1,j}\text{)}$$

The right-hand side (everything at $j{-}1$, $j{+}1$, and the cross-derivative corners):

$$d_i = -\gamma_{i,j}\!\left(x_{i,j+1} + x_{i,j-1}\right) - \frac{J_{i,j}^2\,Q_{i,j}}{2}\!\left(x_{i,j+1} - x_{i,j-1}\right) + \frac{\beta_{i,j}}{2}\!\left(x_{i+1,j+1} - x_{i-1,j+1} - x_{i+1,j-1} + x_{i-1,j-1}\right)$$

Note: the cross-derivative terms in $d_i$ use $x$ values at $(i{\pm}1,\, j{\pm}1)$.
These use the **most recently computed values** (Gauss-Seidel ordering).

Boundary conditions at $i = 1$ and $i = N_i$: these are Dirichlet (fixed), so
$x_{1,j}$ and $x_{N_i,j}$ are known. The tridiagonal system is solved for
$i = 2, \ldots, N_i{-}1$.

For $i = 2$: the $a_2 \cdot x_{1,j}$ term is moved to the right-hand side.
For $i = N_i{-}1$: the $c_{N_i-1} \cdot x_{N_i,j}$ term is moved to the right-hand side.

### 6.2 The Thomas Algorithm

The Thomas algorithm (TDMA) solves a tridiagonal system in $O(n)$ operations:

$$a_i\,u_{i-1} + b_i\,u_i + c_i\,u_{i+1} = d_i \qquad \text{for } i = 2, \ldots, N_i{-}1$$

**Forward sweep** (eliminate sub-diagonal):

$$\text{For } i = 3, 4, \ldots, N_i{-}1: \qquad w = \frac{a_i}{b_{i-1}}, \quad b_i \leftarrow b_i - w\,c_{i-1}, \quad d_i \leftarrow d_i - w\,d_{i-1}$$

**Back substitution**:

$$u_{N_i-1} = \frac{d_{N_i-1}}{b_{N_i-1}}, \qquad \text{For } i = N_i{-}2, \ldots, 2: \quad u_i = \frac{d_i - c_i\,u_{i+1}}{b_i}$$

This is exact (not iterative) and costs $O(N_i)$ operations per $j$-level. The total
cost for one sweep over all $j$-levels is $O(N_i \cdot N_j)$ -- same as point SOR per
iteration, but with far superior convergence properties.

### 6.3 Alternating Direction Sweeps

For even faster convergence, GRAPE alternates:

**Pass 1**: Sweep $j = 2, \ldots, N_j{-}1$ solving tridiagonal systems along $\xi$-lines.

**Pass 2**: Sweep $i = 2, \ldots, N_i{-}1$ solving tridiagonal systems along $\eta$-lines.

In the $\eta$-line sweep, the implicit unknowns are at $(i,\,j{-}1)$, $(i,\,j)$,
$(i,\,j{+}1)$, and the known terms are at $(i{-}1,\,j)$ and $(i{+}1,\,j)$.
The coefficients become:

$$a_j = \gamma_{i,j} - \frac{J_{i,j}^2\,Q_{i,j}}{2}, \qquad b_j = -2\!\left(\alpha_{i,j} + \gamma_{i,j}\right), \qquad c_j = \gamma_{i,j} + \frac{J_{i,j}^2\,Q_{i,j}}{2}$$

$$d_j = -\alpha_{i,j}\!\left(x_{i+1,j} + x_{i-1,j}\right) - \frac{J_{i,j}^2\,P_{i,j}}{2}\!\left(x_{i+1,j} - x_{i-1,j}\right) + \frac{\beta_{i,j}}{2}\!\left(x_{i+1,j+1} - x_{i-1,j+1} - x_{i+1,j-1} + x_{i-1,j-1}\right)$$

Alternating between $\xi$-line and $\eta$-line sweeps provides implicit coupling in both
directions, dramatically accelerating convergence.

### 6.4 SOR Acceleration

After solving the tridiagonal system to get $x_{i,j}^{\text{new}}$, the SOR relaxation
is applied:

$$x_{i,j} = (1 - \omega)\,x_{i,j}^{\text{old}} + \omega\,x_{i,j}^{\text{new}}$$

Because the line solver provides implicit coupling (much more stable than point updates),
**over-relaxation is viable**:

$$\omega \in [1.0,\; 1.8] \qquad \text{(typical range for line Gauss-Seidel)}$$

$$\omega \approx 1.4 \text{--} 1.6 \qquad \text{(commonly optimal for grid generation)}$$

Compare to point SOR which typically requires under-relaxation ($\omega \sim 0.1{-}0.4$)
for this nonlinear problem.

---

## 7. Convergence Properties

### 7.1 Convergence Criterion

Monitor both $x$ and $y$ displacements:

$$\varepsilon_x = \max_{i,j}\left|x_{i,j}^{n+1} - x_{i,j}^{n}\right|, \qquad \varepsilon_y = \max_{i,j}\left|y_{i,j}^{n+1} - y_{i,j}^{n}\right|$$

$$\varepsilon = \max\!\left(\varepsilon_x,\; \varepsilon_y\right)$$

The $L^\infty$ norm is preferred over $L^2$ because it catches localized
non-convergence. Convergence is declared when $\varepsilon < \text{tol}$ (typically
$10^{-6}$ to $10^{-8}$).

### 7.2 Convergence Rates

| Solver | Iterations for $N \times N$ grid | Total work |
|--------|--------------------------------|------------|
| Point Gauss-Seidel ($\omega=1$) | $O(N^2)$ | $O(N^4)$ |
| Point SOR (optimal $\omega$) | $O(N)$ | $O(N^3)$ |
| Point SOR (under-relaxed, $\omega \sim 0.2$) | $O(N^2){-}O(N^3)$ | $O(N^4){-}O(N^5)$ |
| **Line Gauss-Seidel** | $O(N)$ | $O(N^3)$ |
| **Line SOR ($\omega \sim 1.5$)** | $O(\sqrt{N})$ | $O(N^2\sqrt{N})$ |
| Multigrid | $O(1)$ | $O(N^2)$ |

The line Gauss-Seidel with SOR achieves the best practical convergence without the
implementation complexity of multigrid.

### 7.3 Nonlinear Convergence

Because the metric coefficients $\alpha$, $\beta$, $\gamma$ change every iteration,
convergence is not guaranteed by linear SOR theory. The Picard iteration approach (freeze
coefficients, solve, update coefficients, repeat) converges when:

1. The relaxation parameter $\omega$ is not too aggressive
2. The initial grid is reasonable (not severely folded)
3. The forcing terms don't create conflicting constraints

In practice, line Gauss-Seidel with $\omega \sim 1.0{-}1.4$ converges reliably for all
reasonable grid generation problems.

---

## 8. Wall Spacing Prescription

### 8.1 The Spacing Function $s(\xi)$

Along each boundary, the wall-normal spacing $s$ can be specified as any positive
function of the tangential coordinate. Common choices:

**Uniform spacing**: $s(\xi) = s_0$ (constant along the wall)

**Geometric stretching**: $s$ varies to match the tangential spacing:

$$s(\xi) = r \cdot \Delta\xi_{\text{tangential}}(\xi)$$

where $r$ is an aspect ratio parameter.

**Interpolated from corners**: The simplest approach (used in the current
GridGeneration.jl implementation) linearly interpolates between the actual spacings at
the two endpoints:

$$s(\xi) = s_{\text{left}} + \left(s_{\text{right}} - s_{\text{left}}\right)\frac{\xi - \xi_{\text{left}}}{\xi_{\text{right}} - \xi_{\text{left}}}$$

**User-prescribed**: GRAPE allows the user to specify $s$ at each boundary point.

### 8.2 Spacing and the Desired Derivatives

The spacing $s$ enters the desired derivatives (equations 4a, 4b) and thereby controls:
- The distance to the first interior grid line
- The local cell aspect ratio at the wall
- Indirectly (through $P$, $Q$ propagation) the spacing several lines into the interior

---

## 9. Variational Interpretation

### 9.1 Dirichlet Energy

The Laplace equations ($P = Q = 0$) are the Euler-Lagrange equations of the Dirichlet
energy functional:

$$E[x, y] = \frac{1}{2}\iint\!\left(x_\xi^2 + y_\xi^2 + x_\eta^2 + y_\eta^2\right)\,d\xi\,d\eta$$

Minimizing $E$ produces a mapping that minimizes total stretching. In a sense, the grid
"relaxes" to the smoothest possible configuration consistent with boundary constraints.

### 9.2 With Control Functions

Adding $P$ and $Q$ modifies the variational problem:

$$E_{\text{modified}} = E + \iint\!\left(P\,x + Q\,y\right)\,d\xi\,d\eta$$

The control functions bias the energy minimization toward configurations that also
satisfy the orthogonality constraints. The exponential decay ensures the bias is
strongest near walls and fades in the interior.

### 9.3 Physical Analogy

The Laplace-smoothed grid can be thought of as a membrane under tension: the boundary
is pinned, and the interior settles to the minimum-energy (smoothest) shape. The
control functions are like point loads applied near the boundary to tilt the grid lines
toward the desired angle.

---

## 10. Summary of the Complete Steger & Sorenson Method

```
INPUT:  Initial grid x[Ni,Nj], y[Ni,Nj]
        Boundary specifications: which walls have orthogonality forcing
        Wall spacing functions: s(xi) or s(eta) for each enabled wall
        Parameters: omega, decay rates, tolerance, max_iterations

PREPROCESSING:
    Compute P, Q = 0 everywhere (will be overwritten at boundaries)

FOR iter = 1 to max_iterations:

    1. RECOMPUTE CONTROL FUNCTIONS at each enabled boundary:
       - Compute tangent vectors along the wall
       - Compute perpendicular direction with desired spacing (eqs. 4)
       - Compute second derivatives at wall (eq. 8)
       - Solve 2x2 system for P, Q at each boundary point (Section 4.4)
       - Propagate P, Q into interior with exponential decay (eqs. 5)
       - Superpose contributions from all enabled walls

    2. COMPUTE METRICS from current grid:
       alpha, beta, gamma, J at every point (Section 3.3)

    3. LINE GAUSS-SEIDEL SWEEP (xi-lines):
       For j = 2 to Nj-1:
         Form tridiagonal system (Section 6.1)
         Solve with Thomas algorithm (Section 6.2)
         Apply SOR relaxation

    4. LINE GAUSS-SEIDEL SWEEP (eta-lines):
       For i = 2 to Ni-1:
         Form tridiagonal system (Section 6.3)
         Solve with Thomas algorithm
         Apply SOR relaxation

    5. CONVERGENCE CHECK:
       error = max(||x^{n+1} - x^n||_inf, ||y^{n+1} - y^n||_inf)
       If error < tol: STOP

OUTPUT: Smoothed grid, final error, iteration count
```

---

## References

1. Steger, J.L. and Sorenson, R.L. (1979). "Automatic mesh-point clustering near a
   boundary in grid generation with elliptic partial differential equations."
   *J. Comp. Phys.* 33, 405-410.

2. Sorenson, R.L. (1980). "A computer program to generate two-dimensional grids about
   airfoils and other shapes by the use of Poisson's equation." NASA TM-81198. (GRAPE)

3. Thompson, J.F., Thames, F.C., and Mastin, C.W. (1974). "Automatic numerical
   generation of body-fitted curvilinear coordinate system for field containing any
   number of arbitrary two-dimensional bodies." *J. Comp. Phys.* 15, 299-319.

4. Thompson, J.F., Warsi, Z.U.A., and Mastin, C.W. (1985). *Numerical Grid Generation:
   Foundations and Applications.* North-Holland.

5. Winslow, A.M. (1967). "Numerical solution of the quasilinear Poisson equation in a
   nonuniform triangle mesh." *J. Comp. Phys.* 1(2), 149-172.

6. Thomas, P.D. and Middlecoff, J.F. (1980). "Direct control of the grid point
   distribution in meshes generated by elliptic equations." *AIAA J.* 18(6), 652-656.

7. Young, D.M. (1954). "Iterative methods for solving partial difference equations of
   elliptic type." *Trans. Amer. Math. Soc.* 76, 92-111.
