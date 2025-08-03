# Ordinary Differential Equation for Metric-Conforming Structured Grid Generator
The aim is to create a simple metric-conforming structured grid generator by dividing the entire domain $\Omega$ into sub-domains $\hat{\Omega}_i$ known as blocks. Along the sides of the blocks, we solve a 1D grid-spacing problem through the calculus of variation whose solution gives a coordinate mapping for the distribution of points according to the metric tensor. An algebraic method known as Transfinite Interpolation is used to fill in each block by interpolating the distribution of points along the edges into 2D. Supporting mathematical steps are shown in the [Appendix](MathematicalWork.md).

## General Grid-Spacing Problem
We define the metric

$m_\alpha(x_s, x, s) = \sigma_\alpha^2  M x_s^2 - 1,$

where $M(x(s))$ is the metric tensor field designed from the residual posterior estimate, $x$ is the physical domain and $s$ is the computational domain. We will use variable subscripts to notation derivatives $\left( x_s := \frac{dx}{ds} \right)$. Next, we define our loss as 

$L_\text{misfit}(x_s, x, s) = m_\alpha^2,$ 

which we use to define the functional  

$\mathcal{L}[x(s)] = \int_\Omega L(x(s), x_s(s), s) d s.$

From variations we know the general solution is the Euler-Lagrange equation

$\frac{\partial L}{\partial x}- \frac{d}{d s} \frac{\partial L}{\partial x_s} = 0.$

For this choice of kernel $L_\text{misfit}$, the solution in $\R^n$ is a second-order nonlinear PDE of the form

$-\sum_{\alpha=1}^{n} 8 \sigma_\alpha^4 M_{kl}\frac{\partial x_l}{\partial s_\alpha} M_{ij} \frac{\partial x_i}{\partial s_\alpha} \frac{\partial^2 x_j}{\partial s_\alpha^2} - \sum_{\alpha = 1}^n 4 \sigma_\alpha^4 M_{kl} \frac{\partial x_l}{\partial s_\alpha}\frac{\partial M_{ij}}{\partial x_p} \frac{\partial x_p}{\partial s_\alpha} \frac{\partial x_i}{\partial s_\alpha} \frac{\partial x_j}{\partial s_\alpha} - \sum_{\alpha=1}^n m_\alpha \sigma_\alpha^2\left( 4 M_{kj} \frac{\partial^2 x_j}{\partial s_\alpha^2} + 4 \frac{\partial M_{kj}}{\partial x_p}\frac{\partial x_p}{\partial x_s} \frac{\partial x_j}{\partial s_\alpha}  - 2 \frac{\partial M_{ij}}{\partial x_k}\frac{\partial x_i}{\partial s_\alpha}  \frac{\partial x_j}{\partial s_\alpha} \right) = 0, k \in \{1, 2, \dots, n\}$

where Einstein Notation is used for non-Greek subscripts. 

## 1D Grid-Spacing
In $\R$, the above PDE turns into a friendly nonlinear ODE.

### Second Order ODE
Simplifying the above equation with $n=1$ yields a second-order nonlinear ODE of the form 

$\boxed{ 8 \sigma^4  M^2 x_s^2 x_{ss}  + 4 \sigma^4  M M_x x_s^4 + 4\sigma^2 m M x_{ss} + 2 \sigma^2 m M_x x_s^2  = 0, }$

with dirichlet boundary conditions

$x(0) = x_0, \quad x(1) = x_1$ 

### First Order System of ODEs
We can rewrite the second order ODE as a system of first order ODEs by letting $u_1 = x$ and $u_2 = x_s$ to find

$\begin{cases} 
u_1' = u_2 = x',\\
u'_2 = x'' = - \frac{1}{2} \frac{M_x}{M}x_s^2.
\end{cases}$

with dirichlet boundary conditions:

$u_1(1) = x_1, u_1(0) = x_0, u_2(0) = u_2(1) = 0.$

### Optimal Grid Spacing
We also require that 

$\frac{\partial L}{\partial \sigma} = 0,$ 

which gives us the condition

$\boxed{\sigma = \left( \frac{\int_{\hat{\Omega}} M x^2 ds}{\int_{\hat{\Omega}} (M x^2)^2 ds} \right)^{1/2}.}$

where $\sigma = \frac{1}{n_\text{opt}}$. 


We aim to solve the boxed equation numerically to find the distribution of points along the boundary of the blocks.