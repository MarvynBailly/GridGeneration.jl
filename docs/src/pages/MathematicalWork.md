# Mathematical Work
## PDE Formulation
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

## ODE Formulation
### Second Order Nonlinear ODE
In one dimension, we can clean up the above equation notably by removing all the summations and indices:

$\begin{align*}
-8 \sigma^4 M\frac{\partial x}{\partial s} M \frac{\partial x}{\partial s} \frac{\partial^2 x}{\partial s^2} - 4 \sigma^4 M \frac{\partial x}{\partial s}\frac{\partial M}{\partial x} \frac{\partial x}{\partial s}\frac{\partial x}{\partial s} \frac{\partial x}{\partial s} -  m \sigma^2\left( 4 M\frac{\partial^2 x}{\partial s^2} + 4 \frac{\partial M}{\partial x}\frac{\partial x}{\partial s}\frac{\partial x}{\partial s}  - 2 \frac{\partial M}{\partial x}\frac{\partial x}{\partial s}  \frac{\partial x}{\partial s} \right) = 0
\end{align*}$

Let's combine like terms, fully expand, and introduce the subscript notation for partial derivatives

$\begin{align*}
-8 \sigma^4 M^2 x_s^2  x_{ss} - 4 \sigma^4 M M_x x_s^4  -  4  \sigma^2 m M x_{ss}  -  2 \sigma^2 m M_x x_s^2  = 0,
\end{align*}$

which is our second order nonlinear ODE.

Working towards a first order system, let's simplify by dividing both sides by $-2 \sigma^2$ to get

$\begin{align*}
4 \sigma^2 M^2 x_s^2  x_{ss} + 2 \sigma^2 M M_x x_s^4  +  2  m M x_{ss}  +  m M_x x_s^2  = 0,
\end{align*}$

Now let's solve for $x_{ss}$ to find

$\begin{align*}
\left( 4 \sigma^2 M^2 x_s^2 + 2 m M \right) x_{ss} + 2 \sigma^2 M M_x x_s^4 +  m M_x x_s^2  = 0 \iff x_{ss} = - \frac{2 \sigma^2 M M_x x_s^4 +  m M_x x_s^2}{4 \sigma^2 M^2 x_s^2 + 2 m M}
\end{align*}$

Further simplifying we arrive at

$\begin{align*}
x_{ss} = - \frac{ M_x x_s^2 \left( 2\sigma^2 M x_s^2 +  m \right)}{2M \left( 2\sigma^2 M x_s^2 + m\right)} = - \frac{ M_x x_s^2}{2M}.
\end{align*}$


### First Order System

To turn the above into a first order system, let $u_1 = x$ and $u_2 = x_s$, then

$\begin{cases} 
u_1' = u_2 = x',\\
u'_2 = x'' = - \frac{1}{2} \frac{M_x}{M} x_s^2.
\end{cases}$

### Semi-Analytic Solution

We can continue to simplify by letting $u = x_s$, then we arrive at

$\begin{align*}
x_{ss} = - \frac{ M_x x_s^2}{2M} \implies u' = - \frac{1}{2} \frac{M_x}{M} u.
\end{align*}$

If we assume that $x_s \neq 0$, we can divide both sides by $u$ and integrate to find that

$\ln(p) = - \frac{1}{2} \ln(M) \iff p = x_s = \frac{C_1}{\sqrt{M(x)}},$

where $C_1$ comes from the constant of integration. 

We have an initial value boundary value problem (IVBP ODE), so let's prescribe boundary conditions at $x(s=0) = 0$ and $x(s=1) = L$ such that our computational domain is $s \in [0,1]$ and the physical domain is $x \in [0,L] \subset \R$.

Next, let's separate variables to get

$\frac{dx}{ds} = \frac{C_1}{\sqrt{M(x)}} \rightarrow \sqrt{M(x)} dx = C_1 ds.$

Integrating both sides from $s=0$ to $s$, 

$\int_{x(0)}^{x(s)} \sqrt{M(\xi)} d\xi = C_1 s$

Finally let's enforce the boundary condition at $x(1) = L$ to solve for $C_1$ and find that

$C_1 = \int_{0}^{L} \sqrt{M(\xi)} d \xi = I.$

Therefore the final solution becomes

$\int_{0}^{x(s)} \sqrt{M(\xi)} d\xi = I s.$

If we let $I(x) = \int_0^x \sqrt{M(\xi)} d \xi$, we can further clean up the express as

$I(x(s)) = s I(x(s=1))$