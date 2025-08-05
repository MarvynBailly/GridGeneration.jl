# Mapping the 2 Dimensional Problem to 1 Dimension

Since we aim to construct grid in $\Omega \subset \R^2$, we need a method to convert the 2D problem into a 1D problem. While we have the 1D equation

$\begin{align*}
-8 \sigma^4 M^2 x_s^2  x_{ss} - 4 \sigma^4 M M_x x_s^4  -  4  \sigma^2 m M x_{ss}  -  2 \sigma^2 m M_x x_s^2  = 0,
\end{align*}$

we need to reduce the dimension of $x \in \R^2$ and $M \in \R^{2 \times 2}$. This will consistent of four steps:
- [Convert the 2D metric to 1D](./MetricReformulation.md)
- [Convert the 2D curve to 1D](./PointProjection.md)
- [Numerically solve the Spacing ODE](./ODENumericalMethods.md)
- [Convert 1D distribution of points back to curve](./PointProjection.md)

## Results
We look at three metric cases:
- Uniform.
- Clustering at $x=0.0$
- Clustering at $x=1$