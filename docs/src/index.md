```@meta
CurrentModule = GridGeneration
```

# GridGeneration.jl

GridGeneration.jl generates structured multi-block grids whose point spacing follows a
user-supplied metric field. Starting from an initial single- or multi-block grid it

1. **splits** blocks at chosen indices ([`SplitMultiBlock`](@ref), [Single-block splitting](./pages/SingleBlock/splitting.md)),
2. **redistributes** each block's edge points so that the spacing equidistributes the metric, choosing the
   number of points per direction automatically ([`SolveAllBlocks`](@ref)),
3. **rebuilds** each block by transfinite interpolation ([`TFI`](@ref)), and
4. **smooths** the blocks with an elliptic (Winslow-type) solver with wall-orthogonality forcing ([`SmoothBlocks`](@ref)).

[`GenerateGrid`](@ref) runs the whole pipeline. Grids can be read from and written to the
Tortuga (`.grid`) format.

New here? Start with [Getting Started](./pages/GettingStarted.md), then see the
[airfoil example](./pages/Examples/airfoil.md) and the [API Reference](./pages/api.md).

## Overview of the method

A brief description of the underlying ordinary differential equation (ODE) is presented in
[ODE Formulation](./pages/ODE/ODEFormulation.md) with supporting work shown in
[Mathematical Work](./pages/ODE/MathematicalWork.md). The ODE is a nonlinear second-order
boundary value problem, which can be reformulated as a system of first-order ODEs. Two numerical
methods were tried: one for the [First Order System](./pages/NumericalMethods/FirstOrderSystem.md)
using [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/), and one for the
[Second Order BVP ODE](./pages/NumericalMethods/SecondOrderBVP.md) using central differencing and
fixed-point iteration with under-relaxation. Both proved fragile, so a
[semi-analytical method](./pages/NumericalMethods/SemiAnalyticalMethod.md) (semi due to the use of
numerical integration and inversion) is the default (`boundarySolver = :analytic`).

To create 2D and 3D grids, we present a method of [mapping 2D to 1D](./pages/2Dto1D/Mapping2Dto1D.md)
and projecting the 1D solution back to 2D. The elliptic smoother is described in the
[Theory](./pages/Theory/elliptic_smoothing.md) section.

## Acknowledgements

Work is done under the supervision and support of [Dr. Larsson](https://larsson.umd.edu) at the University of Maryland.
