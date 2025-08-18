```@meta
CurrentModule = GridGeneration
```

# GridGeneration.jl Documentation

Welcome to the documentation and formulation for GridGeneration.jl

## Current Example

Just finished adding examples for [airfoil case](./pages/Examples/airfoil.md)


## Overview

A brief description of the underlying ordinary differential equation (ODE) is presented in [ODE Formulation](./pages/ODE/ODEFormulation.md) with supporting work shown in [Mathematical Work](./pages/ODE/MathematicalWork.md). The ODE is nonlinear and second order boundary value problem which can be reformulated as a system system of first order ODEs. Two numerical methods, one for the [First Order System](./pages/NumericalMethods/FirstOrderSystem.md) using Julia's library [DifferentialEquations.jl](./pages/NumericalMethods/FirstOrderSystem.md) and second for the [Second Order BVP ODE](./pages/NumericalMethods/SecondOrderBVP.md) using central differencing and fixed point iteration with under-relaxation. Both methods prove to be rather unstable so a [semi-analytical method](./pages/NumericalMethods/SemiAnalyticalMethod.md) (semi due to the use of numerical integration and inversion) is adopted. 

To create 2D and 3D grids, we present a method of [Mapping 2D to 1D](./pages/2Dto1D/Mapping2Dto1D.md) and 1D back to 2D. 

### 
###
###
###
Work is down under the supervision and support of [Dr. Larsson](https://larsson.umd.edu) at the University of Maryland.