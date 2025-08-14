# Single Block with No Splitting
## 2D Single Block
Let's start by allowing the user to input a single [Tortuga](../GridFormat.md) block in the code. Thus the input will be the initial grid and the boundary information. The code will take in a 2D grid and solve the ODE along each edge of the domain. As the optimal number of points need not be the same for boundary edges across from each other, the code will proceed to resolve the edge that does not have the maximum number of points. Finally, the code will solve for the interior points using Transfinite Interpolation and update the boundary information with the new dimensions of the block. 

## Algorithm
To solve a single block, we can pair up the left and top edges with the right and bottom edges respectively. For each pair, we do the following:
- for each pair:
  - for each edge in pair
    - Compute the 1D metric along the edge using `GridGeneration.Get1DMetric(pnts, metricFunction)`
    - Project the edge to 1D using `GridGeneration.ProjectBoundary2Dto1D(edge)`
    - Solve the ODE problem using the 1D metric and 1D points using `GridGeneration.SolveODE(m,m,N,xs)`
    - Compute the optimal number of points using `GridGeneration.ComputeOptimalNumberofPoints(xs, m, sol)`
  - set optimal number of points for pair to the max of the edge optimal numbers
  - for each edge
    - Solve the ODE problem using the pair optimal number of points using `GridGeneration.ComputeOptimalNumberofPoints(xs, m, sol)`
    - Save edge distribution
- Run TFI on edge distributions to get final grid using `GridGeneration.TFI_2D(boundary)`

Writing this in Julia:

```julia
function ProcessEdgePair(edgeA, edgeB, metricFunc, solver)
    metricA = GridGeneration.Get1DMetric(edgeA, metricFunc)
    metricB = GridGeneration.Get1DMetric(edgeB, metricFunc)

    xsA = GridGeneration.ProjectBoundary2Dto1D(edgeA)
    xsB = GridGeneration.ProjectBoundary2Dto1D(edgeB)

    solA = GridGeneration.SolveODE(metricA, metricA, length(xsA), xsA; method=solver)
    solB = GridGeneration.SolveODE(metricB, metricB, length(xsB), xsB; method=solver)

    optNA = GridGeneration.ComputeOptimalNumberofPoints(xsA, metricA, solA[1, :])
    optNB = GridGeneration.ComputeOptimalNumberofPoints(xsB, metricB, solB[1, :])

    optN = max(optNA, optNB)
    @info "Optimal number of points: $optN"

    solOptA = GridGeneration.SolveODE(metricA, metricA, optN, xsA; method=solver)
    solOptB = GridGeneration.SolveODE(metricB, metricB, optN, xsB; method=solver)

    projectedA = GridGeneration.ProjectBoundary1Dto2D(edgeA, solOptA[1, :], xsA)
    projectedB = GridGeneration.ProjectBoundary1Dto2D(edgeB, solOptB[1, :], xsB)

    return projectedA, projectedB
end

function SolveBlock(block, bndInfo, interInfo, metricFunc; solver="analytic")
    left   = block[:, 1, :]
    right  = block[:, end, :]
    projectedLeft, projectedRight = ProcessEdgePair(left, right, metricFunc, solver)

    bottom = block[:, :, 1]
    top    = block[:, :, end]
    projectedBottom, projectedTop = ProcessEdgePair(bottom, top, metricFunc, solver)

    computedBlock = GridGeneration.TFI_2D([projectedTop', projectedRight', projectedBottom', projectedLeft'])

    bndInfo = GridGeneration.UpdateBndInfo!(bndInfo, computedBlock)

    return computedBlock, bndInfo, interInfo
end
```

### Examples 
Let's use a small grid around an airfoil as an example. Since we are doing single block, let's not have a c-grid yet but keep one block around

