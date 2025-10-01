function GetOptNEdgePair(edgeA, edgeB, M; solver=:analytic)
    metricA = GridGeneration.Get1DMetric(edgeA, M)
    metricB = GridGeneration.Get1DMetric(edgeB, M)

    xsA = GridGeneration.ProjectBoundary2Dto1D(edgeA)
    xsB = GridGeneration.ProjectBoundary2Dto1D(edgeB)

    metricFuncA = GridGeneration.LinearInterpolate(xsA, metricA)
    metricFuncB = GridGeneration.LinearInterpolate(xsB, metricB)

    solA = GridGeneration.SolveODE(metricFuncA, xsA; solver=solver)
    solB = GridGeneration.SolveODE(metricFuncB, xsB; solver=solver)
    
    optNA = GridGeneration.ComputeOptimalNumberofPoints(solA, metricFuncA)
    optNB = GridGeneration.ComputeOptimalNumberofPoints(solB, metricFuncB)

    optN = max(optNA, optNB)
    @info "Optimal N for edge pair: $optN"
    
    return optN
end

function ProcessEdgePairFixedN(edgeA, edgeB, M, N; solver=:analytic)
    metricA = GridGeneration.Get1DMetric(edgeA, M)
    metricB = GridGeneration.Get1DMetric(edgeB, M)

    xsA = GridGeneration.ProjectBoundary2Dto1D(edgeA)
    xsB = GridGeneration.ProjectBoundary2Dto1D(edgeB)
    
    metricFuncA = GridGeneration.LinearInterpolate(xsA, metricA)
    metricFuncB = GridGeneration.LinearInterpolate(xsB, metricB)

    solOptA = GridGeneration.SolveODEFixedN(metricFuncA, xsA, N; solver=solver)
    solOptB = GridGeneration.SolveODEFixedN(metricFuncB, xsB, N; solver=solver)

    projectedA = GridGeneration.ProjectBoundary1Dto2D(edgeA, solOptA)
    projectedB = GridGeneration.ProjectBoundary1Dto2D(edgeB, solOptB)

    return projectedA, projectedB
end


function SolveBlockFixedN(block, bndInfo, interInfo, metricFunc, optNs; solver=:analytic, tfi_method=:TFI)
    Ni = optNs[1]
    Nj = optNs[2]

    left   = block[:, 1, :]
    right  = block[:, end, :]
    projectedLeft, projectedRight = ProcessEdgePairFixedN(left, right, metricFunc, Ni)

    bottom = block[:, :, 1]
    top    = block[:, :, end]
    projectedBottom, projectedTop = ProcessEdgePairFixedN(bottom, top, metricFunc, Nj)

    if tfi_method == :TFI
        computedBlock = GridGeneration.TFI([projectedTop', projectedRight', projectedBottom', projectedLeft'])
    elseif tfi_method == :TFI_2D_Hermite
        computedBlock = GridGeneration.TFI_2D_Hermite([projectedTop', projectedRight', projectedBottom', projectedLeft'])
    else
        error("Unknown TFI method: $tfi_method")
    end

    return computedBlock, bndInfo, interInfo
end