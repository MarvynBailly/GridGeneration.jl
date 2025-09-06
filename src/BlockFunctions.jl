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
    # @info "Optimal number of points: $optN"

    solOptA = GridGeneration.SolveODE(metricA, metricA, optN, xsA; method=solver)
    solOptB = GridGeneration.SolveODE(metricB, metricB, optN, xsB; method=solver)

    projectedA = GridGeneration.ProjectBoundary1Dto2D(edgeA, solOptA[1, :], xsA)
    projectedB = GridGeneration.ProjectBoundary1Dto2D(edgeB, solOptB[1, :], xsB)

    return projectedA, projectedB
end

function GetOptNEdgePair(edgeA, edgeB, metricFunc; solver="analytic")
    metricA = GridGeneration.Get1DMetric(edgeA, metricFunc)
    metricB = GridGeneration.Get1DMetric(edgeB, metricFunc)

    xsA = GridGeneration.ProjectBoundary2Dto1D(edgeA)
    xsB = GridGeneration.ProjectBoundary2Dto1D(edgeB)

    metricADiff = GridGeneration.CentralDiff(xsA, metricA)
    metricBDiff = GridGeneration.CentralDiff(xsB, metricB)

    solA = GridGeneration.SolveODE(metricA, metricADiff, length(xsA), xsA; method=solver)
    solB = GridGeneration.SolveODE(metricB, metricBDiff, length(xsB), xsB; method=solver)

    optNA = GridGeneration.ComputeOptimalNumberofPoints(xsA, metricA, solA[1, :])
    optNB = GridGeneration.ComputeOptimalNumberofPoints(xsB, metricB, solB[1, :])

    optN = max(optNA, optNB)
    
    return optN
end

function ProcessEdgePairFixedN(edgeA, edgeB, metricFunc, solver, N)
    metricA = GridGeneration.Get1DMetric(edgeA, metricFunc)
    metricB = GridGeneration.Get1DMetric(edgeB, metricFunc)

    xsA = GridGeneration.ProjectBoundary2Dto1D(edgeA)
    xsB = GridGeneration.ProjectBoundary2Dto1D(edgeB)

    metricADiff = GridGeneration.CentralDiff(xsA, metricA)
    metricBDiff = GridGeneration.CentralDiff(xsB, metricB)

    optN = N
    # @info "Optimal number of points: $optN"

    solOptA = GridGeneration.SolveODE(metricA, metricADiff, optN, xsA; method=solver)
    solOptB = GridGeneration.SolveODE(metricB, metricBDiff, optN, xsB; method=solver)

    projectedA = GridGeneration.ProjectBoundary1Dto2D(edgeA, solOptA[1, :], xsA)
    projectedB = GridGeneration.ProjectBoundary1Dto2D(edgeB, solOptB[1, :], xsB)

    return projectedA, projectedB
end


function SolveBlockFixedN(block, bndInfo, interInfo, metricFunc, optNs; solver="analytic", tfi_method="TFI")
    Ni = optNs[1]
    Nj = optNs[2]

    left   = block[:, 1, :]
    right  = block[:, end, :]
    projectedLeft, projectedRight = ProcessEdgePairFixedN(left, right, metricFunc, solver, Ni)

    bottom = block[:, :, 1]
    top    = block[:, :, end]
    projectedBottom, projectedTop = ProcessEdgePairFixedN(bottom, top, metricFunc, solver, Nj)

    if tfi_method == "TFI"
        computedBlock = GridGeneration.TFI_2D([projectedTop', projectedRight', projectedBottom', projectedLeft'])
    elseif tfi_method == "TFI_2D_Hermite"
        computedBlock = GridGeneration.TFI_2D_Hermite([projectedTop', projectedRight', projectedBottom', projectedLeft'])
    else
        error("Unknown TFI method: $tfi_method")
    end

    return computedBlock, bndInfo, interInfo
end


function SolveBlock(block, bndInfo, interInfo, metricFunc; solver="analytic", tfi_method="TFI")
    left   = block[:, 1, :]
    right  = block[:, end, :]
    projectedLeft, projectedRight = ProcessEdgePair(left, right, metricFunc, solver)

    bottom = block[:, :, 1]
    top    = block[:, :, end]
    projectedBottom, projectedTop = ProcessEdgePair(bottom, top, metricFunc, solver)

    if tfi_method == "TFI"
        computedBlock = GridGeneration.TFI_2D([projectedTop', projectedRight', projectedBottom', projectedLeft'])
    elseif tfi_method == "TFI_2D_Hermite"
        computedBlock = GridGeneration.TFI_2D_Hermite([projectedTop', projectedRight', projectedBottom', projectedLeft'])
    else
        error("Unknown TFI method: $tfi_method")
    end

    bndInfo = GridGeneration.UpdateBndInfo!(bndInfo, computedBlock)

    return computedBlock, bndInfo, interInfo
end
