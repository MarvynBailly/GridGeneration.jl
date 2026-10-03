"""
    ComputeAngleDeviation(blocks) -> (angleDeviations, maxBlockId, maxDeviation, maxI, maxJ)

Measure grid orthogonality. For every cell `(i, j)` of each `[2, Ni, Nj]` block, compute the
angle at its corner node `(i, j)` between the grid lines towards `(i+1, j)` and `(i, j+1)`, and
return its deviation from 90° in degrees.

- `angleDeviations[b]` is an `Ni×Nj` matrix for block `b` (the last row and column, which
  have no cell, are zero).
- `maxDeviation` is the largest deviation over all blocks, found in block `maxBlockId` at
  cell `(maxI, maxJ)`.

Degenerate cells (a zero-length edge) are reported as `NaN` and ignored for the maximum.
"""
function ComputeAngleDeviation(blocks)
    angleDeviations = Matrix{Float64}[]

    maxBlockId = -1
    maxDeviation = -1.0
    maxI = -1
    maxJ = -1

    for (b, block) in enumerate(blocks)
        x = block[1, :, :]
        y = block[2, :, :]
        nrows, ncols = size(x)
        deviations = zeros(nrows, ncols)

        @inbounds for j in 1:ncols-1, i in 1:nrows-1
            # grid line towards (i, j+1) and towards (i+1, j)
            v1x = x[i, j+1] - x[i, j];  v1y = y[i, j+1] - y[i, j]
            v2x = x[i+1, j] - x[i, j];  v2y = y[i+1, j] - y[i, j]
            n = sqrt(v1x^2 + v1y^2) * sqrt(v2x^2 + v2y^2)
            if n == 0
                deviations[i, j] = NaN
                continue
            end
            # clamp: rounding can push the cosine just outside [-1, 1]
            θ = acosd(clamp((v1x * v2x + v1y * v2y) / n, -1.0, 1.0))
            dev = abs(θ - 90.0)
            deviations[i, j] = dev
            if dev > maxDeviation
                maxDeviation, maxBlockId, maxI, maxJ = dev, b, i, j
            end
        end

        push!(angleDeviations, deviations)
    end

    return angleDeviations, maxBlockId, maxDeviation, maxI, maxJ
end
