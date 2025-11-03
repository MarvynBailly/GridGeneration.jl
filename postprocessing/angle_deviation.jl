function ComputeAngleDeviation(blocks)
    angleDeviations = []
    
    maxblockId = -1
    maxDeviation = -1.0
    maxI = -1
    maxJ = -1

    for (i, block) in enumerate(blocks)
        x = block[1, :, :]
        y = block[2, :, :]
        nrows, ncols = size(x)
        
        # Initialize a matrix to store the deviation for each cell corner.
        # The result will be (nrows-1, ncols-1), so the last row/col will be 0.
        deviations = zeros(nrows, ncols)

        # --- Corrected Vector Calculation ---
        
        # Vector 1: Along the xi-direction (horizontal grid lines)
        # For each cell (i,j), this is the vector from point (i,j) to (i,j+1)
        v1_x = x[1:end-1, 2:end]   - x[1:end-1, 1:end-1]
        v1_y = y[1:end-1, 2:end]   - y[1:end-1, 1:end-1]

        # Vector 2: Along the eta-direction (vertical grid lines)
        # For each cell (i,j), this is the vector from point (i,j) to (i+1,j)
        v2_x = x[2:end,   1:end-1] - x[1:end-1, 1:end-1]
        v2_y = y[2:end,   1:end-1] - y[1:end-1, 1:end-1]

        # --- Vectorized Angle Calculation ---
        
        # Dot product: v1 ⋅ v2
        dot_product = v1_x .* v2_x .+ v1_y .* v2_y

        # Magnitude of v1: ||v1||
        norm1 = sqrt.(v1_x.^2 .+ v1_y.^2)
        
        # Magnitude of v2: ||v2||
        norm2 = sqrt.(v2_x.^2 .+ v2_y.^2)

        # Angle from dot product formula: θ = acos((v1 ⋅ v2) / (||v1|| * ||v2||))
        epsilon = 1e-12
        theta = acosd.(dot_product ./ (norm1 .* norm2))

        # The deviation is the absolute difference from 90 degrees.
        deviations[1:end-1, 1:end-1] = abs.(theta .- 90.0)
        
        # Track the maximum deviation and its block ID
        local_max = maximum(deviations)
        if local_max > maxDeviation
            maxDeviation = local_max
            maxblockId = i
            maxI = findmax(deviations)[2][1]
            maxJ = findmax(deviations)[2][2]
        end

        push!(angleDeviations, deviations)
    end

    return angleDeviations, maxblockId, maxDeviation, maxI, maxJ
end