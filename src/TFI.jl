"""
2D Transfinite interpolation (Coons patch).
Input boundary with order: (top, right, bottom, left), each as an N×2 array.
Returns X', Y' (matching your original orientation).
"""
function TFI_2D(boundary)
    top, right, bottom, left = boundary

    N2 = size(left,   1)   # points along η (vertical)
    N1 = size(bottom, 1)   # points along ξ (horizontal)

    # Parameter lines as row/column for broadcast-friendly shapes
    ξ = reshape(range(0.0, 1.0; length=N1), 1, :)   # 1×N1
    η = reshape(range(0.0, 1.0; length=N2), :, 1)   # N2×1

    # Convenience shorthands
    A = 1 .- ξ          # 1×N1
    B = ξ               # 1×N1
    C = 1 .- η          # N2×1
    D = η               # N2×1

    # Split boundary coordinates
    Lx, Ly = left[:,1],   left[:,2]
    Rx, Ry = right[:,1],  right[:,2]
    Bx, By = bottom[:,1], bottom[:,2]
    Tx, Ty = top[:,1],    top[:,2]

    # Corner points (consistent with your original code using left/right corners)
    x00, y00 = left[1,1],  left[1,2]   # η=0, ξ=0
    x10, y10 = right[1,1], right[1,2]  # η=0, ξ=1
    x01, y01 = left[end,1], left[end,2]   # η=1, ξ=0
    x11, y11 = right[end,1], right[end,2] # η=1, ξ=1

    # X (N2×N1) via broadcasting
    X =  A .* Lx .+ B .* Rx .+ C .* Bx' .+ D .* Tx' .-
         (A .* C .* x00 .+ B .* C .* x10 .+ A .* D .* x01 .+ B .* D .* x11)

    # Y (N2×N1)
    Y =  A .* Ly .+ B .* Ry .+ C .* By' .+ D .* Ty' .-
         (A .* C .* y00 .+ B .* C .* y10 .+ A .* D .* y01 .+ B .* D .* y11)

    return X', Y'
end


"""
2D Transfinite interpolation. Input boundary with (top, right, bottom, left) order. Outputs the X and Y coordinates of the grid points.
"""
function TFI_2D_old(boundary)
    top = boundary[1]
    right = boundary[2]
    bottom = boundary[3]
    left = boundary[4]

    N2 = size(left, 1)   
    N1 = size(bottom, 1)  

    s1 = range(0, 1; length=N1)
    s2 = range(0, 1; length=N2)

    X = zeros(N2, N1)
    Y = zeros(N2, N1)

    for j in 1:N2
        for i in 1:N1
            
            ξ = s1[i]
            η = s2[j]
            
            # Interpolated values from edges
            x = (1 - ξ) * left[j, 1] + ξ * right[j, 1] +
                (1 - η) * bottom[i, 1] + η * top[i, 1] -
                ((1 - ξ)*(1 - η) * left[1, 1] + ξ*(1 - η) * right[1, 1] +
                 (1 - ξ)*η * left[end, 1] + ξ*η * right[end, 1])

            y = (1 - ξ) * left[j, 2] + ξ * right[j, 2] +
                (1 - η) * bottom[i, 2] + η * top[i, 2] -
                ((1 - ξ)*(1 - η) * left[1, 2] + ξ*(1 - η) * right[1, 2] +
                 (1 - ξ)*η * left[end, 2] + ξ*η * right[end, 2])

            X[j, i] = x
            Y[j, i] = y            
        end
    end

    return X', Y'
end