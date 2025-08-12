"""
2D Transfinite interpolation. Input boundary with (top, right, bottom, left) order. Outputs the X and Y coordinates of the grid points.
"""
function TFI_2D(boundary)
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