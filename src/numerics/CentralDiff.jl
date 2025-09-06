function CentralDiff(x, f)
    n = length(f)
    df = similar(f)

    # left boundary (forward diff)
    df[1] = (f[2] - f[1]) / (x[2] - x[1])

    # interior points (vectorized central diff on nonuniform grid)
    df[2:n-1] .= (f[3:n] .- f[1:n-2]) ./ (x[3:n] .- x[1:n-2])

    # right boundary (backward diff)
    df[n] = (f[n] - f[n-1]) / (x[n] - x[n-1])

    return df
end