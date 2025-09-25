function ThomasAlg(a, b, c, d)
    n = length(d)
    cp = similar(c)
    dp = similar(d)

    # forward sweep
    cp[1] = c[1] / b[1]
    dp[1] = d[1] / b[1]

    for i in 2:n
        denom = b[i] - a[i] * cp[i-1]
        cp[i] = c[i] / denom
        dp[i] = (d[i] - a[i] * dp[i-1]) / denom
    end

    # backward substitution
    x = zeros(n)
    x[end] = dp[end]
    for i = n-1:-1:1
        x[i] = dp[i] - cp[i] * x[i+1]
    end

    return x
end
