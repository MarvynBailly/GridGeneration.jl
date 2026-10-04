using LinearAlgebra

"""
Calculates the metric terms α, β, γ and the Jacobian J.
Uses central differences for the interior and one-sided for boundaries.
"""
function calculate_metrics(x::Matrix, y::Matrix)
    Ni, Nj = size(x)
    alpha = zeros(Ni, Nj); beta = zeros(Ni, Nj); gamma = zeros(Ni, Nj); J = zeros(Ni, Nj)
    calculate_metrics!(alpha, beta, gamma, x, y; J=J)
    return alpha, beta, gamma, J
end

"""
In-place version of [`calculate_metrics`](@ref): fills `alpha`, `beta`, `gamma` (and `J` if given).
"""
function calculate_metrics!(alpha, beta, gamma, x::Matrix, y::Matrix; J=nothing)
    Ni, Nj = size(x)
    @inbounds for j in 1:Nj, i in 1:Ni
        # ξ-derivatives
        if i == 1
            x_xi = x[i+1, j] - x[i, j];  y_xi = y[i+1, j] - y[i, j]
        elseif i == Ni
            x_xi = x[i, j] - x[i-1, j];  y_xi = y[i, j] - y[i-1, j]
        else
            x_xi = (x[i+1, j] - x[i-1, j]) / 2.0;  y_xi = (y[i+1, j] - y[i-1, j]) / 2.0
        end
        # η-derivatives
        if j == 1
            x_eta = x[i, j+1] - x[i, j];  y_eta = y[i, j+1] - y[i, j]
        elseif j == Nj
            x_eta = x[i, j] - x[i, j-1];  y_eta = y[i, j] - y[i, j-1]
        else
            x_eta = (x[i, j+1] - x[i, j-1]) / 2.0;  y_eta = (y[i, j+1] - y[i, j-1]) / 2.0
        end

        alpha[i, j] = x_eta^2 + y_eta^2
        beta[i, j]  = x_xi * x_eta + y_xi * y_eta
        gamma[i, j] = x_xi^2 + y_xi^2
        J === nothing || (J[i, j] = x_xi * y_eta - x_eta * y_xi)
    end
    return alpha, beta, gamma
end

"""
Add the forcing for the bottom/top wall (enforcing orthogonality and spacing) to `RHS_x`/`RHS_y`.
wall=1 for bottom (j=1), wall=Nj for top (j=Nj). `decay_x[d+1]`/`decay_y[d+1]` are the
decay factors at distance `d` from the wall.
"""
function add_forcing_eta!(RHS_x, RHS_y, x::Matrix, y::Matrix, decay_x, decay_y; wall::Int = 1)
    Ni, Nj = size(x)
    @assert wall == Nj || wall == 1 "wall must be either 1 (bottom) or Nj (top)"
    dir = wall == 1 ? 1 : -1

    # Desired spacing interpolated between the current spacings at the wall ends
    @inbounds s1 = sqrt((x[1, wall] - x[1, wall + dir*1])^2 + (y[1, wall] - y[1, wall + dir*1])^2)
    @inbounds s2 = sqrt((x[Ni, wall] - x[Ni, wall + dir*1])^2 + (y[Ni, wall] - y[Ni, wall + dir*1])^2)
    s_vec = LinRange(s1, s2, Ni)

    @inbounds for i in 2:Ni-1
        # ξ-derivatives along boundary
        x_xi = (x[i+1, wall] - x[i-1, wall]) / 2.0
        y_xi = (y[i+1, wall] - y[i-1, wall]) / 2.0
        x_xixi = x[i+1, wall] - 2*x[i, wall] + x[i-1, wall]
        y_xixi = y[i+1, wall] - 2*y[i, wall] + y[i-1, wall]

        # Impose orthogonality and spacing to find desired η-derivatives
        ds_inv = 1.0 / sqrt(x_xi^2 + y_xi^2)
        x_eta_desired = -s_vec[i] * y_xi * ds_inv
        y_eta_desired = s_vec[i] * x_xi * ds_inv

        # η second derivatives from a one-sided second-order formula
        x_etaeta = 0.5*(-7*x[i, wall] + 8*x[i, wall+dir*1] - x[i, wall+dir*2]) - dir*3*x_eta_desired
        y_etaeta = 0.5*(-7*y[i, wall] + 8*y[i, wall+dir*1] - y[i, wall+dir*2]) - dir*3*y_eta_desired

        # Boundary metrics: α = |r_η|² = s² (desired wall-normal spacing), γ = |r_ξ|² (along the wall)
        alpha_b = x_eta_desired^2 + y_eta_desired^2
        gamma_b = x_xi^2 + y_xi^2

        rx = -(alpha_b * x_xixi + gamma_b * x_etaeta)
        ry = -(alpha_b * y_xixi + gamma_b * y_etaeta)

        # Propagate into the domain with exponential decay
        for j in 1:Nj
            d = abs(j - wall) + 1
            RHS_x[i, j] += rx * decay_x[d]
            RHS_y[i, j] += ry * decay_y[d]
        end
    end
    return RHS_x, RHS_y
end

"""
Add the forcing for the left/right wall (enforcing orthogonality and spacing) to `RHS_x`/`RHS_y`.
wall=1 for left (i=1), wall=Ni for right (i=Ni).
"""
function add_forcing_xi!(RHS_x, RHS_y, x::Matrix, y::Matrix, decay_x, decay_y; wall::Int = 1)
    Ni, Nj = size(x)
    @assert wall == 1 || wall == Ni "wall must be either 1 (left) or Ni (right)"
    dir = wall == 1 ? 1 : -1

    @inbounds s1 = sqrt((x[wall, 1] - x[wall + dir*1, 1])^2 + (y[wall, 1] - y[wall + dir*1, 1])^2)
    @inbounds s2 = sqrt((x[wall, Nj] - x[wall + dir*1, Nj])^2 + (y[wall, Nj] - y[wall + dir*1, Nj])^2)
    s_vec = LinRange(s1, s2, Nj)

    @inbounds for j in 2:Nj-1
        s = s_vec[j]

        # η-derivatives along boundary
        x_eta   = (x[wall, j+1] - x[wall, j-1]) / 2.0
        y_eta   = (y[wall, j+1] - y[wall, j-1]) / 2.0
        x_etaeta = x[wall, j+1] - 2*x[wall, j] + x[wall, j-1]
        y_etaeta = y[wall, j+1] - 2*y[wall, j] + y[wall, j-1]

        # Impose orthogonality and spacing to find desired ξ-derivatives
        ds_inv = 1.0 / sqrt(x_eta^2 + y_eta^2)
        x_xi_desired = -s * y_eta * ds_inv
        y_xi_desired = s * x_eta * ds_inv

        # ξ second derivatives from a one-sided second-order formula
        x_xixi = 0.5*(-7*x[wall, j] + 8*x[wall+dir*1, j] - x[wall+dir*2, j]) + dir*3*x_xi_desired
        y_xixi = 0.5*(-7*y[wall, j] + 8*y[wall+dir*1, j] - y[wall+dir*2, j]) + dir*3*y_xi_desired

        # Boundary metrics: α = |r_η|² (along the wall), γ = |r_ξ|² = s² (desired wall-normal spacing)
        alpha_b = x_eta^2 + y_eta^2
        gamma_b = s^2

        rx = -(alpha_b * x_xixi + gamma_b * x_etaeta)
        ry = -(alpha_b * y_xixi + gamma_b * y_etaeta)

        for i in 1:Ni
            d = abs(i - wall) + 1
            RHS_x[i, j] += rx * decay_x[d]
            RHS_y[i, j] += ry * decay_y[d]
        end
    end
    return RHS_x, RHS_y
end

# In-place Thomas algorithm for a tridiagonal system with constant-length work vectors.
# Solves a[k] u[k-1] + b[k] u[k] + c[k] u[k+1] = d[k] for k = 1..n (a[1], c[n] unused).
function thomas!(u, a, b, c, d, cp, dp, n)
    @inbounds begin
        cp[1] = c[1] / b[1]
        dp[1] = d[1] / b[1]
        for k in 2:n
            m = b[k] - a[k] * cp[k-1]
            cp[k] = c[k] / m
            dp[k] = (d[k] - a[k] * dp[k-1]) / m
        end
        u[n] = dp[n]
        for k in n-1:-1:1
            u[k] = dp[k] - cp[k] * u[k+1]
        end
    end
    return u
end

# One line Gauss-Seidel sweep of the discrete grid equations along ξ-lines (dir = 1: each
# j-line solved implicitly in i) or η-lines (dir = 2). Uses the frozen metrics/forcing of the
# current iteration and the latest values on neighbouring lines; relaxes with ω.
function line_sweep!(f, alpha, beta, gamma, RHS, ω, dir, w)
    Ni, Nj = size(f)
    a, b, c, d, u, cp, dp = w
    nLines, n = dir == 1 ? (Nj - 2, Ni - 2) : (Ni - 2, Nj - 2)
    @inbounds for line in 1:nLines
        for k in 1:n
            i, j = dir == 1 ? (k + 1, line + 1) : (line + 1, k + 1)
            cross = 0.5 * beta[i,j] * (f[i+1,j+1] - f[i-1,j+1] - f[i+1,j-1] + f[i-1,j-1])
            b[k] = 2.0 * (alpha[i,j] + gamma[i,j])
            if dir == 1   # implicit in i (coefficient α), explicit in j (coefficient γ)
                a[k] = -alpha[i,j]; c[k] = -alpha[i,j]
                d[k] = gamma[i,j] * (f[i,j+1] + f[i,j-1]) - cross + RHS[i,j]
                k == 1 && (d[k] += alpha[i,j] * f[1, j])
                k == n && (d[k] += alpha[i,j] * f[Ni, j])
            else          # implicit in j (coefficient γ), explicit in i (coefficient α)
                a[k] = -gamma[i,j]; c[k] = -gamma[i,j]
                d[k] = alpha[i,j] * (f[i+1,j] + f[i-1,j]) - cross + RHS[i,j]
                k == 1 && (d[k] += gamma[i,j] * f[i, 1])
                k == n && (d[k] += gamma[i,j] * f[i, Nj])
            end
        end
        thomas!(u, a, b, c, d, cp, dp, n)
        for k in 1:n
            i, j = dir == 1 ? (k + 1, line + 1) : (line + 1, k + 1)
            f[i, j] = (1 - ω) * f[i, j] + ω * u[k]
        end
    end
    return f
end

# One point SOR sweep of the discrete grid equations over the interior points.
function point_sweep!(x, y, alpha, beta, gamma, RHS_x, RHS_y, ω)
    Ni, Nj = size(x)
    @inbounds for j in 2:Nj-1
        for i in 2:Ni-1
            denom = 2.0 * (alpha[i,j] + gamma[i,j])

            # --- X Equation ---
            term_xixi_x = alpha[i,j] * (x[i+1,j] + x[i-1,j])
            term_etaeta_x = gamma[i,j] * (x[i,j+1] + x[i,j-1])
            term_xieta_x = 0.5 * beta[i,j] * (x[i+1,j+1] - x[i-1,j+1] - x[i+1,j-1] + x[i-1,j-1])

            x_new = (term_xixi_x + term_etaeta_x - term_xieta_x + RHS_x[i,j]) / denom
            x[i, j] = (1 - ω) * x[i, j] + ω * x_new

            # --- Y Equation ---
            term_xixi_y = alpha[i,j] * (y[i+1,j] + y[i-1,j])
            term_etaeta_y = gamma[i,j] * (y[i,j+1] + y[i,j-1])
            term_xieta_y = 0.5 * beta[i,j] * (y[i+1,j+1] - y[i-1,j+1] - y[i+1,j-1] + y[i-1,j-1])

            y_new = (term_xixi_y + term_etaeta_y - term_xieta_y + RHS_y[i,j]) / denom
            y[i, j] = (1 - ω) * y[i, j] + ω * y_new
        end
    end
    return x, y
end

"""
    EllipticSolver(x, y; params) -> (x, y, error, iterations)

Smooth a block with the elliptic (Winslow-type) grid equations and optional wall forcing.
`x` and `y` are `Ni×Nj` coordinate matrices; boundary nodes are kept fixed and the interior is
updated in place. Iterates until the largest coordinate change `max(‖Δx‖, ‖Δy‖)` drops below
`params.tol` or `params.max_iter` is reached.

`params.sweep` selects the iteration: `:point` (point SOR, the default) or `:line` (alternating
ξ- and η-line Gauss-Seidel with tridiagonal solves). Both solve the same discrete equations;
`:line` typically needs far fewer iterations and tolerates ω ≈ 1.
"""
function EllipticSolver(x::Matrix, y::Matrix; params)
    verbose_print(msg) = params.verbose && println(msg)

    Ni, Nj = size(x)
    err = 0.0
    finalIter = 0

    # work buffers, reused every iteration
    RHS_x = zeros(Ni, Nj); RHS_y = zeros(Ni, Nj)
    alpha = zeros(Ni, Nj); beta = zeros(Ni, Nj); gamma = zeros(Ni, Nj)
    x_old = similar(x); y_old = similar(y); dx = similar(x); dy = similar(y)

    # decay factors exp(-a d) for d = 0, 1, ... (index d+1)
    decay(a, n) = [exp(-a * d) for d in 0:n-1]
    dTop    = (decay(params.a_decay_top, Nj),    decay(params.b_decay_top, Nj))
    dBottom = (decay(params.a_decay_bottom, Nj), decay(params.b_decay_bottom, Nj))
    dLeft   = (decay(params.a_decay_left, Ni),   decay(params.b_decay_left, Ni))
    dRight  = (decay(params.a_decay_right, Ni),  decay(params.b_decay_right, Ni))
    ω = params.ω
    line = params.sweep == :line
    nmax = max(Ni, Nj)
    lineWork = ntuple(_ -> zeros(nmax), 7)   # a, b, c, d, u, cp, dp for the Thomas solves

    verbose_print("Starting $(line ? "line Gauss-Seidel" : "SOR") iterations...")
    for iter in 1:params.max_iter
        finalIter = iter

        # wall forcing from the current grid
        fill!(RHS_x, 0.0); fill!(RHS_y, 0.0)
        params.useTopWall    && add_forcing_eta!(RHS_x, RHS_y, x, y, dTop...;    wall = Nj)
        params.useBottomWall && add_forcing_eta!(RHS_x, RHS_y, x, y, dBottom...; wall = 1)
        params.useLeftWall   && add_forcing_xi!(RHS_x, RHS_y, x, y, dLeft...;    wall = 1)
        params.useRightWall  && add_forcing_xi!(RHS_x, RHS_y, x, y, dRight...;   wall = Ni)

        copyto!(x_old, x); copyto!(y_old, y)

        # metric terms from the current grid
        calculate_metrics!(alpha, beta, gamma, x, y)

        if line
            # alternating ξ- and η-line sweeps with the same frozen coefficients
            line_sweep!(x, alpha, beta, gamma, RHS_x, ω, 1, lineWork)
            line_sweep!(y, alpha, beta, gamma, RHS_y, ω, 1, lineWork)
            line_sweep!(x, alpha, beta, gamma, RHS_x, ω, 2, lineWork)
            line_sweep!(y, alpha, beta, gamma, RHS_y, ω, 2, lineWork)
        else
            point_sweep!(x, y, alpha, beta, gamma, RHS_x, RHS_y, ω)
        end


        # Check for convergence (change in both coordinates)
        dx .= x .- x_old; dy .= y .- y_old
        err = max(norm(dx), norm(dy))

        isfinite(err) || error("EllipticSolver diverged at iteration $iter (non-finite update); try a smaller ω or weaker wall forcing")

        if iter % 500 == 0
            verbose_print("Iter: $iter, Error: $err")
        end
        if err < params.tol
            verbose_print("Convergence reached at iteration $iter with error $err.")
            break
        end
        if iter == params.max_iter
            @warn "EllipticSolver reached max_iter=$(params.max_iter) without converging (error $err > tol $(params.tol))"
        end
    end

    return x, y, err, finalIter
end
