using LinearAlgebra

# ---------- helpers: projection to segments / polyline (with closest point) ----------
@inline function project_point_to_segment(x, y, x1, y1, x2, y2)
    vx, vy = x2 - x1, y2 - y1
    wx, wy = x  - x1, y  - y1
    c2 = vx*vx + vy*vy
    if c2 == 0.0
        # degenerate segment; project to endpoint 1
        return x1, y1, 0.0
    end
    t = (vx*wx + vy*wy) / c2
    t_clamped = clamp(t, 0.0, 1.0)
    px, py = x1 + t_clamped*vx, y1 + t_clamped*vy
    return px, py, t_clamped
end

"""
Return (d, px, py) where (px,py) is the closest point on the polyline to (x,y).
`airfoil` can be 2×N or N×2. If `closed=true` we also check the closing segment.
"""
function dist_to_polyline_with_projection(x, y, airfoil; closed::Bool=true)
    if size(airfoil,1)==2
        N = size(airfoil,2)
        best_d2 = Inf; px=NaN; py=NaN
        for i in 1:(N-1)
            pxi, pyi, _ = project_point_to_segment(x, y, airfoil[1,i], airfoil[2,i], airfoil[1,i+1], airfoil[2,i+1])
            dx = x - pxi; dy = y - pyi; d2 = dx*dx + dy*dy
            if d2 < best_d2
                best_d2 = d2; px = pxi; py = pyi
            end
        end
        if closed
            pxi, pyi, _ = project_point_to_segment(x, y, airfoil[1,N], airfoil[2,N], airfoil[1,1], airfoil[2,1])
            dx = x - pxi; dy = y - pyi; d2 = dx*dx + dy*dy
            if d2 < best_d2
                best_d2 = d2; px = pxi; py = pyi
            end
        end
        return sqrt(best_d2), px, py
    elseif size(airfoil,2)==2
        N = size(airfoil,1)
        best_d2 = Inf; px=NaN; py=NaN
        for i in 1:(N-1)
            pxi, pyi, _ = project_point_to_segment(x, y, airfoil[i,1], airfoil[i,2], airfoil[i+1,1], airfoil[i+1,2])
            dx = x - pxi; dy = y - pyi; d2 = dx*dx + dy*dy
            if d2 < best_d2
                best_d2 = d2; px = pxi; py = pyi
            end
        end
        if closed
            pxi, pyi, _ = project_point_to_segment(x, y, airfoil[N,1], airfoil[N,2], airfoil[1,1], airfoil[1,2])
            dx = x - pxi; dy = y - pyi; d2 = dx*dx + dy*dy
            if d2 < best_d2
                best_d2 = d2; px = pxi; py = pyi
            end
        end
        return sqrt(best_d2), px, py
    else
        error("airfoil must be 2×N or N×2")
    end
end

# ---------- metric with gradient ----------
"""
    make_getMetric_with_gradient(airfoil; kwargs...)

Build a closure `(x,y) -> (M_vec, gradM)` where:
- `M_vec = [M, M]` (isotropic metric, as in your code)
- `gradM = (dMdx, dMdy)`

Keyword args mirror `make_getMetric`:
- A_airfoil, ℓ_airfoil, p_airfoil
- A_origin,  ℓ_origin,  p_origin, origin_center
- floor, profile∈{:rational,:gauss}, closed_airfoil
"""
function make_getMetric_with_gradient(airfoil;
    A_airfoil::Real = 400.0,  ℓ_airfoil::Real = 0.5, p_airfoil::Real = 2,
    A_origin::Real  = 10000.0, ℓ_origin::Real  = 0.05, p_origin::Real  = 10,
    origin_center::Tuple{<:Real,<:Real} = (0.0, 0.0),
    floor::Real = 1e-4,
    profile::Symbol = :rational,
    closed_airfoil::Bool = true,
)

    ax, ay = origin_center

    # radial profile and its radial derivative dw/dd (both support :rational and :gauss with general p)
    @inline function w_and_dwdd(d::Real, A::Real, ℓ::Real, p::Real)
        if A == 0
            return (0.0, 0.0)
        end
        if profile === :rational
            u = (d/ℓ)^p
            w = A / (1 + u)
            # dw/dd = -A * p/ℓ^p * d^(p-1) / (1+u)^2
            dwdd = -A * p / (ℓ^p) * (d^(p-1)) / (1 + u)^2
            return (w, dwdd)
        elseif profile === :gauss
            u = (d/ℓ)^p
            e = exp(-u)
            w = A * e
            # dw/dd = A * e * (- p/ℓ^p * d^(p-1))
            dwdd = -A * p / (ℓ^p) * (d^(p-1)) * e
            return (w, dwdd)
        else
            error("profile must be :rational or :gauss")
        end
    end

    return function (x::Real, y::Real)
        # ----- distances & their gradients -----
        # airfoil polyline
        d_air, px, py = dist_to_polyline_with_projection(x, y, airfoil; closed=closed_airfoil)
        if d_air > 0
            gx_air = (x - px)/d_air
            gy_air = (y - py)/d_air
        else
            gx_air = 0.0; gy_air = 0.0   # undefined at exact projection; choose 0
        end

        # hotspot at (ax, ay)
        dxh = x - ax; dyh = y - ay
        d_hot = hypot(dxh, dyh)
        if d_hot > 0
            gx_hot = dxh/d_hot
            gy_hot = dyh/d_hot
        else
            gx_hot = 0.0; gy_hot = 0.0
        end

        # ----- contributions and chain rule -----
        w_air, dwdd_air = w_and_dwdd(d_air, A_airfoil, ℓ_airfoil, p_airfoil)
        w_hot, dwdd_hot = w_and_dwdd(d_hot, A_origin,  ℓ_origin,  p_origin)

        M = floor + w_air + w_hot
        dMdx = dwdd_air * gx_air + dwdd_hot * gx_hot
        dMdy = dwdd_air * gy_air + dwdd_hot * gy_hot

        return ([M, M], (dMdx, dMdy))
    end
end


using LinearAlgebra

# Smooth, positive, decaying profiles
@inline w_rational(d, A, ℓ, p) = A / (1 + (d/ℓ)^p)
@inline w_gauss(d, A, ℓ)      = A * exp(-(d/ℓ)^2)


# --- geometry helpers (handle airfoil as 2×N or N×2) ---
@inline function dist_point_segment(x, y, x1, y1, x2, y2)
    vx, vy = x2 - x1, y2 - y1
    wx, wy = x  - x1, y  - y1
    c1 = vx*wx + vy*wy
    if c1 <= 0
        return hypot(wx, wy)
    end
    c2 = vx*vx + vy*vy
    if c2 <= c1
        return hypot(x - x2, y - y2)
    end
    t = c1 / c2
    px, py = x1 + t*vx, y1 + t*vy
    return hypot(x - px, y - py)
end

function dist_to_polyline(x, y, airfoil::AbstractMatrix; closed::Bool=true)
    if size(airfoil, 1) == 2         # 2×N
        N = size(airfoil, 2)
        dmin = Inf
        for i in 1:(N-1)
            dmin = min(dmin, dist_point_segment(x, y,
                airfoil[1,i], airfoil[2,i], airfoil[1,i+1], airfoil[2,i+1]))
        end
        if closed
            dmin = min(dmin, dist_point_segment(x, y,
                airfoil[1,N], airfoil[2,N], airfoil[1,1], airfoil[2,1]))
        end
        return dmin
    elseif size(airfoil, 2) == 2     # N×2
        N = size(airfoil, 1)
        dmin = Inf
        for i in 1:(N-1)
            dmin = min(dmin, dist_point_segment(x, y,
                airfoil[i,1], airfoil[i,2], airfoil[i+1,1], airfoil[i+1,2]))
        end
        if closed
            dmin = min(dmin, dist_point_segment(x, y,
                airfoil[N,1], airfoil[N,2], airfoil[1,1], airfoil[1,2]))
        end
        return dmin
    else
        error("airfoil must be 2×N or N×2")
    end
end


# --- metric factory with movable hotspot ---
"""
    make_getMetric(boundary; closed=true,
                   A_airfoil=1.0, ℓ_airfoil=0.05, p_airfoil=2,
                   A_origin=1.0,  ℓ_origin=0.10, p_origin=2,
                   floor=1e-4, profile=:rational)

Create a closure `getMetric(x,y) => (M11, M22)` where the scalar field
w(x,y) = floor + w_airfoil(dist_to_airfoil) + w_origin(dist_to_origin)
is applied isotropically: M11 = M22 = w.

- `boundary` is a 2×N polyline for the airfoil (assumed closed by default).
- `A_*` set the peak amplitudes.
- `ℓ_*` set decay lengths (units = your coordinate units).
- `p_*` control tail sharpness for the rational profile.
- `floor` prevents degeneracy far away.
- `profile=:rational` or `:gauss`.
"""
function make_getMetric(airfoil;
    A_airfoil::Real = 400.0,  ℓ_airfoil::Real = 0.5, p_airfoil::Real = 2,
    A_origin::Real  = 10000.0, ℓ_origin::Real  = 0.05, p_origin::Real  = 10,
    origin_center::Tuple{<:Real,<:Real} = (0.0, 0.0),
    floor::Real = 1e-4,
    profile::Symbol = :rational,
    closed_airfoil::Bool = true,
)
    ax, ay = origin_center
    @inline contrib(d, A, ℓ, p) = profile === :rational ? A / (1 + (d/ℓ)^p) :
                                                    profile === :gauss    ? A * exp(-(d/ℓ)^p) :
                                                    error("profile must be :rational or :gauss")


    return function (x::Real, y::Real)
        # distances
        d_air = dist_to_polyline(x, y, airfoil; closed=closed_airfoil)
        d_hot = hypot(x - ax, y - ay)                # distance to (a,b)

        # contributions
        c_air = A_airfoil == 0 ? 0.0 : contrib(d_air, A_airfoil, ℓ_airfoil, p_airfoil)
        c_hot = A_origin  == 0 ? 0.0 : contrib(d_hot, A_origin,  ℓ_origin,  p_origin)

        M = floor + c_air + c_hot
        return [M, M]   # isotropic metric: diag(M, M)
    end
end
