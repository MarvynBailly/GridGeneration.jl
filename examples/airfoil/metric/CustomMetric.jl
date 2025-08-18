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
