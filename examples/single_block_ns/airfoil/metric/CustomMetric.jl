using LinearAlgebra

# Smooth, positive, decaying profiles
@inline w_rational(d, A, ℓ, p) = A / (1 + (d/ℓ)^p)
@inline w_gauss(d, A, ℓ)      = A * exp(-(d/ℓ)^2)

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
function make_getMetric(boundary; closed=true,
                        A_airfoil=1.0, ℓ_airfoil=0.05, p_airfoil=2,
                        A_origin=1.0,  ℓ_origin=0.10, p_origin=2,
                        floor=1e-4, profile=:rational)

    bx = @view boundary[1, :]
    by = @view boundary[2, :]
    N  = length(bx)
    @assert N ≥ 2 "boundary must have at least 2 points"

    last = closed ? N : N-1
    x1 = bx[1:last]
    y1 = by[1:last]
    x2 = [bx[2:last]; closed ? bx[1] : bx[end]]
    y2 = [by[2:last]; closed ? by[1] : by[end]]

    weight_airfoil(d) = profile === :gauss ? w_gauss(d, A_airfoil, ℓ_airfoil) : w_rational(d, A_airfoil, ℓ_airfoil, p_airfoil)
    weight_origin(d)  = profile === :gauss ? w_gauss(d, A_origin,  ℓ_origin) : w_rational(d, A_origin,  ℓ_origin,  p_origin)

    # point-to-polyline (segment) distance
    function d_to_polyline(x::Real, y::Real)
        d2min = Inf
        @inbounds for i in eachindex(x1)
            dx = x2[i] - x1[i]; dy = y2[i] - y1[i]
            denom = dx*dx + dy*dy + eps()   # handle possible zero-length segment
            t = ((x - x1[i])*dx + (y - y1[i])*dy) / denom
            t = clamp(t, 0.0, 1.0)
            px = x1[i] + t*dx; py = y1[i] + t*dy
            d2 = (x - px)^2 + (y - py)^2
            if d2 < d2min; d2min = d2; end
        end
        return sqrt(d2min)
    end

    return function getMetric(x::Real, y::Real)
        d_airfoil = d_to_polyline(x, y)  # distance to airfoil
        d_origin  = hypot(x, y)          # distance to (0,0)
        w = floor + weight_airfoil(d_airfoil) + weight_origin(d_origin)
        return (w, w)  # isotropic
    end
end
