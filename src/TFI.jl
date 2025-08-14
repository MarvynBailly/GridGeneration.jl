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

    block = permutedims(cat(X', Y'; dims=3), (3, 1, 2))

    return block
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


# ---------------------------------------------
# Hermite TFI (
# ---------------------------------------------

# Cubic Hermite basis (scalars)
@inline h00(t) =  2t^3 - 3t^2 + 1
@inline h10(t) =    t^3 - 2t^2 + t
@inline h01(t) = -2t^3 + 3t^2
@inline h11(t) =    t^3 -   t^2

# 2D Hermite between endpoints P0,P1 with end-derivatives V0,V1 (all as 2-tuples)
@inline function hermite2(P0::Tuple, P1::Tuple, V0::Tuple, V1::Tuple, t)
    x = h00(t)*P0[1] + h01(t)*P1[1] + h10(t)*V0[1] + h11(t)*V1[1]
    y = h00(t)*P0[2] + h01(t)*P1[2] + h10(t)*V0[2] + h11(t)*V1[2]
    return (x, y)
end

# rotate (vx,vy) 90° CCW
@inline rot90ccw(v::Tuple) = (-v[2], v[1])

# normalize a 2-tuple; returns (0,0) if tiny
@inline function normalize2(v::Tuple)
    n = hypot(v[1], v[2])
    n > 0 ? (v[1]/n, v[2]/n) : (0.0, 0.0)
end

# centered/one-sided tangent along an N×2 edge at row i, returned as 2-tuple unit vector
function edge_tangent(edge::AbstractMatrix{<:Real}, i::Int)
    N = size(edge, 1)
    if i == 1
        tx = edge[2,1] - edge[1,1]; ty = edge[2,2] - edge[1,2]
    elseif i == N
        tx = edge[N,1] - edge[N-1,1]; ty = edge[N,2] - edge[N-1,2]
    else
        tx = edge[i+1,1] - edge[i-1,1]; ty = edge[i+1,2] - edge[i-1,2]
    end
    return normalize2((tx, ty))
end

# local edge step length Δs at node i
function edge_step(edge::AbstractMatrix{<:Real}, i::Int)
    N = size(edge, 1)
    if i == 1
        return hypot(edge[2,1]-edge[1,1], edge[2,2]-edge[1,2])
    elseif i == N
        return hypot(edge[N,1]-edge[N-1,1], edge[N,2]-edge[N-1,2])
    else
        f = hypot(edge[i+1,1]-edge[i,1], edge[i+1,2]-edge[i,2])
        b = hypot(edge[i,1]-edge[i-1,1], edge[i,2]-edge[i-1,2])
        return 0.5*(f + b)
    end
end

# inward derivatives along an edge; returns N×2 matrix (∂X/∂η at bottom/top, ∂X/∂ξ at left/right)
# where ∈ (:bottom, :top, :left, :right)
function inward_derivatives(edge::AbstractMatrix{<:Real}, where::Symbol;
                            Δcomp::Real, damping::Real=1.0, αcap::Union{Nothing,Real}=nothing)
    N = size(edge, 1)
    d = zeros(eltype(edge), N, 2)
    # sign mapping so rot90ccw(tangent) points inward
    sgn = where === :bottom ? +1 :
          where === :top    ? -1 :
          where === :left   ? -1 : +1   # right => +1

    # (optional) cap to prevent huge corner jets
    # compute a robust global scale if capping
    medΔs = αcap === nothing ? 0.0 : median([edge_step(edge,i) for i in 1:N])

    for i in 1:N
        t  = edge_tangent(edge, i)
        n  = normalize2((sgn*rot90ccw(t)[1], sgn*rot90ccw(t)[2]))
        Δs = edge_step(edge, i)
        α  = damping * (Δs / Δcomp)
        if αcap !== nothing
            α = min(α, αcap * (medΔs / Δcomp))
        end
        d[i,1] = α * n[1]
        d[i,2] = α * n[2]
    end
    return d
end

"""
    TFI_2D_Hermite(boundary; db=nothing, dt=nothing, dl=nothing, dr=nothing, damping=1.0, αcap=nothing)

Hermite (C¹) transfinite interpolation that carries boundary spacing into the interior.

Inputs:
- `boundary = [top, right, bottom, left]`
    top::(N1×2), right::(N2×2), bottom::(N1×2), left::(N2×2)
- Optional inward derivatives:
    db (N1×2) = ∂X/∂η at bottom, dt (N1×2) = ∂X/∂η at top,
    dl (N2×2) = ∂X/∂ξ at left,   dr (N2×2) = ∂X/∂ξ at right.
  If omitted, they’re auto-built from edge normals with magnitude ≈ (Δs / Δcomp).
- `damping` ∈ (0,1] scales the inward derivatives.
- `αcap` (optional) caps derivative magnitudes to `αcap * median(Δs)/Δcomp` to avoid blow-ups at corners.

Returns:
- `block::Array{T,3}` with shape (2, N1, N2).
"""
function TFI_2D_Hermite(boundary; db=nothing, dt=nothing, dl=nothing, dr=nothing, damping=1.0, αcap=nothing)
    top, right, bottom, left = boundary
    @assert size(top,2)==2 && size(bottom,2)==2 && size(left,2)==2 && size(right,2)==2

    N1 = size(bottom, 1)  # along ξ (left→right)
    N2 = size(left,   1)  # along η (bottom→top)

    ξs = range(0.0, 1.0; length=N1)
    ηs = range(0.0, 1.0; length=N2)

    # Auto inward derivatives if not provided
    Δξ = 1.0 / (N1 - 1)
    Δη = 1.0 / (N2 - 1)
    db === nothing && (db = inward_derivatives(bottom, :bottom; Δcomp=Δη, damping=damping, αcap=αcap))
    dt === nothing && (dt = inward_derivatives(top,    :top;    Δcomp=Δη, damping=damping, αcap=αcap))
    dl === nothing && (dl = inward_derivatives(left,   :left;   Δcomp=Δξ, damping=damping, αcap=αcap))
    dr === nothing && (dr = inward_derivatives(right,  :right;  Δcomp=Δξ, damping=damping, αcap=αcap))

    # Corners (positions)
    P00 = (left[1,1],   left[1,2])     # ξ=0, η=0
    P10 = (right[1,1],  right[1,2])    # ξ=1, η=0
    P01 = (left[end,1], left[end,2])   # ξ=0, η=1
    P11 = (right[end,1],right[end,2])  # ξ=1, η=1

    # Corner jets (take first/last entries on each edge derivative)
    Bη0 = (db[1,1],   db[1,2]);  Bη1 = (db[end,1], db[end,2])
    Tη0 = (dt[1,1],   dt[1,2]);  Tη1 = (dt[end,1], dt[end,2])
    Lξ0 = (dl[1,1],   dl[1,2]);  Lξ1 = (dl[end,1], dl[end,2])
    Rξ0 = (dr[1,1],   dr[1,2]);  Rξ1 = (dr[end,1], dr[end,2])

    X = zeros(eltype(bottom), N2, N1)
    Y = zeros(eltype(bottom), N2, N1)

    for j in 1:N2
        η = ηs[j]
        Lp = (left[j,1],  left[j,2])
        Rp = (right[j,1], right[j,2])
        Ld = (dl[j,1],    dl[j,2])
        Rd = (dr[j,1],    dr[j,2])

        for i in 1:N1
            ξ  = ξs[i]
            Bp = (bottom[i,1], bottom[i,2])
            Tp = (top[i,1],    top[i,2])
            Bd = (db[i,1],     db[i,2])
            Td = (dt[i,1],     dt[i,2])

            # Two Hermite strips
            Sb = hermite2(Bp, Tp, Bd, Td, η)   # vary η at fixed ξ
            Sl = hermite2(Lp, Rp, Ld, Rd, ξ)   # vary ξ at fixed η

            # Corner correction (bi-Hermite "Coons-like" patch)
            Cξ0 = hermite2(P00, P01, Bη0, Tη0, η)
            Cξ1 = hermite2(P10, P11, Bη1, Tη1, η)
            Dη0 = hermite2(P00, P10, Lξ0, Rξ0, ξ)
            Dη1 = hermite2(P01, P11, Lξ1, Rξ1, ξ)
            Hη  = hermite2(Cξ0, Cξ1, (0.0,0.0), (0.0,0.0), ξ)
            Hξ  = hermite2(Dη0, Dη1, (0.0,0.0), (0.0,0.0), η)

            BLx = (1-ξ)*(1-η)*P00[1] + ξ*(1-η)*P10[1] + (1-ξ)*η*P01[1] + ξ*η*P11[1]
            BLy = (1-ξ)*(1-η)*P00[2] + ξ*(1-η)*P10[2] + (1-ξ)*η*P01[2] + ξ*η*P11[2]

            Cx = Hη[1] + Hξ[1] - BLx
            Cy = Hη[2] + Hξ[2] - BLy

            Px = Sb[1] + Sl[1] - Cx
            Py = Sb[2] + Sl[2] - Cy

            X[j,i] = Px;  Y[j,i] = Py
        end
    end

    # Return layout: (2, N1, N2) like your TFI_2D
    return permutedims(cat(X', Y'; dims=3), (3, 1, 2))
end
