using NearestNeighbors

struct PtRef
    block::Int
    i::Int
    j::Int
end

"""
    setup_metric_tree(data) -> (tree, refs)
    setup_metric_tree(blocks::Vector{Array{Float64,3}}) -> (tree, refs)

Build a KD-tree over the points of a gridded metric field for nearest-neighbour lookup.
`data` is a dictionary with per-block coordinate matrices `data["x"]` and `data["y"]`;
alternatively pass the blocks directly as `[2, Ni, Nj]` arrays.
`refs[k]` records the (block, i, j) index of the k-th point in the tree.
Use with [`find_nearest_kd`](@ref).
"""
function setup_metric_tree(data)
    refs = PtRef[]
    coords = Float64[]
    for (b, (Xb,Yb)) in enumerate(zip(data["x"], data["y"]))
    Ny, Nz = size(Xb)
    for i in 1:Ny, j in 1:Nz
        push!(coords, Xb[i,j])
        push!(coords, Yb[i,j])
        push!(refs, PtRef(b,i,j))
    end
    end
    coords = reshape(coords, (2, length(refs)))  # now 2×N_points
    
    tree = KDTree(coords)
    return tree, refs
end

function setup_metric_tree(blocks::Array{Array{Float64,3},1})
    refs = PtRef[]
    coords = Float64[]
    for (b, blk) in enumerate(blocks)
        Xb, Yb = blk[1,:,:], blk[2,:,:]
        Ny, Nz = size(Xb)
        for i in 1:Ny, j in 1:Nz
            push!(coords, Xb[i,j])
            push!(coords, Yb[i,j])
            push!(refs, PtRef(b,i,j))
        end
    end
    coords = reshape(coords, (2, length(refs)))  # now 2×N_points
    
    tree = KDTree(coords)
    return tree, refs
end


"""
    find_nearest_kd(data, tree, refs, xq, yq) -> [M11, M22]

Return the diagonal metric components `[M11, M22]` stored in `data["M11"]`/`data["M22"]`
at the grid point nearest to `(xq, yq)`. `tree` and `refs` come from [`setup_metric_tree`](@ref).
"""
function find_nearest_kd(data, tree::KDTree, refs, xq, yq)
    idxs, dists = knn(tree, [xq,yq], 1)   # 1‐NN
    ref = refs[idxs[1]]
    blk = ref.block
    i,j = ref.i, ref.j

    if length(size(data["M11"][1])) == 3
        M11_val = data["M11"][blk][i, j, 1]
        M22_val = data["M22"][blk][i, j, 1]
        return [M11_val, M22_val]
    elseif length(size(data["M11"][1])) == 2
        M11_val = data["M11"][blk][i, j]
        M22_val = data["M22"][blk][i, j]
        return [M11_val, M22_val]
    else
        error("Unsupported metric data format")
    end
end
