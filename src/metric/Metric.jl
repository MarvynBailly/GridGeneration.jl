using NearestNeighbors

struct PtRef
    block::Int
    i::Int
    j::Int
end

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


function find_nearest_kd(data, tree::KDTree, refs, xq, yq)
    idxs, dists = knn(tree, [xq,yq], 1)   # 1‐NN
    ref = refs[idxs[1]]
    blk = ref.block
    i,j = ref.i, ref.j
    M11_val = data["M11"][blk][i, j]
    M22_val = data["M22"][blk][i, j]
    return [M11_val, M22_val]
end