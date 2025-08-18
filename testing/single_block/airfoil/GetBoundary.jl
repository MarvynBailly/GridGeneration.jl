function getBoundaryConditions(grid; trailingEdge = 0)
    boundaryConditions = Dict(
        "BCInflow" => [],
        "BCOutflow" => [],
        "BCWall" => [],
        "BCInterface" => []
    )

    blk = 1
    ni = size(grid, 2)
    nj = size(grid, 3)
    nk = 1

    # left Outflow
    push!(boundaryConditions["BCOutflow"], Dict(
        "block" => blk,
        "start" => [1, 1, 1],
        "end" => [1, nj, nk]
    ))

    #right Outflow
    push!(boundaryConditions["BCOutflow"], Dict(
        "block" => blk,
        "start" => [ni, 1, 1],
        "end" => [ni, nj, nk]
    ))

    # top inflow
    push!(boundaryConditions["BCInflow"], Dict(
        "block" => blk,
        "start" => [1, nj, 1],
        "end" => [ni, nj, nk]
    ))

    if trailingEdge == 0
        push!(boundaryConditions["BCWall"], Dict(
            "block" => blk,
            "start" => [1, 1, 1],
            "end" => [ni, 1, nk]
        ))
    else
        # bottom interface
        push!(boundaryConditions["BCInterface"], Dict(
            "block" => blk,
            "start" => [1, 1, 1],
            "end" => [horzN+1, 1, nk]
        ))
        # bottom wall
        push!(boundaryConditions["BCWall"], Dict(
            "block" => blk,
            "start" => [horzN+1, 1, 1],
            "end" => [horzN+airfoilN, 1, nk]
        ))

        push!(boundaryConditions["BCInterface"], Dict(
            "block" => blk,
            "start" => [horzN+airfoilN, 1, 1],
            "end" => [ni, 1, nk]
        ))
    end
    bndInfo = []
    for (name, faces) in boundaryConditions
        if !isempty(faces)
            push!(bndInfo, Dict("name"=>name, "faces"=>faces))
        end
    end
    return bndInfo
end