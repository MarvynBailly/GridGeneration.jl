"""
Get all (transitively) neighboring block IDs reachable from `blockId`
through interfaces aligned with `dir` (1=horizontal, 2=vertical).
- `interInfo` entries must have keys: "blockA","blockB","start_blkA","end_blkA".
- Set `include_start=true` to include `blockId` in the result.
"""
function GetNeighbors(blockId::Int, interInfo, dir::Int; include_start::Bool=false)
    seen = Set{Int}()

    function visit(bid::Int)
        # already visited this block
        bid ∈ seen && return
        push!(seen, bid)

        for info in interInfo
            if info["blockA"] == bid || info["blockB"] == bid
                # Orientation is defined by the A-side span (shared for both sides)
                startA = info["start_blkA"]
                endA   = info["end_blkA"]
                is_vertical   = (startA[1] != endA[1])
                is_horizontal = (startA[2] != endA[2])

                if (dir == 2 && is_vertical) || (dir == 1 && is_horizontal)
                    other = (info["blockA"] == bid) ? info["blockB"] : info["blockA"]
                    visit(other)
                end
            end
        end
    end

    visit(blockId)

    return include_start ? collect(seen) : [b for b in seen if b != blockId]
end