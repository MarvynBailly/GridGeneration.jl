include("C:\\Users\\admin\\Documents\\GitHub\\GridGeneration\\src\\GridGeneration.jl")
using .GridGeneration

include("setup/domain_setup.jl")

println("\n" * "="^70)
println("TESTING UNEVEN INTERFACE FIX")
println("="^70)

blocks, bndInfo, interfaceInfo, M = setup_turtle_grid_domain(
    "step/BFstepTest_entropy.metric", 
    "step/coarseGrids";
    fix_uneven_interfaces=true
)

println("\n" * "-"^70)
println("BLOCKS")
println("-"^70)
println("Number of blocks: ", length(blocks))
for i in 1:length(blocks)
    println("  Block $i size: (i=$(size(blocks[i], 2)), j=$(size(blocks[i], 3)))")
end

println("\n" * "-"^70)
println("INTERFACES")
println("-"^70)
println("Number of interfaces: ", length(interfaceInfo))
for (i, itf) in enumerate(interfaceInfo)
    println("\nInterface $i:")
    println("  Connects: Block $(itf["blockA"]) <-> Block $(itf["blockB"])")
    println("  Block $(itf["blockA"]): start=$(itf["start_blkA"]), end=$(itf["end_blkA"])")
    println("  Block $(itf["blockB"]): start=$(itf["start_blkB"]), end=$(itf["end_blkB"])")
end

println("\n" * "-"^70)
println("BOUNDARIES")
println("-"^70)
println("Number of boundary types: ", length(bndInfo))
for bnd in bndInfo
    println("\nBoundary '$(bnd["name"])':")
    println("  Number of faces: $(length(bnd["faces"]))")
    for face in bnd["faces"]
        println("    Block $(face["block"]): start=$(face["start"]), end=$(face["end"])")
    end
end

println("\n" * "="^70)


