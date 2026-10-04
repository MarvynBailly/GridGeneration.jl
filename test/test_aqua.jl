using Test
using Aqua
using GridGeneration

@testset "Aqua quality checks" begin
    Aqua.test_all(GridGeneration)
end
