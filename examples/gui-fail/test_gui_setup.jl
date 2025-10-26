"""
Test script to verify GUI dependencies and basic functionality

Run this to check if your GUI setup is working before launching.

Usage:
    julia test_gui_setup.jl
"""

using Pkg

println("="^60)
println("  GridGeneration GUI Setup Test")
println("="^60)

# Test 1: Core GridGeneration
println("\n[Test 1/5] Checking GridGeneration.jl...")
try
    include("../../src/GridGeneration.jl")
    using .GridGeneration
    println("  ✅ GridGeneration.jl loads successfully")
    
    # Test basic exports
    @assert isdefined(GridGeneration, :GenerateGrid)
    @assert isdefined(GridGeneration, :TFI)
    @assert isdefined(GridGeneration, :SimParams)
    println("  ✅ Core functions available")
catch e
    println("  ❌ Error loading GridGeneration.jl")
    println("     $e")
    exit(1)
end

# Test 2: Example helpers
println("\n[Test 2/5] Checking example helpers...")
try
    include("../airfoil/airfoil.jl")
    include("../rectangle/rectangle.jl")
    println("  ✅ Airfoil helpers load")
    println("  ✅ Rectangle helpers load")
catch e
    println("  ❌ Error loading example helpers")
    println("     $e")
    exit(1)
end

# Test 3: GLMakie dependencies
println("\n[Test 3/5] Checking GLMakie GUI dependencies...")
glmakie_ok = true
for pkg in ["GLMakie", "DelimitedFiles", "Printf"]
    try
        eval(Meta.parse("using $pkg"))
        println("  ✅ $pkg available")
    catch
        println("  ⚠️  $pkg not installed")
        glmakie_ok = false
    end
end

try
    using MAT
    println("  ✅ MAT.jl available")
catch
    println("  ⚠️  MAT.jl not installed")
    glmakie_ok = false
end

if glmakie_ok
    println("  ✅ All GLMakie GUI dependencies satisfied")
else
    println("  ℹ️  Run: Pkg.add([\"GLMakie\", \"DelimitedFiles\", \"Printf\"])")
    println("         Pkg.add(url=\"https://github.com/JuliaIO/MAT.jl\")")
end

# Test 4: Web GUI dependencies
println("\n[Test 4/5] Checking Web GUI dependencies...")
web_ok = true
for pkg in ["Interact", "Blink", "Plots", "WebIO", "DelimitedFiles"]
    try
        eval(Meta.parse("using $pkg"))
        println("  ✅ $pkg available")
    catch
        println("  ⚠️  $pkg not installed")
        web_ok = false
    end
end

if web_ok
    println("  ✅ All Web GUI dependencies satisfied")
    
    # Check Blink Electron
    try
        using Blink
        if Blink.AtomShell.isinstalled()
            println("  ✅ Blink Electron installed")
        else
            println("  ⚠️  Blink Electron not installed")
            println("     Run: using Blink; Blink.AtomShell.install()")
            web_ok = false
        end
    catch e
        println("  ⚠️  Issue checking Blink: $e")
    end
else
    println("  ℹ️  Run: Pkg.add([\"Interact\", \"Blink\", \"Plots\", \"WebIO\"])")
end

# Test 5: Basic TFI functionality
println("\n[Test 5/5] Testing basic grid generation...")
try
    using .GridGeneration
    
    # Create simple square
    N = 5
    top = [range(0, 1, length=N) ones(N)]
    right = [ones(N) range(1, 0, length=N)]
    bottom = [range(1, 0, length=N) zeros(N)]
    left = [zeros(N) range(0, 1, length=N)]
    
    grid = TFI([top, right, bottom, left])
    
    @assert size(grid) == (2, N, N)
    println("  ✅ TFI grid generation works")
    
    # Test simple metric
    M(x, y) = [1.0, 1.0]
    params = SimParams(
        useSplitting=false,
        useEdgeSolver=false,
        useSmoothing=false
    )
    
    result = GenerateGrid(grid, [], [], M; params=params)
    println("  ✅ GenerateGrid works")
    
catch e
    println("  ❌ Error in basic grid generation")
    println("     $e")
    exit(1)
end

# Summary
println("\n" * "="^60)
println("  Summary")
println("="^60)

println("\nGUI Status:")
if glmakie_ok
    println("  ✅ GLMakie GUI ready to use")
    println("     Launch: julia GridGenerationGUI.jl")
else
    println("  ⚠️  GLMakie GUI needs dependencies")
    println("     Install and rerun this test")
end

if web_ok
    println("  ✅ Web GUI ready to use")
    println("     Launch: julia GridGenerationWebGUI.jl")
else
    println("  ⚠️  Web GUI needs dependencies")
    println("     Install and rerun this test")
end

if glmakie_ok || web_ok
    println("\n  🎉 You can use the launcher:")
    println("     julia launcher.jl")
else
    println("\n  ⚠️  Install dependencies first, then rerun this test")
end

println("\n" * "="^60)
println("  Test Complete!")
println("="^60)
