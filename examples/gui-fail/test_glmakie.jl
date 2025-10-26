"""
Simple test to verify GLMakie is working correctly

Run this first to ensure GLMakie can display windows before trying the full GUI.

Usage:
    julia test_glmakie.jl
"""

println("Testing GLMakie installation and display...")

try
    using GLMakie
    println("✅ GLMakie loaded successfully")
catch e
    println("❌ Error loading GLMakie:")
    println(e)
    println("\nInstall with: using Pkg; Pkg.add(\"GLMakie\")")
    exit(1)
end

println("Creating a simple test plot...")

try
    # Create a simple figure
    fig = Figure(resolution=(800, 600))
    ax = Axis(fig[1, 1], title="GLMakie Test Plot")
    
    # Plot a simple line
    x = range(0, 10, length=100)
    y = sin.(x)
    lines!(ax, x, y, color=:blue, linewidth=2)
    
    println("✅ Figure created successfully")
    println("📊 Displaying figure...")
    
    screen = display(fig)
    
    println("\n" * "="^60)
    println("SUCCESS! GLMakie window should now be visible.")
    println("="^60)
    println("\nIf you can see a sine wave plot:")
    println("  ✅ GLMakie is working correctly")
    println("  ✅ You can proceed to use GridGenerationGUI.jl")
    println("\nIf you CANNOT see the window:")
    println("  ⚠️  Check that you're not in a headless environment")
    println("  ⚠️  Update graphics drivers")
    println("  ⚠️  Try the Web GUI instead (GridGenerationWebGUI.jl)")
    println("\nClose the window or press Ctrl+C to exit...")
    
    # Keep window open
    try
        while isopen(screen)
            sleep(0.1)
        end
        println("\n✅ Window closed successfully")
    catch e
        if isa(e, InterruptException)
            println("\n✅ Interrupted by user")
        else
            rethrow(e)
        end
    end
    
catch e
    println("❌ Error creating or displaying figure:")
    println(e)
    println("\nThis might indicate:")
    println("  - Missing OpenGL support")
    println("  - Graphics driver issues")
    println("  - Headless environment (no display)")
    println("\nTry the Web GUI instead: julia GridGenerationWebGUI.jl")
    exit(1)
end
