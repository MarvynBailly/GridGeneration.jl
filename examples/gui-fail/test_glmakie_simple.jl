"""
Simple GLMakie window test - Alternative method

This uses GLMakie's screen management more explicitly.
"""

using GLMakie

println("Testing GLMakie with explicit screen management...")
println()

# Create a simple figure
fig = Figure(resolution=(800, 600))
ax = Axis(fig[1, 1], title="GLMakie Test - If you see this, it works!")

# Plot something
x = range(0, 10, length=100)
lines!(ax, x, sin.(x), color=:blue, linewidth=3, label="sin(x)")
lines!(ax, x, cos.(x), color=:red, linewidth=3, label="cos(x)")
axislegend(ax)

println("Displaying window...")
println()

# Method 1: Get screen and wait for it
screen = display(fig)

println("="^60)
println("  Window should now be visible!")
println("="^60)
println()
println("What you should see:")
println("  - A window titled 'GLMakie...'")
println("  - Blue and red sine/cosine waves")
println("  - Legend in top-right")
println()
println("To close: Close the window or press Ctrl+C")
println("="^60)
println()

# Wait for window to close
try
    while isopen(screen)
        sleep(0.1)
    end
    println("\n✅ Success! Window was displayed and closed properly.")
catch e
    if isa(e, InterruptException)
        println("\n✅ Interrupted - GLMakie is working!")
    else
        println("\n❌ Error: $e")
        rethrow(e)
    end
end
