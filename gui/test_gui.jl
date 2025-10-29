# Simple test script to verify GUI loads
println("Loading GUI...")
try
    include("GridGenerationGUI.jl")
    println("GUI loaded successfully!")
    println("Waiting 3 seconds...")
    sleep(3)
    println("Test complete.")
catch e
    println("ERROR: Failed to load GUI")
    println(e)
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
end
