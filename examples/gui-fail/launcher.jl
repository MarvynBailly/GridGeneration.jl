"""
GridGeneration.jl GUI Launcher

Simple script to launch either GUI version with dependency checking.

Usage:
    julia launcher.jl [gui_type]
    
Where gui_type is:
    - makie (default): Launch GLMakie GUI
    - web: Launch web-based GUI
    - both: Show selection menu
"""

using Pkg

function check_dependencies(deps::Vector{String})
    """Check if required packages are installed"""
    missing_deps = String[]
    
    for dep in deps
        try
            eval(Meta.parse("using $dep"))
        catch
            push!(missing_deps, dep)
        end
    end
    
    return missing_deps
end

function install_dependencies(deps::Vector{String})
    """Install missing dependencies"""
    println("\n📦 Installing missing dependencies...")
    for dep in deps
        println("  Installing $dep...")
        try
            Pkg.add(dep)
        catch e
            @warn "Failed to install $dep" exception=e
        end
    end
end

function launch_makie_gui()
    """Launch GLMakie-based GUI"""
    println("\n🚀 Launching GLMakie GUI...")
    
    required = ["GLMakie", "DelimitedFiles"]
    missing = check_dependencies(required)
    
    if !isempty(missing)
        println("\n⚠️  Missing dependencies: $(join(missing, ", "))")
        print("Install now? (y/n): ")
        response = strip(readline())
        
        if lowercase(response) == "y"
            install_dependencies(missing)
        else
            println("❌ Cannot launch without dependencies")
            return
        end
    end
    
    # Check for MAT.jl
    try
        eval(Meta.parse("using MAT"))
    catch
        println("\n📦 Installing MAT.jl...")
        Pkg.add(url="https://github.com/JuliaIO/MAT.jl")
    end
    
    println("\n✅ All dependencies satisfied")
    println("🎨 Opening GLMakie window...")
    
    include("GridGenerationGUI.jl")
end

function launch_web_gui()
    """Launch web-based GUI"""
    println("\n🌐 Launching Web GUI...")
    
    required = ["Interact", "Blink", "Plots", "WebIO", "DelimitedFiles"]
    missing = check_dependencies(required)
    
    if !isempty(missing)
        println("\n⚠️  Missing dependencies: $(join(missing, ", "))")
        print("Install now? (y/n): ")
        response = strip(readline())
        
        if lowercase(response) == "y"
            install_dependencies(missing)
        else
            println("❌ Cannot launch without dependencies")
            return
        end
    end
    
    # Check for MAT.jl
    try
        eval(Meta.parse("using MAT"))
    catch
        println("\n📦 Installing MAT.jl...")
        Pkg.add(url="https://github.com/JuliaIO/MAT.jl")
    end
    
    # Ensure Blink has Electron
    println("\n🔧 Checking Blink.jl Electron installation...")
    try
        # Check if Blink is installed and has Electron
        if "Blink" in [pkg.name for pkg in values(Pkg.dependencies())]
            eval(Meta.parse("using Blink"))
            # Now check for Electron using eval
            has_electron = eval(Meta.parse("Blink.AtomShell.isinstalled()"))
            if !has_electron
                println("📦 Installing Electron for Blink...")
                eval(Meta.parse("Blink.AtomShell.install()"))
            end
        end
    catch e
        @warn "Issue with Blink setup" exception=e
    end
    
    println("\n✅ All dependencies satisfied")
    println("🌐 Opening browser window...")
    
    include("GridGenerationWebGUI.jl")
end

function show_menu()
    """Interactive menu for GUI selection"""
    println("\n" * "="^60)
    println("  GridGeneration.jl - Interactive GUI Launcher")
    println("="^60)
    println("\nSelect GUI version:")
    println("  1. GLMakie GUI (recommended for desktop)")
    println("     - Native OpenGL rendering")
    println("     - Fast, responsive")
    println("     - Best for large grids")
    println()
    println("  2. Web GUI (recommended for simplicity)")
    println("     - Browser-based interface")
    println("     - Easier setup")
    println("     - Familiar controls")
    println()
    print("Enter choice (1 or 2): ")
    
    choice = strip(readline())
    
    if choice == "1"
        launch_makie_gui()
    elseif choice == "2"
        launch_web_gui()
    else
        println("❌ Invalid choice. Please enter 1 or 2.")
        show_menu()
    end
end

function main()
    """Main launcher function"""
    
    # Check if we're in the right directory
    if !isfile("GridGenerationGUI.jl") || !isfile("GridGenerationWebGUI.jl")
        println("❌ Error: GUI files not found")
        println("   Please run this script from the examples/gui/ directory")
        return
    end
    
    # Parse command line argument
    gui_type = length(ARGS) > 0 ? lowercase(ARGS[1]) : "both"
    
    if gui_type == "makie"
        launch_makie_gui()
    elseif gui_type == "web"
        launch_web_gui()
    elseif gui_type in ["both", "menu", "select"]
        show_menu()
    else
        println("❌ Unknown GUI type: $gui_type")
        println("\nUsage: julia launcher.jl [gui_type]")
        println("  gui_type: makie, web, or both (default)")
    end
end

# Run launcher
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
