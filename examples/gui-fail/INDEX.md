# GridGeneration.jl GUI - Complete Guide

## 📁 Files in This Directory

| File | Purpose | Use When |
|------|---------|----------|
| `GridGenerationGUI.jl` | GLMakie-based native GUI | Desktop app, fast rendering |
| `GridGenerationWebGUI.jl` | Web-based GUI | Browser interface, simpler setup |
| `launcher.jl` | Interactive launcher script | Starting either GUI |
| `launch_gui.bat` | Windows batch launcher | Windows users (double-click) |
| `launch_gui.sh` | Unix/Mac shell launcher | Mac/Linux users |
| `test_gui_setup.jl` | Dependency checker | Before first run |
| `test_glmakie.jl` | GLMakie display test | Troubleshooting display issues |
| `README.md` | Main documentation | Complete reference |
| `QUICKSTART.md` | Fast tutorial | First-time users |
| `RUNNING.md` | Launch troubleshooting | When GUI won't appear |
| `ARCHITECTURE.md` | System design docs | Understanding internals |
| `TROUBLESHOOTING.md` | Problem solving | When things break |

## 🚀 Quick Start (30 seconds)

### Option 1: Interactive Launcher
```julia
cd("examples/gui")
julia launcher.jl
```

### Option 2: Direct Launch
```julia
# GLMakie GUI (recommended):
julia GridGenerationGUI.jl

# OR Web GUI:
julia GridGenerationWebGUI.jl
```

### Option 3: One-Click (Windows)
Double-click `launch_gui.bat`

### Option 4: One-Click (Mac/Linux)
```bash
chmod +x launch_gui.sh
./launch_gui.sh
```

## 📋 Prerequisites

### For GLMakie GUI
```julia
using Pkg
Pkg.add(["GLMakie", "DelimitedFiles", "Printf"])
Pkg.add(url="https://github.com/JuliaIO/MAT.jl")
```

### For Web GUI
```julia
using Pkg
Pkg.add(["Interact", "Blink", "Plots", "WebIO", "DelimitedFiles"])
Pkg.add(url="https://github.com/JuliaIO/MAT.jl")

# Then install Electron:
using Blink
Blink.AtomShell.install()
```

### Test Installation
```julia
include("test_gui_setup.jl")
```

## 🎯 Recommended Workflow

1. **First Time Setup**
   - Run `test_gui_setup.jl`
   - Install missing dependencies
   - Read `QUICKSTART.md`

2. **Using the GUI**
   - Launch via `launcher.jl`
   - Load a domain
   - Experiment with parameters
   - Generate and visualize grids

3. **If Problems Occur**
   - Check terminal for errors
   - Consult `TROUBLESHOOTING.md`
   - Run `test_gui_setup.jl` again

4. **Advanced Usage**
   - Read `README.md` for all features
   - See `ARCHITECTURE.md` for customization
   - Modify GUI files for custom needs

## 📚 Documentation Map

```
Start Here ────────> QUICKSTART.md
                           │
                           ▼
                     Try the GUI
                           │
        ┌──────────────────┼──────────────────┐
        │                  │                  │
        ▼                  ▼                  ▼
   Works Great      Need Details      Problems?
        │                  │                  │
        ▼                  ▼                  ▼
   Use & Enjoy       README.md      TROUBLESHOOTING.md
                           │
                           ▼
                    Want to Customize?
                           │
                           ▼
                   ARCHITECTURE.md
```

## 🎨 GUI Comparison

| Feature | GLMakie GUI | Web GUI |
|---------|-------------|---------|
| **Setup** | ⭐⭐⭐ | ⭐⭐⭐⭐ |
| **Performance** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ |
| **Graphics Quality** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ |
| **Large Grids** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ |
| **Interactivity** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ |
| **Cross-platform** | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| **Dependencies** | 3 packages | 5 packages |
| **Window Type** | Native | Browser |

**Recommendation**: 
- **GLMakie** for regular use, large grids, best performance
- **Web** for quick demos, teaching, or if OpenGL issues

## 🔧 Common Tasks

### Load an Airfoil Grid
```
1. Domain Type: Airfoil
2. Load Domain
3. Set splits: I=300,400  J=30
4. Generate Grid
```

### Simple Rectangle Test
```
1. Domain Type: Rectangle
2. Load Domain
3. Generate Grid (no splits needed)
```

### High-Quality Smoothing
```
1. Load domain
2. Max Iterations: 10000
3. Tolerance: 1e-7
4. All wall forcing: ON
5. Generate Grid
```

### Debugging Parameters
```
1. Disable splitting
2. Disable smoothing
3. Use edge solver only
4. Generate and verify
5. Add features back one at a time
```

## 💡 Tips & Tricks

### Performance
- Start with small grids for testing
- Use analytic solver (default)
- Disable verbose mode for speed
- GLMakie is faster for visualization

### Split Locations
- Indices from initial grid (before splitting)
- I-splits = vertical lines
- J-splits = horizontal lines
- Check grid size first!

### Smoothing
- ω = 0.2 is a safe default
- Lower ω = more stable, slower
- Higher ω = faster, may oscillate
- Wall forcing only where needed

### Visualization
- Use mouse to pan/zoom in GLMakie
- Watch terminal for error messages
- Status label shows current state
- Three panels show pipeline stages

## 🐛 Quick Troubleshooting

| Problem | Solution |
|---------|----------|
| GUI won't start | Run `test_gui_setup.jl` |
| Domain not found | Run from `examples/gui/` directory |
| Grid not generated | Click "Load Domain" first |
| Split error | Check indices < grid size |
| Slow performance | Reduce grid size or iterations |
| Blank plots | Check zoom level, resize window |
| Dependencies missing | See Prerequisites above |

## 📖 Learning Path

### Beginner
1. Read `QUICKSTART.md`
2. Run `test_gui_setup.jl`
3. Launch GUI with `launcher.jl`
4. Try rectangle example
5. Try airfoil example

### Intermediate
1. Read full `README.md`
2. Experiment with split locations
3. Try different metrics
4. Adjust smoothing parameters
5. Compare solver types

### Advanced
1. Read `ARCHITECTURE.md`
2. Modify GUI for custom domains
3. Add new visualization features
4. Export grids in custom formats
5. Integrate with workflows

## 🔗 Related Resources

- **Main Package**: `../../README.md`
- **Examples**: `../generalExample.jl`
- **Documentation**: `../../docs/`
- **Source Code**: `../../src/`
- **Tests**: `../../test/`

## 📝 Version Info

- GridGeneration.jl: v3.0.0
- Julia: 1.10.4
- Created: 2025
- Platform: Cross-platform (Windows, Mac, Linux)

## 🤝 Contributing

To improve the GUI:

1. Test changes with both GUIs
2. Update relevant documentation
3. Verify on multiple platforms
4. Add troubleshooting entries if needed

## 📄 License

Same as GridGeneration.jl package.

---

**Ready to start?**

```julia
cd("examples/gui")
julia launcher.jl
```

Enjoy interactive grid generation! 🎉
