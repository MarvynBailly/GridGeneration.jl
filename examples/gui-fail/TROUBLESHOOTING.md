# GUI Troubleshooting Guide

## Installation Issues

### "GLMakie not found"

**Problem**: `ERROR: ArgumentError: Package GLMakie not found`

**Solution**:
```julia
using Pkg
Pkg.add("GLMakie")
```

If that fails, try:
```julia
Pkg.update()
Pkg.add("GLMakie")
```

### "Blink not found" or Electron Issues

**Problem**: Web GUI won't start, Electron missing

**Solution**:
```julia
using Pkg
Pkg.add("Blink")

using Blink
Blink.AtomShell.install()
```

### MAT.jl Installation Fails

**Problem**: `Pkg.add("MAT")` doesn't work

**Solution**: Install from GitHub:
```julia
using Pkg
Pkg.add(url="https://github.com/JuliaIO/MAT.jl")
```

## Runtime Errors

### "initialGrid must be a 3D array"

**Problem**: Grid has wrong dimensions

**Cause**: Didn't load domain before generating

**Solution**: 
1. Click "Load Domain" button first
2. Verify initial grid appears in left panel
3. Then click "Generate Grid"

### "Domain Not Loaded" / Grid is 2×2×2

**Problem**: Trying to generate before loading

**Solution**:
```
1. Select domain type (Airfoil/Rectangle)
2. Click "Load Domain" ← MUST DO THIS
3. See grid appear
4. Configure parameters
5. Click "Generate Grid"
```

### "File not found: A-airfoil.txt"

**Problem**: Working directory is wrong

**Solution**: Run from correct directory:
```julia
# From GridGeneration root:
cd("examples/gui")
julia launcher.jl

# OR use absolute paths in domain loading
```

**Alternative**: Check file exists:
```julia
isfile("examples/airfoil/data/A-airfoil.txt")  # Should be true
```

### Split Index Out of Bounds

**Problem**: Error about split indices

**Cause**: Split location exceeds grid dimensions

**Example**:
```
Grid is 600×50
J-split: 100  ← ERROR! Max is 50
```

**Solution**:
- Check initial grid size (displayed after loading)
- I-splits must be ≤ Ni
- J-splits must be ≤ Nj
- Use smaller values

### "Method SolveODE not found"

**Problem**: GridGeneration.jl not loading properly

**Solution**:
```julia
# From examples/gui/:
include("../../src/GridGeneration.jl")
using .GridGeneration

# Verify it worked:
GridGeneration.GenerateGrid  # Should show function signature
```

## Display Issues

### Blank/Empty Plots

**Problem**: Grid generated but nothing shows

**GLMakie Solution**:
- Try resizing window
- Check zoom level (scroll to zoom out)
- Look for error messages in terminal

**Web GUI Solution**:
- Refresh browser window
- Check browser console (F12) for errors
- Try clicking widget again

### Grid Looks Distorted

**Problem**: Aspect ratio is wrong

**Solution**:
```julia
# GLMakie: Should have aspect=DataAspect()
# Check this is set in plot_grid!()

# Manual fix if needed:
ax.aspect = DataAspect()
```

### Can't See Split Lines

**Problem**: Entered splits but no red lines appear

**Causes & Solutions**:
1. **Wrong format**: Use commas, no spaces
   - ✗ `300 400`
   - ✓ `300,400`

2. **Not updating**: Press Enter after typing
   - Type in textbox
   - Press Enter/Return
   - Lines should appear

3. **Splits outside grid**: 
   - Check grid dimensions first
   - Use valid indices

### GLMakie Window Won't Open

**Problem**: Script runs but no window appears

**Solution 1 - Test GLMakie First**:
```julia
# Run the test script:
include("test_glmakie.jl")

# Should show a simple sine wave plot
# If this works, GLMakie is fine
# If not, there's an OpenGL/display issue
```

**Solution 2 - Check OpenGL**:
```julia
using GLMakie
GLMakie.activate!()

# Test simple plot:
fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, 1:10, 1:10)
display(fig)
```

**Solution 3 - Try different backend**:
```julia
# Use CairoMakie instead (slower but more compatible)
using CairoMakie
# (but you'll need to modify GUI code)
```

**Solution 3 - Use Web GUI instead**:
```julia
julia launcher.jl web
```

## Parameter Issues

### Smoothing Not Converging

**Problem**: "Max iterations reached" warning

**Symptoms**:
- Status shows warning
- Final error is large (>1e-3)
- Smoothed grid looks same as generated

**Solutions**:
1. **Reduce ω (relaxation)**:
   - Default: 0.2
   - Try: 0.1 or even 0.05
   
2. **Increase iterations**:
   - Default: 5000
   - Try: 10000 or 20000
   
3. **Check initial grid**:
   - Bad initial grid won't smooth well
   - Try without splitting first
   
4. **Adjust decay parameters**:
   - Too high (>0.7) can cause instability
   - Try 0.3-0.5 range

### Grid Generation Takes Forever

**Problem**: GUI freezes, nothing happens

**Causes**:
1. **Large grid size**: 
   - 500×500+ grids take 10-30 seconds
   - 1000×1000 can take minutes
   
2. **Too many iterations**:
   - Max_iter = 10000 takes longer
   - Reduce to 1000 for testing

**Solutions**:
- Check terminal for progress messages
- Wait patiently (it's computing, not frozen)
- Reduce grid size for testing
- Disable smoothing for faster tests
- Use analytic solver (faster than numeric)

### Wrong Boundary Solver Selected

**Problem**: Results look incorrect

**When to use each**:
- **analytic**: Smooth metrics, typical use (RECOMMENDED)
- **numeric**: Complex boundary conditions, debugging

**If unsure**: Use analytic

### Parameters Reset After Load

**Problem**: Changes don't stick

**Cause**: Loading domain doesn't preserve parameters

**Solution**: 
1. Load domain
2. THEN set parameters
3. Generate

(Parameters are independent of domain loading)

## Interaction Issues

### Sliders Not Responding

**Problem**: Moving slider doesn't update value

**GLMakie**: Should update label immediately
- If not, try clicking elsewhere first
- Check console for errors

**Web GUI**: Might need to release slider
- Move slider
- Release mouse button
- Value should update

### Buttons Don't Work

**Problem**: Clicking button does nothing

**Debug Steps**:
```julia
# Check if handler is defined:
# Should see "on(button.clicks) do n" in code

# Check for errors in terminal
# Click button and watch console
```

**Common Cause**: Event handler crashed
- Look for stack trace in terminal
- Fix error and restart GUI

### Text Input Not Parsing

**Problem**: Entered "300,400" but not recognized

**Solutions**:
1. **Press Enter** after typing
2. **Check format**:
   - ✓ `300,400`
   - ✗ `300, 400` (space)
   - ✗ `300.400` (period)
   - ✗ `300;400` (semicolon)

3. **Verify in terminal**:
```julia
# Check what was parsed:
println(state.i_splits)  # Should show [300, 400]
```

## Performance Issues

### Slow Visualization

**Problem**: Plotting takes long time

**GLMakie Solutions**:
- Reduce line width (default 0.5 is good)
- Don't show_points unless needed
- Use fewer grid lines for preview

**Web GUI Solutions**:
- Plots.jl is slower for large grids
- Consider switching to GLMakie GUI
- Reduce plot complexity

### Memory Issues

**Problem**: "Out of memory" or computer slows

**Cause**: Grid too large

**Solution**:
- Reduce grid size
- Close other applications
- Restart Julia session
- Consider running overnight for very large grids

## Platform-Specific Issues

### Windows: OpenGL Error

**Problem**: "OpenGL 3.3+ required"

**Solution**:
- Update graphics drivers
- Use Web GUI instead
- Run on different machine

### macOS: Security Warning

**Problem**: "App from unidentified developer"

**Solution**:
- Right-click launcher, select "Open"
- Or use System Preferences > Security

### Linux: Display Issues

**Problem**: X11 errors, display not working

**Solution**:
```bash
# Set display:
export DISPLAY=:0

# Or use web GUI (doesn't need X11)
```

## Getting More Help

### Enable Verbose Mode

```julia
# In GUI, toggle "Verbose Output"
# Or in code:
verbose=true

# This shows:
# - Iteration progress
# - Convergence info
# - Detailed errors
```

### Check Terminal Output

Most errors appear in terminal, not GUI!
- Always watch terminal while using GUI
- Look for red error messages
- Stack traces show what went wrong

### Minimal Test Case

If still having issues, try minimal example:

```julia
using GridGeneration

# Simple 5×5 grid
N = 5
top = [range(0, 1, length=N) ones(N)]
right = [ones(N) range(1, 0, length=N)]
bottom = [range(1, 0, length=N) zeros(N)]
left = [zeros(N) range(0, 1, length=N)]

grid = TFI([top, right, bottom, left])
M(x, y) = [1.0, 1.0]

params = SimParams(
    useSplitting=false,
    useEdgeSolver=false,
    useSmoothing=false
)

result = GenerateGrid(grid, [], [], M; params=params)
println("Success!")
```

If this works, issue is in GUI setup, not core library.

### Report Issue

If nothing works:

1. **Run test script**:
```julia
include("test_gui_setup.jl")
```

2. **Collect info**:
   - Julia version: `versioninfo()`
   - Error message (full text)
   - Steps to reproduce

3. **File GitHub issue** with:
   - Test script output
   - Error messages
   - Platform (Windows/Mac/Linux)
   - What you tried

## Quick Diagnostics

Run this to check everything:

```julia
# From examples/gui/:
include("test_gui_setup.jl")

# Should show:
# ✅ for working components
# ⚠️  for missing dependencies
# ❌ for errors
```

Fix any ⚠️  or ❌ items first, then try GUI again.
