# Quick Start Guide - GridGeneration GUI

## 30-Second Start

```julia
# Navigate to GUI directory
cd("examples/gui")

# Run launcher (will prompt for GUI choice)
julia launcher.jl
```

## 2-Minute Tutorial

### Step 1: Launch GUI

```julia
julia launcher.jl makie    # or 'web' for browser version
```

### Step 2: Load a Domain

1. Select "Airfoil" from Domain Type dropdown
2. Click "Load Domain" button
3. See the initial grid displayed

### Step 3: Configure Splits (Optional)

1. Check "Use Block Splitting"
2. Enter split locations:
   - I-splits: `300,400`
   - J-splits: `30`
3. See split lines highlighted on grid

### Step 4: Generate Grid

1. Keep default solver settings (analytic)
2. Click "Generate Grid"
3. View three panels:
   - Initial grid
   - Generated grid (after splitting + edge solving)
   - Smoothed grid (after elliptic smoothing)

## What the Parameters Mean

### Quick Reference

| Parameter | What It Does | Typical Value |
|-----------|-------------|---------------|
| **I-splits** | Vertical split lines | `300,400` |
| **J-splits** | Horizontal split lines | `30` |
| **Boundary Solver** | How edges are solved | `analytic` |
| **Max Iterations** | SOR iteration limit | `5000` |
| **Tolerance** | Convergence criterion | `1e-6` |
| **ω (omega)** | Relaxation factor | `0.2` |
| **Decay Parameter** | Wall forcing strength | `0.4` |

### Understanding Splits

```
Initial Grid: 600×50
I-splits: 300,400
J-splits: 30

Creates blocks:
  ┌─────┬─────┬─────┐
  │  3  │  4  │  5  │  ← top (j > 30)
  ├─────┼─────┼─────┤
  │  0  │  1  │  2  │  ← bottom (j ≤ 30)
  └─────┴─────┴─────┘
    ↑     ↑     ↑
  i≤300 300<i≤400 i>400
```

## Common Workflows

### Airfoil C-Grid with Clustering

```julia
# In GUI:
Domain: Airfoil
Metric Problem: 6
Metric Scale: 0.05
I-splits: 300,400
J-splits: 30
Boundary Solver: analytic
Use Smoothing: ✓
```

### Simple Rectangle Test

```julia
# In GUI:
Domain: Rectangle
Metric Problem: 1
I-splits: (leave empty)
J-splits: (leave empty)
Use Smoothing: ✗
```

### High-Quality Smooth Grid

```julia
# In GUI:
Max Iterations: 10000
Tolerance: 1e-7
ω: 0.2
All Wall Forcing: ✓
Decay Parameter: 0.5
```

## Keyboard Shortcuts (GLMakie GUI)

- **Mouse Drag**: Pan view
- **Scroll**: Zoom in/out
- **Right Click + Drag**: Rotate (3D views)
- **Ctrl+C**: Copy view to clipboard (when focus on plot)

## Tips for Best Results

### 🎯 Splitting Strategy
- Start with no splits, verify base grid works
- Add splits one at a time to isolate issues
- Split indices must be within initial grid bounds

### ⚡ Performance
- Large grids (>500×500) take 10-30 seconds
- Disable verbose mode for speed
- Use analytic solver when possible

### 🎨 Visualization
- Zoom in to inspect grid quality
- Look for:
  - Smooth transitions between blocks
  - Orthogonality at walls
  - Point clustering near geometry

### 🔧 Troubleshooting

**Grid looks wrong?**
- Check boundary order: [top, right, bottom, left]
- Verify split indices are correct
- Try simpler configuration first

**Smoothing not converging?**
- Reduce ω (try 0.1)
- Increase max iterations
- Check initial grid quality

**GUI frozen?**
- Large grid generation is computing
- Check terminal for progress messages
- Wait or reduce grid size

## Example Session

```julia
# 1. Start GUI
cd("examples/gui")
julia launcher.jl makie

# 2. In GUI window:
#    - Select Airfoil
#    - Click "Load Domain"
#    - Enter splits: I=300,400  J=30
#    - Click "Generate Grid"

# 3. Observe:
#    - Initial grid with red split lines
#    - Blue generated blocks
#    - Green smoothed final grid

# 4. Experiment:
#    - Change decay parameter
#    - Toggle wall forcing
#    - Regenerate and compare
```

## Next Steps

- **Read full documentation**: `examples/gui/README.md`
- **Try custom metrics**: Modify metric problem/scale
- **Explore source**: See how grid generation works
- **Script version**: Use `examples/generalExample.jl` for batch processing

## Getting Help

- Check `examples/gui/README.md` for detailed docs
- See `examples/generalExample.jl` for scripting examples
- Review `.github/copilot-instructions.md` for architecture
- File issues on GitHub for bugs
