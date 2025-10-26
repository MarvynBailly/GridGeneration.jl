# GridGeneration.jl Interactive GUI

This directory contains interactive graphical user interfaces for GridGeneration.jl, allowing you to visualize grids and configure parameters through an intuitive interface.

## Available GUIs

### 1. **GLMakie GUI** (Recommended for Desktop)
- **File**: `GridGenerationGUI.jl`
- **Backend**: GLMakie.jl (native OpenGL-based rendering)
- **Pros**: 
  - Fast, responsive, native performance
  - Beautiful high-quality graphics
  - Better for large grids
  - Interactive zoom/pan built-in
- **Cons**: 
  - Requires OpenGL support
  - Slightly larger dependency footprint

### 2. **Web GUI** (Recommended for Simplicity)
- **File**: `GridGenerationWebGUI.jl`
- **Backend**: Interact.jl + Blink.jl (web-based)
- **Pros**:
  - Browser-based interface
  - Easier to set up
  - Familiar web UI controls
- **Cons**:
  - Slightly slower for very large grids
  - Requires Chrome/Electron

## Installation

### Prerequisites

First, ensure GridGeneration.jl is set up:

```julia
using Pkg
Pkg.activate(".")
```

### For GLMakie GUI

```julia
Pkg.add(["GLMakie", "DelimitedFiles"])
Pkg.add(url="https://github.com/JuliaIO/MAT.jl")
```

### For Web GUI

```julia
Pkg.add(["Interact", "Blink", "Plots", "WebIO", "DelimitedFiles"])
Pkg.add(url="https://github.com/JuliaIO/MAT.jl")
```

## Usage

### Launching the GLMakie GUI

From the project root directory:

```julia
include("examples/gui/GridGenerationGUI.jl")
```

Or directly:

```julia
cd("examples/gui")
julia GridGenerationGUI.jl
```

### Launching the Web GUI

From the project root directory:

```julia
include("examples/gui/GridGenerationWebGUI.jl")
```

Or directly:

```julia
cd("examples/gui")
julia GridGenerationWebGUI.jl
```

## GUI Features

### Domain Setup
1. **Select Domain Type**: Choose between Airfoil and Rectangle geometries
2. **Metric Configuration**: 
   - Problem number (1-10 for different clustering patterns)
   - Scale parameter for metric strength
3. **Load Domain**: Click to load the selected domain and visualize initial grid

### Block Splitting
- **Toggle Splitting**: Enable/disable multi-block splitting
- **I-splits**: Enter comma-separated indices for vertical split lines (e.g., `300,400`)
- **J-splits**: Enter comma-separated indices for horizontal split lines (e.g., `30`)
- Split lines are highlighted in red on the initial grid view

### Solver Selection
- **Edge Solver**: Enable metric-adaptive point distribution on block edges
- **Boundary Solver Type**:
  - `analytic`: Semi-analytical integration (faster, smooth metrics)
  - `numeric`: Thomas algorithm (more robust)

### Smoothing Options
- **Enable Smoothing**: Toggle elliptic smoothing
- **Method**: Currently supports `ellipticSS` (SOR iteration)
- **Elliptic Parameters**:
  - **Max Iterations**: Maximum SOR iterations (100-10000)
  - **Tolerance**: Convergence tolerance (1e-8 to 1e-4)
  - **Relaxation ω**: SOR relaxation parameter (0.1-2.0, typically 0.2)
  - **Decay Parameter**: Wall forcing decay rate (0.1-1.0)
  - **Wall Forcing**: Toggle forcing for each wall (top/bottom/left/right)
  - **Verbose**: Print convergence information

### Grid Generation
1. Configure all parameters
2. Click **Generate Grid**
3. View results in three panels:
   - **Initial Grid**: TFI-generated grid with split lines
   - **Generated Grid**: After splitting and edge solving
   - **Smoothed Grid**: After elliptic smoothing (if enabled)

## Workflow Example

### Typical Airfoil Grid Generation

1. **Load Domain**
   - Domain Type: Airfoil
   - Metric Problem: 6
   - Metric Scale: 0.05
   - Click "Load Domain"

2. **Configure Splitting**
   - Use Block Splitting: ✓
   - I-splits: `300,400`
   - J-splits: `30`

3. **Configure Solvers**
   - Use Edge Solver: ✓
   - Boundary Solver: analytic

4. **Configure Smoothing**
   - Use Smoothing: ✓
   - Max Iterations: 5000
   - Tolerance: 1e-6
   - ω: 0.2
   - Decay: 0.4
   - All walls: ✓

5. **Generate**
   - Click "Generate Grid"
   - Wait for processing
   - View results in all three panels

### Rectangle Grid Generation

1. **Load Domain**
   - Domain Type: Rectangle
   - Metric Problem: 1
   - Click "Load Domain"

2. **Simple Configuration**
   - Leave splitting disabled for uniform grid
   - Or enable with custom splits

3. **Generate**
   - Click "Generate Grid"

## Tips & Best Practices

### Split Locations
- Indices refer to the **initial TFI grid**
- I-splits create vertical split lines (vary in ξ-direction)
- J-splits create horizontal split lines (vary in η-direction)
- Splits automatically propagate across block interfaces

### Metric Configuration
- Lower scale values → tighter clustering
- Higher problem numbers → different clustering patterns
- For airfoils: problems 5-7 work well with scale 0.05-0.1

### Elliptic Smoothing
- Start with ω = 0.2 for stability
- Increase for faster convergence (may become unstable)
- Decay 0.4-0.5 gives moderate wall orthogonality
- Enable wall forcing only where needed (typically bottom wall for airfoils)

### Performance
- Large grids (>500×500) may take several seconds
- Disable verbose mode for faster generation
- Use analytic solver when possible (faster)
- GLMakie GUI is faster for visualization of large grids

## Troubleshooting

### GUI Won't Launch

**GLMakie**: Check OpenGL support
```julia
using GLMakie
GLMakie.activate!()
```

**Web GUI**: Ensure Blink.jl can find Electron
```julia
using Blink
Blink.AtomShell.install()
```

### "Domain Not Loaded" Error
- Make sure you clicked "Load Domain" before "Generate Grid"
- Check that airfoil data file exists: `examples/airfoil/data/A-airfoil.txt`

### Grid Generation Fails
- Verify split indices are within grid bounds
- Check metric scale is reasonable (not too small/large)
- Try simpler configuration (disable splitting first)

### Plots Don't Update
- For Web GUI: Click the widget again or refresh
- For GLMakie: Try resizing the window

## Advanced Customization

### Adding New Domain Types

Edit the domain loading section in either GUI file:

```julia
if state.case == :custom
    # Your domain setup here
    state.initialGrid = custom_grid_function()
    state.bndInfo = custom_boundary_info()
    state.interInfo = []
    state.metric = custom_metric_function()
end
```

### Custom Metrics

Modify metric selection:

```julia
# Add to metric selection logic
state.metric = make_getMetric(
    custom_geometry,
    scale=scale_slider.value[],
    profile=:gauss  # or :rational
)
```

### Export Functionality

Add save buttons to export grids:

```julia
# GLMakie version
save_btn = Button(controls[27, 1:2], label="Save Grid")
on(save_btn.clicks) do n
    # Use GridGeneration's Tortuga export or custom format
    writeTurtleGrid(state.smoothBlocks, ...)
end
```

## API Reference

### Main GUI Functions

#### GLMakie Version
- `launch_gui()`: Launch the GLMakie-based GUI
  - Returns: `(fig, state)` - Makie Figure and GUIState
  
- `plot_grid!(ax, block; kwargs...)`: Plot single grid block on axis
- `plot_blocks!(ax, blocks; kwargs...)`: Plot multiple blocks
- `highlight_splits!(ax, block, i_splits, j_splits; kwargs...)`: Show split lines

#### Web Version
- `launch_web_gui()`: Launch web-based GUI
  - Returns: Blink window object
  
- `create_web_gui()`: Create UI without launching window
  - Returns: Interact UI object
  
- `plot_grid(block; kwargs...)`: Plot single grid using Plots.jl
- `plot_blocks(blocks; kwargs...)`: Plot multiple blocks

### GUIState Structure

```julia
mutable struct GUIState
    case::Symbol                              # :airfoil or :rectangle
    initialGrid::Array{Float64, 3}            # Initial TFI grid [2, Ni, Nj]
    bndInfo::Vector{Any}                      # Boundary conditions
    interInfo::Vector{Any}                    # Interface information
    metric::Function                          # Metric field M(x,y)
    i_splits::Vector{Int}                     # Vertical split indices
    j_splits::Vector{Int}                     # Horizontal split indices
    # ... solver parameters ...
    resultBlocks::Union{Nothing, Vector}      # Generated blocks
    smoothBlocks::Union{Nothing, Vector}      # Smoothed blocks
end
```

## Contributing

To add features or improve the GUI:

1. Follow the existing code structure
2. Test with both airfoil and rectangle cases
3. Ensure error handling is in place
4. Update this README with new features

## See Also

- Main package documentation: [GridGeneration.jl](../../docs/)
- Example scripts: [examples/](../)
- Core implementation: [src/](../../src/)
