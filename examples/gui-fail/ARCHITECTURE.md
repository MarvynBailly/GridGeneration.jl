# GridGeneration GUI Architecture

## Application Flow

```
┌─────────────────────────────────────────────────────────────┐
│                     User Launches GUI                        │
│                   (julia launcher.jl)                        │
└────────────────────┬────────────────────────────────────────┘
                     │
        ┌────────────┴────────────┐
        │   Dependency Check      │
        │  - GLMakie or Blink     │
        │  - GridGeneration.jl    │
        │  - Example helpers      │
        └────────────┬────────────┘
                     │
        ┌────────────▼────────────┐
        │    GUI Initializes      │
        │  - Create controls      │
        │  - Setup visualization  │
        │  - Bind event handlers  │
        └────────────┬────────────┘
                     │
                     │
┌────────────────────▼────────────────────────┐
│              USER INTERACTIONS               │
└──────────────────────────────────────────────┘

Step 1: LOAD DOMAIN
┌─────────────────────────────────────┐
│  User selects:                      │
│  • Domain type (Airfoil/Rectangle)  │
│  • Metric problem (1-10)            │
│  • Metric scale                     │
│  Clicks: "Load Domain"              │
└────────────┬────────────────────────┘
             │
    ┌────────▼─────────┐
    │ Load boundaries  │
    │ Create TFI grid  │
    │ Setup metric M   │
    └────────┬─────────┘
             │
    ┌────────▼─────────┐
    │ Display initial  │
    │ grid in viewer   │
    └──────────────────┘

Step 2: CONFIGURE SPLITS (Optional)
┌─────────────────────────────────────┐
│  User enters:                       │
│  • I-splits: "300,400"              │
│  • J-splits: "30"                   │
└────────────┬────────────────────────┘
             │
    ┌────────▼─────────┐
    │ Parse indices    │
    │ Highlight splits │
    │ on initial grid  │
    └──────────────────┘

Step 3: SET PARAMETERS
┌─────────────────────────────────────┐
│  Solver Settings:                   │
│  • Boundary solver (analytic)       │
│  • Use edge solver (✓)              │
│                                     │
│  Smoothing Settings:                │
│  • Enable smoothing (✓)             │
│  • Max iterations (5000)            │
│  • Tolerance (1e-6)                 │
│  • Omega (0.2)                      │
│  • Wall forcing toggles             │
└─────────────────────────────────────┘

Step 4: GENERATE GRID
┌─────────────────────────────────────┐
│  User clicks: "Generate Grid"       │
└────────────┬────────────────────────┘
             │
    ┌────────▼─────────────────────────────┐
    │     GridGeneration Pipeline          │
    ├──────────────────────────────────────┤
    │  1. Split blocks (if enabled)        │
    │     └─> Multi-block with propagation │
    │                                      │
    │  2. Solve edges (if enabled)         │
    │     └─> Metric-adaptive distribution │
    │                                      │
    │  3. Elliptic smooth (if enabled)     │
    │     └─> SOR iteration with forcing   │
    └────────┬─────────────────────────────┘
             │
    ┌────────▼─────────┐
    │ Display results: │
    │ • Initial grid   │
    │ • Generated grid │
    │ • Smoothed grid  │
    └──────────────────┘
```

## GUI Component Architecture

### GLMakie GUI Structure

```
┌───────────────────────────────────────────────────────┐
│                    Makie Figure                       │
├───────────────────┬───────────────────────────────────┤
│                   │                                   │
│   CONTROL PANEL   │      VISUALIZATION PANELS         │
│   (Left Column)   │      (Right Columns)              │
│                   │                                   │
│  ┌─────────────┐  │  ┌─────────┐  ┌──────────┐       │
│  │Domain Select│  │  │ Initial │  │Generated │       │
│  │  Dropdown   │  │  │  Grid   │  │  Grid    │       │
│  └─────────────┘  │  │ Viewer  │  │  Viewer  │       │
│  ┌─────────────┐  │  │         │  │          │       │
│  │Load Button  │  │  └─────────┘  └──────────┘       │
│  └─────────────┘  │                                   │
│  ┌─────────────┐  │  ┌──────────────────────┐        │
│  │Split Inputs │  │  │    Smoothed Grid     │        │
│  │  Textboxes  │  │  │       Viewer         │        │
│  └─────────────┘  │  │                      │        │
│  ┌─────────────┐  │  └──────────────────────┘        │
│  │Solver Menus │  │                                   │
│  └─────────────┘  │                                   │
│  ┌─────────────┐  │    Interactive Features:         │
│  │Param Sliders│  │    • Pan with drag               │
│  └─────────────┘  │    • Zoom with scroll            │
│  ┌─────────────┐  │    • Auto-scaling                │
│  │Wall Toggles │  │                                   │
│  └─────────────┘  │                                   │
│  ┌─────────────┐  │                                   │
│  │Generate Btn │  │                                   │
│  └─────────────┘  │                                   │
│  ┌─────────────┐  │                                   │
│  │Status Label │  │                                   │
│  └─────────────┘  │                                   │
└───────────────────┴───────────────────────────────────┘
```

### Web GUI Structure

```
┌────────────────────────────────────────────────────┐
│              Blink Window (Browser)                │
├────────────────────┬───────────────────────────────┤
│                    │                               │
│   CONTROL PANEL    │    VISUALIZATION PANELS       │
│   (HTML Layout)    │    (Plots.jl)                 │
│                    │                               │
│  • Domain dropdown │    ┌──────────────────┐       │
│  • Load button     │    │  Initial Grid    │       │
│  • Metric sliders  │    │  (Plots.jl)      │       │
│  • Split textboxes │    └──────────────────┘       │
│  • Solver dropdown │    ┌──────────────────┐       │
│  • Param sliders   │    │ Generated Grid   │       │
│  • Toggle switches │    │  (Plots.jl)      │       │
│  • Generate button │    └──────────────────┘       │
│  • Status div      │    ┌──────────────────┐       │
│                    │    │  Smoothed Grid   │       │
│                    │    │  (Plots.jl)      │       │
│                    │    └──────────────────┘       │
└────────────────────┴───────────────────────────────┘
      │                           │
      └────── WebIO/Interact ─────┘
              (Reactive Binding)
```

## Data Flow

```
User Input
    │
    ▼
Observable/Signal
    │
    ▼
Event Handler
    │
    ├─> Update GUIState
    │
    ├─> Call GridGeneration
    │   │
    │   ├─> TFI
    │   ├─> SplitBlock
    │   ├─> SolveAllBlocks
    │   └─> SmoothBlocks
    │
    └─> Update Visualization
        │
        ├─> Clear axes
        ├─> Plot new data
        └─> Refresh display
```

## State Management

```julia
GUIState (Mutable Struct)
├── Domain Data
│   ├── case::Symbol              # :airfoil or :rectangle
│   ├── initialGrid::Array        # [2, Ni, Nj]
│   ├── bndInfo::Vector           # Boundary conditions
│   ├── interInfo::Vector         # Interface info
│   └── metric::Function          # M(x,y) -> [M11, M22]
│
├── User Parameters
│   ├── i_splits::Vector{Int}     # Vertical split indices
│   ├── j_splits::Vector{Int}     # Horizontal split indices
│   ├── useSplitting::Bool
│   ├── boundarySolver::Symbol    # :analytic or :numeric
│   ├── useSmoothing::Bool
│   └── elliptic params...
│
└── Results
    ├── resultBlocks::Vector      # After edge solving
    └── smoothBlocks::Vector      # After smoothing
```

## Event Handling Pattern

```
┌──────────────┐
│ User Action  │ (click button, change slider, etc.)
└──────┬───────┘
       │
       ▼
┌──────────────┐
│  Observable  │ (Makie) or Signal (Interact)
│   Triggers   │
└──────┬───────┘
       │
       ▼
┌──────────────┐
│   Handler    │ on(observable) do value
│   Function   │     # Update state
└──────┬───────┘     # Call GridGen
       │             # Update display
       │         end
       ▼
┌──────────────┐
│Update Display│ (clear axes, plot new data)
└──────────────┘
```

## File Dependencies

```
GridGenerationGUI.jl
├── GLMakie (UI framework)
├── Printf (formatting)
├── DelimitedFiles (file I/O)
├── MAT (MATLAB files)
├── GridGeneration.jl
│   └── [All core modules]
├── airfoil.jl
│   ├── GetAirfoilSetup()
│   └── GetAirfoilMetric()
└── rectangle.jl
    ├── GetRectangleDomain()
    └── GetRectangleMetric()

GridGenerationWebGUI.jl
├── Interact (widgets)
├── Blink (window)
├── Plots (visualization)
├── WebIO (reactive)
└── [Same dependencies as above]
```

## Extension Points

### Adding New Domain Types

```julia
# In domain_select handler:
elseif state.case == :custom
    state.initialGrid = LoadCustomDomain()
    state.metric = CustomMetric()
    # etc.
```

### Adding New Visualizations

```julia
# Add new axis:
ax_custom = Axis(viz_layout[5, 1:2])

# Update in generate handler:
plot_custom_analysis!(ax_custom, state.resultBlocks)
```

### Custom Export

```julia
# Add export button:
export_btn = Button(controls[28, 1:2], label="Export Grid")

on(export_btn.clicks) do n
    write_custom_format(state.smoothBlocks, "output.dat")
end
```

## Performance Considerations

```
Grid Size     │ TFI    │ Splitting │ Solving │ Smoothing │ Total
──────────────┼────────┼───────────┼─────────┼───────────┼──────
50×50         │ <0.1s  │ <0.1s     │ 0.5s    │ 1s        │ ~2s
200×200       │ <0.1s  │ 0.5s      │ 2s      │ 5s        │ ~8s
500×500       │ 0.2s   │ 2s        │ 10s     │ 20s       │ ~32s
1000×1000     │ 1s     │ 8s        │ 40s     │ 80s       │ ~130s
```

Recommendations:
- Use analytic solver (faster than numeric)
- Reduce max_iter for testing
- Start with smaller grids
- GLMakie faster for large grid visualization
