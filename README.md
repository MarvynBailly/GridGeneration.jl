# GridGeneration

See the documentation and development for [GridGeneration.jl](https://marvyn.com/GridGeneration/dev/)

## To Do List
- [x] Fix grid spacing algorithm to normalize local center difference vectors. 
- [x] Bug with the ODE 1D to 2D projection. The orientation of the boundary is throwing off the distribution of points.
  - [x] Fix Projection Algorithm and add clarity to figures.
- [x] Write Semi analytical solver 
  - [x] Derived full analytic solution when $M(x)$ is explicitly provided. 
  - [x] Implement Solver
- [x] Tidy up "2D to 1D Reformulation"
- [x] Rewrite Numerical Solver for nonlinear second order BVP ODE. Use thomas algorithm to invert tridigaonal system with dirichlet boundary values.
  - [x] Create pages and files
  - [ ] Verification of Solver
    - [x] passes eye norm
    - [ ] (non-important) Check via MMS?
  - [x] Implement Solver
    - [x] FIX BUG in NumSolver causing points to cluster the wrong way. 
- [x] Compare analytic and semi-analytical results.
  - [x] For problems without analytic solution, passes "eye norm"
- [x] Pick a method to use for code  
  - [x]  Picked semi-analytical solution
- [ ] Test analytic solution with real metric data
  - [ ] Double check if metric nearest neighbors is worth the package overhead.
- [ ] Add section on Grid format 
- [ ] Add single block input with no splitting. User can upload a block, code with solve along the edges of the block, fill in with TFI, and generate the .grid files
- [ ] (hit here by 15th) Add support for single block input with user splitting. User can upload a block and split locations, code will creat multiblocks, solve along the edges of the block, fill in with TFI, and generate the .grid files
- [ ]  Add support for multiblock input with user splitting.
- [ ] Try custom splitting idea.


## Notes
### include package via github

Add the package to your project:

```julia
using Pkg
Pkg.add(url="https://github.com/MarvynBailly/GridGeneration.jl")
```

and include in your script:

```julia
using GridGeneration
```

### run build script website
julia --project=docs/ docs/make.jl

### to launch local server
julia --project=docs -e "using LiveServer; LiveServer.servedocs()"
