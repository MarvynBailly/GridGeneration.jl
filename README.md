# GridGeneration

See the documentation and development for [GridGeneration.jl](https://marvyn.com/GridGeneration/dev/)

## To Do List
- [x] Fix grid spacing algorithm to normalize local center difference vectors. 
- [ ] Rewrite Numerical Solver for nonlinear second order BVP ODE. Use thomas alogrithm to invert tridigaonal system with dirichlet boundary values.
- [ ] Tidy up "2D to 1D Reformulation"
- [ ] Add section on Grid format
- [ ] Add single block input with no splitting. User can upload a block, code with solve along the edges of the block, fill in with TFI, and generate the .grid files
- [ ] Add support for single block input with user splitting. User can upload a block and split locations, code will creat multiblocks, solve along the edges of the block, fill in with TFI, and generate the .grid files
- [ ] Add support for multiblock input with user splitting.
- [ ] Try custom splitting idea.
