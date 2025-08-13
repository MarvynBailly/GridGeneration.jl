# Single Block with No Splitting
## 2D Single Block
Let's start by allowing the user to input a single [Tortuga](../GridFormat.md) block in the code. Thus the input will be the initial grid and the boundary information. The code will take in a 2D grid and solve the ODE along each edge of the domain. As the optimal number of points need not be the same for boundary edges across from each other, the code will proceed to resolve the edge that does not have the maximum number of points. Finally, the code will solve for the interior points using Transfinite Interpolation and update the boundary information with the new dimensions of the block. 


### Examples 
Let's use a small grid around an airfoil as an example. Since we are doing single block, let's not have a c-grid yet but keep one block around