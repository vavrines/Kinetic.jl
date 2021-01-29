# Physical space

A struct `set <: AbstractPhysicalSpace` defines the geometric setup of a simulation.
For the structured topology, structs for 1 and 2 dimensional physical space are built.
```@docs
PSpace1D
PSpace2D
```
It contains
- x0 (y0): location of starting point.
- x1 (y1): location of ending point.
- nx (ny): number of cells in one direction.
- x (y): locations of middle points of all cells.
- dx (dy): intervals of all cell points.

Besides, a unstrctured mesh struct is built, which supports 1-3 dimensional geometries.
```@docs
UnstructMesh
```
It can be created by the built-in mesh reader.
```@docs
read_mesh
```