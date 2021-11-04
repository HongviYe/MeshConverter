# MeshConverter
A format converter for surface mesh intergrated with muli tools.
## Supported fileformat
Including ACSCII based `vtk`,`pls`,`facet`,`msh`,`obj`. 
## Converter file format
example
```shell
MeshConverter -i example.pls -t
```
Here `-i` followed by the input filename, with `-t` for output type.
type`MeshConveter --help` for more info.
## Supported tools
### Reorient the mesh
reset the oritation by DFS, the connected componment with biggest volume will be set pointed to the internal.
### Reverse the orientation
Reverse all the orientaion of facets
### Add bounding box
Add a bounding box by example mesh
### Remove degenerated mesh element
Repair vtk file for the area is equal to zero.

