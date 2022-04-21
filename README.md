# SIMP-topology-optimization-with-Isoperimetric-Quadrilateral-Finite-Element
 Isoperimetric Quadrilateral  Finite Element mesh with SIMP topology optimization in 2D
 
The code uses  Isoperimetric Quadrilateral Finite Element analysis
All of shape functions are embeded, so the code can be used for any general problem
The mesh can be generated or imported from other meshing software like ANSYS
The mesh information is the node and element information, el4.txt and nd4.txt are two examples from ANSYS
The elemnts information is four locations of each node. If the mesh is triangular, it will be ignored in the TO
The node information is the name of node and location in x-y plane
loading and boundry condition can be defined in the main code TO3
THe passive elements are the ones to be rigid or void all the time

Run the TO3 and change the parameters if necessary
