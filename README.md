# SIMP-topology-optimization-with-Isoperimetric-Quadrilateral-Finite-Element
 Isoperimetric Quadrilateral  Finite Element mesh with SIMP topology optimization in 2D
 
The code uses  Isoperimetric Quadrilateral Finite Element analysis for minimizing compliance for restricted volme fraction
All of shape functions are embeded, so the code can be used for any general problem in 2D
The optimizer is optimally criteria and also uses sensitivity filter ( this part is from famous 99-line code sigmund) (https://link.springer.com/article/10.1007/s001580050176) 

The mesh can be generated or imported from other meshing software like ANSYS
The mesh information is the node and element information, el4.txt and nd4.txt are two examples from ANSYS
The elemnts information is four locations of each node. If the mesh is triangular, it will be ignored in the TO
The node information is the name of node and location in x-y plane
loading and boundry condition can be defined in the main code TO3
THe passive elements are the ones to be rigid or void all the time

Run the TO3 and change the parameters if necessary
