### Organization

## graphs.py file

#### BaseGraph class
Base class for the graph classes. 

    Objects
    - config
    - edges
    - numEdges

    - numNodes
    - nodes

    - weight_matrix
    - degree_matrix
    - laplacian
    - pinv_laplacian
    
    Get functions:
    - get_edges
    - get_nodes
    - get_laplacian
    - get_pinv_laplacian
    - get_edge_weights
    
    Private methods:
    - construct_weight_matrix
    - construct_degree_matrix
    - construct_laplacian
    - construct_pinv_laplacian
    

#### ConnectionGraph(BaseGraph) class
Connection graphs. Inherits from base graph

    Objects
    - connectionDim
    - edgeRotations

    - connection_incidence_matrix
    - connection_incidence_matrix_chung2014
    - connection_laplacian
    - connection_laplacian_chung2014
    - pinv_connection_laplacian
    - pinv_connection_laplacian_chung2014
    
    Private methods:
    - construct_connection_incidence_matrix
    - construct_connection_laplacian
    - construct_connection_incidence_matrix_chung2014: Note the different implementation 
    of this matrix from the connection incidence matrix used in our definition 6.7
    - construct_connection_laplacian_chung2014
    - construct_pinv_connection_laplacian_chung2014
    
    Get methods:
    - get_connectionDim
    - get_numEdges
    - get_numNodes
    - get_connection_laplacian
    - get_connection_laplacian_chung2014
    - get_connection_incidence_matrix
    - get_pinv_connection_laplacian
    - get_pinv_connection_laplacian_chung2014
    - get_connection_incidence_matrix_chung2014
    - get_connection_laplacian_chung2014
    - get_pinv_connection_laplacian_chung2014
    
    
## connection_resistance.py file

#### ConnectionResistance class
The ConnectionResistance class contains functionality for calculating the connection resistance 
in a connection resistance graph

    Objects:
    - connection_graph: Connection graph
    - config: Configuration file of connection graph
    
    Get functions:
    - get_edges: Get edges from connection graph
    
    Public methods:
    - calc_effective_resistance
    - calc_connection_resistance_ij
    - calc_connection_resistance
    - calc_connection_resistance_chung2014
    
    
## rotation_matrices.py file
Contains functions that calculate the rotation matrices in 1D, 2D and 3D around the x,y and z axis

## util.py file
#### schur_comp function
Calculates the schur complement of the input matrix