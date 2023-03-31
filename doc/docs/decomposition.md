<style>
p.function{
    color:#357ab7;
    background:#e0e8f3;
    border-top:4px solid #5897cf;
    padding:0.3em 0.3em 0.3em 0.3em;
    margin-bottom:0em;
    font-weight: bold;
}
</style>

# Decomposition

This chapter describes functions for performing Haar wavelet transform.  
Some of these type structure are based on **gsl** or **Armadillo** library.  
The wavelet functions are declared in the header file **`decomposition.h`**.
The mapping functions are declared in the header files **`mapping_BFS.h`** and **`mapping_SVD.h`**.
<hr/>

## Organization

This module is internally divided into 3 submodules :  
 - **`decomposition.h`** for wavelet decomposition with different mapping  
 - **`mapping_BFS.h`** to run **B**readth-**Fi**rst **S**earch mapping  
 - **`mapping_SVD.h`** to run **S**ingular **V**alue **D**ecomposition mapping
<hr/>

## Initialization

### Transform workspace

**Type sp_workspace**  
This structure contains a sparse matrix to hold intermediate results the transform.

<p class="function">sp_workspace * sp_workspace_graph_alloc (size_t size1, size_t size2 = 1)</p>  
This function allocates memory for structure variable of sp_workspace.


<p class="function">void sp_workspace_delete (sp_workspace * work)</p>   
This function releases allocated memory of sp_workspace structure variable.

### Decomposition Result Matrix

The wavelet decomposition will segment data into 2 parts : **A**ccuracy and **D**ifference.
To hold these results, we use a dense matrix A for the **A**ccuracy part and a sparse matrix for the **D**ifference part.  
The first dimension of these matrices depend on the size of the data and number of level you run. The second dimension is the number of graph in your graph sequence (1 for a simple graph).

<p class="function">size_t size1A (size_t n, uint level)</p>   
This function returns the size of the first dimension for a dense matrix A to hold result from `level` wavelet transform of size n. (use for allocation)

<p class="function">size_t size1D (size_t n, uint level)</p>   
This function returns the size of the first dimension for a sparse matrix D to hold result from `level` wavelet transform of size n. (use for allocation)

<p class="function">gsl_matrix * matrixA_alloc (size_t dim1, size_t dim2, uint level)</p>  
<p class="function">gsl_matrix * matrixA_alloc (map_edge *m, uint level, size_t n_graphs = 1)</p>  
<p class="function">gsl_matrix * matrixA_alloc (graph *g, uint level)</p>  
These functions allocate memory for dense matrix A.


<p class="function">gsl_spmatrix * matrixD_alloc (size_t dim1, size_t dim2, uint level)</p>  
<p class="function">gsl_spmatrix * matrixD_alloc (map_edge *m, uint level, size_t n_graphs = 1)</p>  
<p class="function">gsl_spmatrix * matrixD_alloc (graph *g, uint level)</p>  
These functions allocate memory for sparse matrix D.

<p class="function">void matrixA_delete (gsl_matrix *A)</p>   
This function releases allocated memory of dense matrix A.

<p class="function">void matrixD_delete (gsl_spmatrix *D)</p>   
This function releases allocated memory of sparse matrix D.
<hr/>

## Mapping

There is map to store unique if of each edge (**[map_edge](data_processing.md#linkstream-types)**) and other for each node (**[map_node](data_processing.md#linkstream-types)**).

### BFS

This mapping follows **B**readth-**Fi**rst **S**earch order of a graph to give id to each edges and hold the result in a **[map_edge](data_processing.md#linkstream-types)**.

<p class="function">map_edge * map_edge_alloc ( )</p>  
This function allocates memory for map_edge structure variable.

<p class="function">void map_edge_delete (map_edge *m)</p>  
This function releases allocated memory of map_edge structure variable.

<p class="function">int mapping_BFS (graph *g, map_edge *m)</p>  
This function fills **`m`** (**[map_edge](data_processing.md#linkstream-types)**) by processing BFS algorithm on **`g`** (**[graph](data_processing.md#linkstream-types)**).

### SVD

This mapping computes the adjacency matrix of a graph and does **S**ingular **V**alue **D**ecomposition recursively on it to relabelized each node. Then we can use the function fct_edge_graph who returns an unique id ( |VxV| ) for each edge by giving the label of the twice nodes.

<p class="function">map_node * map_node_alloc ( )</p>  
This function allocates memory for map_node structure variable.

<p class="function">void map_node_delete (map_node *m)</p>  
This function releases allocated memory of map_node structure variable.

<p class="function">int mapping_SVD (graph *g, map_node *m)</p>  
This function fills **`m`** (**[map_node](data_processing.md#linkstream-types)**) by processing recursively SVD on the adjacency matrix of **`g`** (**[graph](data_processing.md#linkstream-types)**).

<p class="function">int fct_map_edge (int i, int j)</p>  
This mapping function returns an unique edge id (|VxV|) by giving 2 nodes labels (i->j).

<hr/>

## Transform graph

<p class="function">int max_lvl_decomposition (graph *g)</p>  
<p class="function">int max_lvl_decomposition (map_edge *m)</p>  
These functions return the maximum level of  wavelet decomposition you can compute.


<p class="function">int wavelet_transform_by_level (const gsl_wavelet * w, gsl_spmatrix *data, size_t n, gsl_wavelet_direction dir, sp_workspace * work, uint nb_level = 1, uint current_lvl = 0)</p>   
This function computes **`nb_level`** haar wavelet transform from **`current_lvl`** in the direction **`dir`**. A sparse matrix of data and a sparse worspace **`work`** with a size **`n`** has to be provided. To perform wavelet transform, gsl_wavelet **`w`** has be to initialized as gsl_wavelet_haar and k = 2.


<p class="function">int decomposition_BFS (graph *g, map_edge *m, uint level, gsl_matrix *A, gsl_matrix *D)</p>  
<p class="function">int decomposition_SVD (graph *g, map_node *m, uint level, gsl_matrix *A, gsl_matrix *D)</p>  
These functions compute decomposition of graph with haar wavelet transform, following the indicated mapping (**[BFS](#bfs)** or **[SVD](#svd)**) and fill matrix A (resp. D) with the accuracy (resp. difference) part. Memory for matrices A et D has to be allocated with appropriate sizes (see functions : **[size1A](#decomposition-result-matrix)**,  **[size1D](#decomposition-result-matrix)**).

<hr/>

## Transform linkstream

<p class="function">int aggregate_link_stream (float_lien *fl, graph *g)</p>
This function fills **`g`** (**[graph](data_processing#linkstream-types)**) as the aggregation of graphs for all time in **`fl`** (**[float_lien](data_processing#linkstream-types)**).  
Memory has to be allocated for **`g`**. (Use for graph decomposition)


<p class="function">int graph_sequence_decomposition_BFS (graph_set *gs, map_edge *m, uint level, gsl_matrix *A, gsl_spmatrix *D)</p>  
This function computes the decomposition of graph sequence with haar wavelet transform, following BFS mapping **`m`** and fill matrix A (resp. D) with the accuracy (resp. difference) part. Memory for matrices A et D has to be allocated with size2 equal to number of graph in the sequence and size1 with the appropriate value (see functions : **[size1A](#decomposition-result-matrix)**,  **[size1D](#decomposition-result-matrix)**).
<hr/>

## Example

** Wavelet Decomposition with SVD mapping : **
```c
// example.cpp
#include "lib/data_processing.hpp"
#include "lib/decomposition.h"

int main(int argc, char **argv)
{
    /* Initialisation */
    (void)(argc); /* avoid unused parameter warning */
    std::string path = argv[1];
    float_lien *fl = float_lien_alloc();
    map_node *m = map_node_alloc();
    graph *g = new graph();
    gsl_matrix *A;
    gsl_spmatrix *D;

    std::cout << "===== Data processing =====\n";
    read_link_stream(path, fl);
    std::cout << "===== End Data processing =====\n";

    /* Aggregation of graphs for all time in fl */
    aggregate_link_stream(fl, g);

    std::cout << "===== MAPPING =====\n";
    mapping_SVD(g, m);
    std::cout << "===== End Mapping =====\n";

    std::cout << "===== Decomposition =====\n";
    uint level = max_lvl_decomposition(g);
    A = matrixA_alloc(g, level);
    D = matrixD_alloc(g, level);
    decomposition_SVD(g, m, level, A, D);
    std::cout << "===== FIN Decomposition =====\n";

    /* Release memory */
    matrixA_delete(A);
    matrixD_delete(D);
    float_lien_delete(fl);
    map_node_delet(m);
    delete g;

    return 0;
}
```

Compile with :  
``g++ ...``

Run with :  
``./example relative_path/file_linkstream.txt``