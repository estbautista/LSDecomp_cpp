# Examples

This chapter shows examples of program to used the library.
<hr/>

## Data processing

```c
//example.cpp
#include "lib/data_processing.hpp"

int main(int argc, char **argv){
    /* Initialization */
    (void)(argc); /* avoid unused parameter warning */
    std::string path = argv[1];
    float_lien *fl = float_lien_alloc();
    read_link_stream(path,fl); 
    
    /* Your code */
    int print_limit = 100;
    
    for (time_edgelist::const_iterator t_g = fl->begin(); t_g != fl->end(); ++t_g)
    {
        int cpt = 1;
        std::cout << "Time = " << t_g->first << " :\n";
        for(weigthed_graph::const_iterator we = t_g->second.begin(); we != t_g->second.end(); ++we)
        {
            std::cout << "(" << we->first.first << "," << we->first.second << "):" << we->second << "   ";
            if (cpt % 5 == 0)
                std::cout << std::endl;
            
            if (cpt % print_limit == 0)
                break;
            cpt++;
        }
        if (cpt % 5 != 1)
            std::cout << std::endl;   
    }

    /* Release memory */
    float_lien_delete(fl);

    return 0;
}
```

## Decomposition

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