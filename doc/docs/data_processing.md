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

# Data Processing

This chapter describes function to process data into float lien and mapping structure. 
The data processing functions are declared in the header file **`data_processing.hpp`**.
Linkstream types are declared in the header file **`linkstream_types.hpp`**.
<hr/>

## Linkstream types
| **Prefix**  | **Type**  | 
| :--------------- | :--------------- |
| time_id | **`int`** |
| node | **`int`** |
| weigth | **`double`** |
| edge | **`std::pair<node,node>`** |
| weighted_graph | **`std::vector< std::pair<edge,weigth> >`** |
| graph | **`std::vector<edge>`** |
| time_edgelist | **`std::map< time_id, weigthed_graph >`** |
| node_set | **`std::set<node>`** |
| map_edge | **`std::map<edge, int>`** |
| map_node | **`std::map<node, int>`** |
| adjacency_list | **`std::map< node, std::vector<node> >`** |
| graph_set | **`std::vector<graph>`** |

<p class="function">Class float_lien</p>
This class contain informations relative to float lien. It can be used to load a float lien from a file (see **[read_link_stream](#functions)**).
Data members :  
- time_graph (**[time_edgelist](#initialization)**) : map where each time refers to a list(graph) of pairs with an edge and its weigth  

<hr/>

## Initialization Functions

<p class="function">float_lien * float_lien_alloc ( )</p>
This function allocates memory for an object of class float_lien.

<p class="function">void float_lien_delete (float_lien* FL)</p>
This function releases allocated memory of the float_lien class object.

<hr/>

## Functions

<p class="function">int read_link_stream (std::string edgelist_path, float_lien *float_lien_to_fill, char delimiter = ' , ')</p>
This function read the file located at **`edgelist_path`** and load information into **`float_lien_to_fill`**. You can add a fourth column to specify the weigth of each edge otherwise it will be considered as one. You can change the **`delimiter`** between columns (default = ' , ').  
The text format has to be :  
**nodeA, nodeB, timeX  
nodeC, nodeD, timeU  
...&emsp;&emsp;,&emsp;&emsp;...&emsp;&emsp;,&emsp;&emsp;...   
nodeY, nodeZ, timeV**

<p class="function">int weigthedGraph_to_Graph (weigthed_graph * wg, graph *g)</p>
This function fills **`g`** (**[graph](#linkstream-types)**) with edges in **`wg`** (**[weighted_graph](#linkstream-types)**).

<hr/>

## Example

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

Compile with :  
``g++ ...``

Run with :  
``./example relative_path/file_linkstream.txt``