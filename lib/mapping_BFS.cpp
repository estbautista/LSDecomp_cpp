/* mapping_BFS.cpp
 * 
 * Copyright (C) 2022 Bastien Guillemare
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "mapping_BFS.h"
//#include <chrono> // For test purpose

map_edge *map_edge_alloc()
{
    map_edge *m = new map_edge();
    if (m == NULL)
    {
        GSL_ERROR_VAL("failed to allocate space for map_edge", GSL_ENOMEM, 0);
    }

    return (m);
}

void map_edge_delete(map_edge *m)
{
    if (m != NULL)
    {
        delete m;
    }
}

/* Return if x is in the array */
int is_in(std::vector<node> *array, node x)
{
    for (std::vector<node>::const_iterator it = array->begin(); it != array->end(); ++it)
    {
        if ((*it) == x)
        {
            return 1;
        }
    }

    return 0;
}

/* Print function for debugging purpose */
void print_graph_set(graph_set *gs)
{
    std::cout << "Graph_set :\n";

    for (graph_set::const_iterator g = gs->begin(); g != gs->end(); ++g)
    {
        std::cout << "{ ";
        int cpt = 1;
        for (graph::const_iterator e = g->begin(); e != g->end(); ++e)
        {
            std::cout << "(" << e->first << "," << e->second << ")\t";
            if (cpt % 5 == 0)
            {
                std::cout << std::endl;
            }
            cpt++;
        }
        std::cout << "}";
        if (cpt % 5 != 1)
        {
            std::cout << std::endl;
        }
    }

    std::cout << std::endl;
}

graph_set *rec_BFS(graph_set *gset)
{
    // std::cout << "BFS size:" << gset->begin()->size() << std::endl;
    // print_graph_set(gset);
    graph_set *partition = new graph_set();

    for (graph_set::const_iterator it = gset->begin(); it != gset->end(); ++it)
    {
        std::set<node> visited;
        std::queue<node> file;
        size_t h_size, size = 0;
        size_t cpt = 0;

        // init adjacency list

        std::vector<node> node_order; // keep node order of graph because map are reordered in tree search
        node_order.reserve(it->size());
        adjacency_list adj_list;
        adjacency_list::iterator item;
        for (graph::const_iterator e = it->begin(); e != it->end(); ++e)
        {
            item = adj_list.find(e->first);
            if (item != adj_list.end())
            {
                if (!is_in(&(item->second), e->second))
                {
                    item->second.push_back(e->second);
                    size++;
                }
            }
            else
            {
                adj_list.insert({e->first, {e->second}});
                node_order.push_back(e->first);
                size++;
            }
        }
        h_size = is_power_of_two(size) ? size / 2 : next_power_of_2(size) / 2;

        std::set<node>::iterator is_visited;
        while (cpt != size)
        {

            if (file.empty())
            {
                // Fill "file" with the first element in the graph that we have not visited yet
                for (std::vector<node>::const_iterator n = node_order.begin(); n != node_order.end(); ++n)
                {
                    is_visited = visited.find((*n));
                    if (is_visited == visited.end())
                    { // if not visited
                        file.push((*n));
                        node_order.erase(n);
                        break;
                    }
                }
            }

            // Visit the first element in the "file" and add it to "visited"
            node first = file.front();
            is_visited = visited.find(first);
            if (is_visited == visited.end())
            {
                item = adj_list.find(first);
                if (item != adj_list.end())
                {
                    for (std::vector<node>::const_iterator succ = item->second.begin(); succ != item->second.end(); ++succ)
                    {
                        if (cpt == 0 || cpt == h_size)
                        {
                            partition->push_back({{first, (*succ)}});
                        }
                        else
                        {
                            partition->back().push_back({first, (*succ)});
                        }
                        cpt++;
                        file.push((*succ));
                    }
                    visited.insert(first);
                    // adj_list.erase(item);
                }
            }
            file.pop();
        }
    }

    delete gset;
    if (partition->begin()->size() > 1)
    {
        partition = rec_BFS(partition);
    }

    return partition;
}

int mapping_BFS(graph *g, map_edge *m)
{
    // Initialisation graph to graph_set
    graph_set *gset = new graph_set({(*g)});
    if (gset == NULL)
    {
        GSL_ERROR_VAL("Failed to allocate memory for graph_set", GSL_ENOMEM, 0);
    }

    // BFS
    gset = rec_BFS(gset);

    // copy in map edge
    int id = 1;
    for (graph_set::const_iterator it = gset->begin(); it != gset->end(); ++it)
    {
        for (graph::const_iterator e = it->begin(); e != it->end(); ++e)
        {
            m->insert({(*e), id});
            id++;
        }
    }

    delete gset;
    return EXIT_SUCCESS;
}