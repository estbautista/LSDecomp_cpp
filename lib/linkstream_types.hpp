/* linkstream_types.hpp
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
#ifndef __LINKSTREAM_TYPE_H
#define __LINKSTREAM_TYPE_H
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <gsl/gsl_spmatrix.h>

typedef int time_id;
typedef int node;
typedef double weigth;
typedef std::pair<node, node> edge;                          // pair of nodes
typedef std::vector<std::pair<edge, weigth>> weigthed_graph; // list of pairs with an edge and its weigth
typedef std::vector<edge> graph;                             // list of edges
typedef std::map<time_id, weigthed_graph> time_edgelist;     // map where each time refers to a graph
typedef std::set<node> node_set;                          // list of nodes
typedef std::map<edge, int> map_edge;                        // Map : return an unique id for each edge in it
typedef std::map<node, int> map_node;                        // Map : return an unique id for each node in it
typedef std::map<node, std::vector<node>> adjacency_list;    // Adjacency list of each node of graph
typedef std::vector<graph> graph_set;                        // List of graph

class float_lien
{
public:
  time_edgelist time_graph; // Map : list edges+weigth by time
  //node_set V;               // List all nodes
};

typedef struct
{
  gsl_spmatrix *data; // mapping data
  size_t size1;       // data dim1 size
  size_t size2;       // data dim2 size
} mapping_wavelet;

/* sparse workspace for wavelet transform */
typedef struct
{
  gsl_spmatrix *m;
  size_t size1;
  size_t size2;
} sp_workspace;

#endif /* __LINKSTREAM_TYPE_H */