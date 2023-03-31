/* mapping_SVD.h
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
#ifndef __MAPPING_SVD_H
#define __MAPPING_SVD_H

#include "linkstream_types.hpp"
#include "fct_utils.hpp"
#include <armadillo>

/* Allocate memory for structure variable of map_node */
map_node *map_node_alloc();

/* Release allocated memory of map_node struct variable */
void map_node_delete(map_node *m);

/* Create an adjacency matrix of edges of the graph and
apply SVD (Singular Value Decomposition) recursively to rank all nodes. 
Fill map_node 'm' with ids of each node in graph 'g' */
int mapping_SVD(graph *g, map_node *m);

/* Mapping function : return an unique edge id (|VxV|) by giving 2 nodes id's */
int fct_map_edge(int i, int j);

#endif /* __MAPPING_SVD_H */