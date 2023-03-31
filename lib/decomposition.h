/* decomposition.h
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
#ifndef __DECOMPOSITION_H
#define __DECOMPOSITION_H
#include <gsl/gsl_errno.h>
#include <gsl/gsl_wavelet.h>

#include "mapping_BFS.h"
#include "mapping_SVD.h"
#include "linkstream_types.hpp"
#include "fct_utils.hpp"

/* Allocate memory for structure variable of sp_workspace */
sp_workspace *sp_workspace_alloc(size_t size1, size_t size2 = 1);

/* Release sparse workspace memory */
void sp_workspace_delete(sp_workspace *work);

/* Wavelet transform of data */
int wavelet_transform_by_level(const gsl_wavelet *w,
                               gsl_spmatrix *data, size_t n,
                               gsl_wavelet_direction dir, sp_workspace *work,
                               uint nb_level = 1, uint current_lvl = 0);

/* Return size1 of A (Accuracy dense matrix) for allocation */
size_t size1A(size_t n, uint level);

/* Return size1 of D (Difference sparse matrix) for allocation () */
size_t size1D(size_t n, uint level);

/* Memory allocation for dense matrix A.
(use to store accuracy part of level-wavelet transform) */
gsl_matrix *matrixA_alloc(size_t dim1, size_t dim2, uint level);
gsl_matrix *matrixA_alloc(map_edge *m, uint level, size_t n_graphs = 1);
gsl_matrix *matrixA_alloc(graph *g, uint level);

/* Memory allocation for sparse matrix D.
(use to store difference part of level-wavelet transform) */
gsl_spmatrix *matrixD_alloc(size_t dim1, size_t dim2, uint level);
gsl_spmatrix *matrixD_alloc(map_edge *m, uint level, size_t n_graphs = 1);
gsl_spmatrix *matrixD_alloc(graph *g, uint level);

/* Release memory for dense matrix A */
void matrixA_delete(gsl_matrix *A);

/* Release memory for sparse matrix D */
void matrixD_delete(gsl_spmatrix *D);

/* Return the maximum level of decomposition you can compute */
int max_lvl_decomposition(graph *g);
int max_lvl_decomposition(map_edge *m);

/* Wavelet decomposition of graph with different mapping
  Store results in A (Accuracy dense matrix) and D (Difference sparse matrix) */
int decomposition_BFS(graph *g, map_edge *m, uint level, gsl_matrix *A, gsl_spmatrix *D);
int decomposition_SVD(graph *g, map_node *m, uint level, gsl_matrix *A, gsl_spmatrix *D);

/* Aggregation of graphs for all time in fl */
int aggregate_link_stream (float_lien *fl, graph *g);

/* Wavelet decomposition of graph_sequence
Store results in A (Accuracy dense matrix) and D (Difference sparse matrix) */
int graph_sequence_decomposition_BFS(graph_set *gs, map_edge *m, uint level, gsl_matrix *A, gsl_spmatrix *D);

/*===== Print functions for debugging purposes =====*/

void print_sp_mat(arma::sp_mat *M);

void print_graph(graph *g, int limit);

void print_weigthed_graph(weigthed_graph *g, int limit);

void print_map_edge(map_edge *mg, int limit);

void print_map_node(map_node *m, int limit);

#endif /* __DECOMPOSITION_H */