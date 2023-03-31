/* decomposition.cpp
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
#include "decomposition.h"
#include <chrono> // For test purpose

sp_workspace *sp_workspace_alloc(size_t size1, size_t size2)
{
  sp_workspace *work = new sp_workspace;

  if (work == NULL)
  {
    GSL_ERROR_VAL("failed to allocate struct", GSL_ENOMEM, 0);
  }

  work->size1 = size1;
  work->size2 = size2;
  work->m = gsl_spmatrix_alloc(work->size1, work->size2);

  if (work->m == NULL)
  {
    /* error in constructor, prevent memory leak */
    free(work);
    GSL_ERROR_VAL("failed to allocate space", GSL_ENOMEM, 0);
  }
  return work;
}

void sp_workspace_delete(sp_workspace *work)
{
  if (work != NULL)
  {
    /* release space */
    if(work->m != NULL)
      gsl_spmatrix_free(work->m);
    work->m = NULL;
    free(work);
  }
}

// jth non zero of m
#define ELEMENT(m, j) (gsl_spmatrix_get(m, m->i[j], m->p[j]))
// Add x into sparse matrix m[i,j]
#define ADD(m, i, j, x) (gsl_spmatrix_set(m, i, j, x + gsl_spmatrix_get(m, i, j)))

/* Return log2(n) or -1 if n is not pair*/
static int binary_logn(const size_t n)
{
  size_t ntest;
  size_t logn = 0;
  size_t k = 1;

  while (k < n)
  {
    k *= 2;
    logn++;
  }

  ntest = ((size_t)1 << logn);

  if (n != ntest)
  {
    return -1; /* n is not a power of 2 */
  }

  return logn;
}

/* 1 level decomposition for index 0 to size-1 in data */
static void dwt_step(const gsl_wavelet *w, gsl_spmatrix *data, size_t n,
                     gsl_wavelet_direction dir, sp_workspace *work)
{
  int nh;
  size_t j;

  /* size_t nmod = w->nc * n;
  nmod -= w->offset;            // center support
  */

  int size = n; // To avoid : comparison of integer expressions of different signedness
  nh = n >> 1;

  if (work->m->nz != 0)
    gsl_spmatrix_set_zero(work->m);

  if (dir == gsl_wavelet_forward)
  {
    for (j = 0; j < data->nz; j++)
    {
      if (data->i[j] < size)
      {
        if (data->i[j] % 2 == 0)
        {
          ADD(work->m, data->i[j] / 2, data->p[j], w->h1[0] * ELEMENT(data, j));
          ADD(work->m, data->i[j] / 2 + nh, data->p[j], w->g1[0] * ELEMENT(data, j));
        }
        else
        {
          ADD(work->m, data->i[j] / 2, data->p[j], w->h1[1] * ELEMENT(data, j));
          ADD(work->m, data->i[j] / 2 + nh, data->p[j], w->g1[1] * ELEMENT(data, j));
        }
      }
      else
      {
        gsl_spmatrix_set(work->m, data->i[j], data->p[j], ELEMENT(data, j));
      }
    }
  }
  else
  {
    for (j = 0; j < data->nz; j++)
    {
      if (data->i[j] < size)
      {
        if (data->i[j] < nh)
        {
          ADD(work->m, 2 * data->i[j], data->p[j], w->h2[0] * ELEMENT(data, j));
          ADD(work->m, 2 * data->i[j] + 1, data->p[j], w->h2[1] * ELEMENT(data, j));
        }
        else
        {
          ADD(work->m, 2 * (data->i[j] - nh), data->p[j], w->g2[0] * ELEMENT(data, j));
          ADD(work->m, 2 * (data->i[j] - nh) + 1, data->p[j], w->g2[1] * ELEMENT(data, j));
        }
      }
      else
      {
        gsl_spmatrix_set(work->m, data->i[j], data->p[j], ELEMENT(data, j));
      }
    }
  }

  gsl_spmatrix_set_zero(data);
  for (j = 0; j < work->m->nz; j++)
  {
    gsl_spmatrix_set(data, work->m->i[j], data->p[j], ELEMENT(work->m, j));
  }
  gsl_spmatrix_set_zero(work->m);
}

size_t size1A(size_t n, uint level)
{
  return n * pow(0.5, level);
}

size_t size1D(size_t n, uint level)
{
  return n;
}

gsl_matrix *matrixA_alloc(size_t dim1, size_t dim2, uint level)
{
  if (dim1 < 1 || dim2 < 1)
  {
    GSL_ERROR_VAL("dimensions has to be positive integer", GSL_EBADLEN, 0);
  }
  gsl_matrix *A = gsl_matrix_alloc(size1A(dim1, level), dim2);

  if (A == NULL)
  {
    GSL_ERROR_VAL("failed to allocate space for dense matrix A", GSL_ENOMEM, 0);
  }

  return A;
}

gsl_matrix *matrixA_alloc(map_edge *m, uint level, size_t n_graphs)
{
  if (n_graphs < 1)
  {
    GSL_ERROR_VAL("n_graphs has to be positive integer", GSL_EBADLEN, 0);
  }
  size_t size = m->size(); // |E|
  size = is_power_of_two(size) ? size : next_power_of_2(size);
  gsl_matrix *A = gsl_matrix_alloc(size1A(size, level), n_graphs);

  if (A == NULL)
  {
    GSL_ERROR_VAL("failed to allocate space for dense matrix A", GSL_ENOMEM, 0);
  }

  return A;
}

gsl_matrix *matrixA_alloc(graph *g, uint level)
{
  std::set<edge> edge_set;
  for (graph::const_iterator it = g->begin(); it != g->end(); it++)
  {
    edge_set.insert((*it));
  }
  size_t size = edge_set.size() * edge_set.size(); //  |VxV|
  size = is_power_of_two(size) ? size : next_power_of_2(size);
  gsl_matrix *A = gsl_matrix_alloc(size1A(size, level), 1);

  if (A == NULL)
  {
    GSL_ERROR_VAL("failed to allocate space for dense matrix A", GSL_ENOMEM, 0);
  }

  return A;
}

gsl_spmatrix *matrixD_alloc(size_t dim1, size_t dim2, uint level)
{
  if (dim1 < 1 || dim2 < 1)
  {
    GSL_ERROR_VAL("dimensions has be positive integer", GSL_EBADLEN, 0);
  }
  gsl_spmatrix *D = gsl_spmatrix_alloc(size1D(dim1, level), dim2);

  if (D == NULL)
  {
    GSL_ERROR_VAL("failed to allocate space for dense matrix A", GSL_ENOMEM, 0);
  }

  return D;
}

gsl_spmatrix *matrixD_alloc(map_edge *m, uint level, size_t n_graphs)
{
  if (n_graphs < 1)
  {
    GSL_ERROR_VAL("n_graphs has be positive integer", GSL_EBADLEN, 0);
  }
  size_t size = m->size();// |E|
  size = is_power_of_two(size) ? size : next_power_of_2(size);
  gsl_spmatrix *D = gsl_spmatrix_alloc(size1D(size, level), n_graphs);

  if (D == NULL)
  {
    GSL_ERROR_VAL("failed to allocate space for dense matrix A", GSL_ENOMEM, 0);
  }

  return D;
}

gsl_spmatrix *matrixD_alloc(graph *g, uint level)
{
  std::set<edge> edge_set;
  for (graph::const_iterator it = g->begin(); it != g->end(); it++)
  {
    edge_set.insert((*it));
  }
  size_t size = edge_set.size() * edge_set.size(); // |VxV|
  size = is_power_of_two(size) ? size : next_power_of_2(size);
  gsl_spmatrix *D = gsl_spmatrix_alloc(size1D(size, level), 1);

  if (D == NULL)
  {
    GSL_ERROR_VAL("failed to allocate space for dense matrix A", GSL_ENOMEM, 0);
  }

  return D;
}

void matrixA_delete(gsl_matrix *A)
{
  gsl_matrix_free(A);
}

void matrixD_delete(gsl_spmatrix *D)
{
  gsl_spmatrix_free(D);
}

int max_lvl_decomposition(graph *g)
{
  std::set<edge> edge_set;
  for (graph::const_iterator it = g->begin(); it != g->end(); it++)
  {
    edge_set.insert((*it));
  }
  size_t size = edge_set.size() * edge_set.size(); //  |VxV|
  size = is_power_of_two(size) ? size : next_power_of_2(size);
  return binary_logn(edge_set.size());
}

int max_lvl_decomposition(map_edge *m)
{
  size_t size = m->size();
  size = is_power_of_two(size) ? size : next_power_of_2(size);
  return binary_logn(m->size());
}

int wavelet_transform_by_level(const gsl_wavelet *w,
                               gsl_spmatrix *data, size_t n,
                               gsl_wavelet_direction dir, sp_workspace *work,
                               uint nb_level, uint current_lvl)
{
  size_t i;

  if (work->size1 < n)
  {
    GSL_ERROR("not enough workspace provided", GSL_EINVAL);
  }

  if (binary_logn(n) == -1)
  {
    GSL_ERROR("n is not a power of 2", GSL_EINVAL);
  }

  if (dir == gsl_wavelet_forward)
  {
    for (i = n * pow(0.5, current_lvl); (i > n * pow(0.5, current_lvl + nb_level)) && (i >= 2); i >>= 1)
    {
      dwt_step(w, data, i, dir, work);
    }
  }
  else
  {
    for (i = n * pow(0.5, current_lvl - 1); (i < n * pow(0.5, current_lvl - 1 - nb_level)) && (i <= n); i <<= 1)
    {
      dwt_step(w, data, i, dir, work);
    }
  }

  return GSL_SUCCESS;
}

int decomposition_BFS(graph *g, map_edge *m, uint level, gsl_matrix *A, gsl_spmatrix *D)
{
  if (level == 0)
  {
    GSL_ERROR_VAL("level must be greater than zero", GSL_EINVAL, 0);
  }

  size_t size = m->size(); // |E|
  size = is_power_of_two(size) ? size : next_power_of_2(size);

  if ((int)level > binary_logn(size))
  {
    GSL_ERROR_VAL("level must be less than log2(number of e)", GSL_EINVAL, 0);
  }

  if (A == NULL || D == NULL)
  {
    GSL_ERROR_VAL("Memory has to be allocated for A and D", GSL_EFAULT, 0);
  }

  if (A->size1 < size1A(size, level))
  {
    GSL_ERROR_VAL("A size is not conform", GSL_EBADLEN, 0);
  }

  if (D->size1 < size1D(size, level))
  {
    GSL_ERROR_VAL("D size is not conform", GSL_EBADLEN, 0);
  }

  gsl_wavelet *w = gsl_wavelet_alloc(gsl_wavelet_haar, 2);
  sp_workspace *work = sp_workspace_alloc(size);
  gsl_spmatrix *data = gsl_spmatrix_alloc(size, 1);

  for (graph::const_iterator e = g->begin(); e != g->end(); ++e)
  {
    ADD(data, (*m)[(*e)], 0, 1);
  }

  wavelet_transform_by_level(w, data, size, gsl_wavelet_forward, work, level, 0);

  for (size_t j = 0; j < data->nz; j++)
  {
    if (data->i[j] < (int)A->size1)
      gsl_matrix_set(A, data->i[j], 0, gsl_spmatrix_get(data, data->i[j], 0));
    else
      gsl_spmatrix_set(D, data->i[j], 0, gsl_spmatrix_get(data, data->i[j], 0));
  }

  sp_workspace_delete(work);
  gsl_wavelet_free(w);
  return GSL_SUCCESS;
}

int decomposition_SVD(graph *g, map_node *m, uint level, gsl_matrix *A, gsl_spmatrix *D)
{
  if (level == 0)
  {
    GSL_ERROR_VAL("level must be greater than zero", GSL_EINVAL, 0);
  }
  if (A == NULL || D == NULL)
  {
    GSL_ERROR_VAL("Memory has to be allocated for A and D", GSL_EFAULT, 0);
  }

  std::set<edge> edge_set;
  for (graph::const_iterator ed = g->begin(); ed != g->end(); ++ed)
  {
    edge_set.insert((*ed));
  }
  size_t size = edge_set.size() * edge_set.size(); // |VxV|
  size = is_power_of_two(size) ? size : next_power_of_2(size);

  if ((int)level > binary_logn(size))
  {
    GSL_ERROR_VAL("level must be less than or equal to log2(number of e)", GSL_EINVAL, 0);
  }

  if (A->size1 < size1A(size, level))
  {
    GSL_ERROR_VAL("A size is not conform", GSL_EBADLEN, 0);
  }

  if (D->size1 < size1D(size, level))
  {
    GSL_ERROR_VAL("D size is not conform", GSL_EBADLEN, 0);
  }

  gsl_wavelet *w = gsl_wavelet_alloc(gsl_wavelet_haar, 2);
  sp_workspace *work = sp_workspace_alloc(size);
  gsl_spmatrix *data = gsl_spmatrix_alloc(size, 1);
  for (graph::const_iterator ed = g->begin(); ed != g->end(); ++ed)
  {
    ADD(data, fct_map_edge((*m)[ed->first], (*m)[ed->second]) - 1, 0, 1);
  }

  wavelet_transform_by_level(w, data, size, gsl_wavelet_forward, work, level, 0);
  
  for (size_t j = 0; j < data->nz; j++)
  {
    if (data->i[j] < (int)A->size1)
      gsl_matrix_set(A, data->i[j], 0, gsl_spmatrix_get(data, data->i[j], 0));
    else
      gsl_spmatrix_set(D, data->i[j], 0, gsl_spmatrix_get(data, data->i[j], 0));
  }

  sp_workspace_delete(work);
  gsl_wavelet_free(w);
  return GSL_SUCCESS;
}

int aggregate_link_stream(float_lien *fl, graph *g)
{
    if (fl == NULL)
    {
        GSL_ERROR_VAL("float_lien is NULL", GSL_EFAULT, 0);
    }
    if (g == NULL)
    {
        GSL_ERROR_VAL("Memory has to be allocated for graph g", GSL_EFAULT, 0);
    }

    for (time_edgelist::const_iterator it = fl->time_graph.begin(); it != fl->time_graph.end(); ++it)
    {
        for (weigthed_graph::const_iterator e_w = it->second.begin(); e_w != it->second.end(); ++e_w)
        {
            g->push_back(e_w->first);
        }
    }

    return EXIT_SUCCESS;
}

int graph_sequence_decomposition_BFS(graph_set *gs, map_edge *m, uint level, gsl_matrix *A, gsl_spmatrix *D)
{
  if (level == 0)
  {
    GSL_ERROR_VAL("level must be greater than zero", GSL_EINVAL, 0);
  }

  size_t size = m->size();
  size = is_power_of_two(size) ? size : next_power_of_2(size);

  if ((int)level > binary_logn(size))
  {
    GSL_ERROR_VAL("level must be less than or equal to log2(number of e)", GSL_EINVAL, 0);
  }

  if (A->size1 < size1A(size, level) || A->size2 != gs->size())
  {
    GSL_ERROR_VAL("A size is not conform", GSL_EBADLEN, 0);
  }

  if (D->size1 < size1D(size, level) || D->size2 != gs->size())
  {
    GSL_ERROR_VAL("D size is not conform", GSL_EBADLEN, 0);
  }

  gsl_wavelet *w = gsl_wavelet_alloc(gsl_wavelet_haar, 2);
  sp_workspace *work = sp_workspace_alloc(size, gs->size());
  gsl_spmatrix *data = gsl_spmatrix_alloc(size, gs->size());

  int num_graph = 0;
  for (graph_set::const_iterator g = gs->begin(); g != gs->end(); ++g)
  {
    for (graph::const_iterator e = g->begin(); e != g->end(); ++e)
    {
      ADD(data, (*m)[(*e)], num_graph, 1);
    }
    num_graph++;
  }

  wavelet_transform_by_level(w, data, size, gsl_wavelet_forward, work, level, 0);

  for (size_t j = 0; j < data->nz; j++)
  {
    if (data->i[j] < (int)A->size1)
      gsl_matrix_set(A, data->i[j], data->p[j], gsl_spmatrix_get(data, data->i[j], data->p[j]));
    else
      gsl_spmatrix_set(D, data->i[j], data->p[j], gsl_spmatrix_get(data, data->i[j], data->p[j]));
  }

  sp_workspace_delete(work);
  gsl_wavelet_free(w);
  return GSL_SUCCESS;
}

/*===== Debugging functions =====*/

void print_sp_mat(arma::sp_mat *M)
{
  for (arma::u64 row = 0; row < M->n_rows; row++)
  {
    for (arma::u64 col = 0; col < M->n_cols - 1; col++)
    {
      std::cout << (*M)(row, col) << " ";
    }
    std::cout << (*M)(row, M->n_cols - 1) << "\n";
  }
}

void print_graph(graph *g, int limit)
{
  std::cout << "Graph :\n";
  int cpt = 1;
  for (graph::const_iterator e = g->begin(); e != g->end(); e++)
  {
    std::cout << "(" << e->first << "," << e->second << ")    ";
    if (cpt % 5 == 0)
    {
      std::cout << std::endl;
    }
    if (cpt % limit == 0)
      break;
    cpt++;
  }
  if (cpt % 5 != 1)
  {
    std::cout << std::endl;
  }
}

void print_weigthed_graph(weigthed_graph *g, int limit)
{
  std::cout << "Graph : (ni,nj):w\n";
  int cpt = 1;
  for (weigthed_graph::const_iterator we = g->begin(); we != g->end(); we++)
  {
    std::cout << "(" << we->first.first << "," << we->first.second << "):" << we->second << "   ";
    if (cpt % 5 == 0)
    {
      std::cout << std::endl;
    }
    if (cpt % limit == 0)
      break;
    cpt++;
  }
  if (cpt % 5 != 1)
  {
    std::cout << std::endl;
  }
}

void print_map_edge(map_edge *mg, int limit)
{
  std::cout << "Map graph : (ni,nj):id\n";
  int cpt = 1;
  for (map_edge::const_iterator e_id = mg->begin(); e_id != mg->end(); e_id++)
  {
    std::cout << "(" << e_id->first.first << "," << e_id->first.second << "):" << e_id->second << "   ";
    if (cpt % 5 == 0)
    {
      std::cout << std::endl;
    }
    if (cpt % limit == 0)
      break;
    cpt++;
  }
  if (cpt % 5 != 1)
  {
    std::cout << std::endl;
  }
}

void print_map_node(map_node *m, int limit)
{
  std::cout << "Map graph : (ni):id\n";
  int cpt = 1;
  for (map_node::const_iterator n_id = m->begin(); n_id != m->end(); n_id++)
  {
    std::cout << "(" << n_id->first << "):" << n_id->second << "    ";
    if (cpt % 5 == 0)
    {
      std::cout << std::endl;
    }
    if (cpt % limit == 0)
      break;
    cpt++;
  }
  if (cpt % 5 != 1)
  {
    std::cout << std::endl;
  }
}
