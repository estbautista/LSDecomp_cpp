/* mapping_SVD.cpp
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
#include "mapping_SVD.h"
#include <chrono> // For test purpose

map_node *map_node_alloc()
{
  map_node *m = new map_node();
  if (m == NULL)
  {
    GSL_ERROR_VAL("failed to allocate space for map_node", GSL_ENOMEM, 0);
  }

  return (m);
}

void map_node_delete(map_node *m)
{
  if (m != NULL)
  {
    delete m;
  }
}

void rec_SVD(arma::sp_mat *adj_Matrix, arma::u64 n_rows, uint k, arma::uvec *alphak)
{ 
  std::chrono::duration<double, std::milli> interval_SVD;
  auto t1 = std::chrono::high_resolution_clock::now();
  // Init for svd
  arma::mat U;
  arma::vec s;
  arma::mat V;
  size_t tol = 1e-4;
  // SVD : adj_Matrix = U*D*V and s = diag(D)
  arma::svds(U, s, V, (*adj_Matrix), k,tol);

  auto t2 = std::chrono::high_resolution_clock::now();
  interval_SVD = t2 - t1;
  
  if (n_rows < 3) // Base case
  {
    std::chrono::duration<double, std::milli> interval_tri;
    t1 = std::chrono::high_resolution_clock::now();
    // sort value : alphak[rank SVD] = ID
    (*alphak) = arma::sort_index(U.col(1).rows(0, n_rows - 1), "descend").as_col();
    t2 = std::chrono::high_resolution_clock::now();
    interval_tri = t2 - t1;
    /* std::cout << n_rows <<" " << interval_SVD.count() / 1000 << " ";
    std::cout << interval_tri.count() / 1000 << "\n"; */
    
  }
  else if (n_rows > 3)
  {
    std::chrono::duration<double, std::milli> interval_tri;
    t1 = std::chrono::high_resolution_clock::now();
    // sort value : alphak[rank SVD] = ID
    arma::uvec rank = arma::sort_index(U.col(1).rows(0, n_rows - 1), "descend").as_col();
    t2 = std::chrono::high_resolution_clock::now();
    interval_tri = t2 - t1;
    

    std::chrono::duration<double, std::milli> interval_split;
    t1 = std::chrono::high_resolution_clock::now();
    // Split in 2 parts : init
    arma::u64 h_size = adj_Matrix->n_rows / 2;
    arma::sp_mat *part1 = new arma::sp_mat(h_size, adj_Matrix->n_cols);
    int size_part2 = n_rows - h_size;
    size_part2 = is_power_of_two(size_part2) ? size_part2 : next_power_of_2(size_part2);
    arma::sp_mat *part2 = new arma::sp_mat(size_part2, adj_Matrix->n_cols);

    // Compute row order : order[id] = SVD rank
    arma::uvec order(n_rows); 
    for (arma::u64 row = 0; row < n_rows; ++row)
    {
      order(rank(row)) = row;
    }

    // Fill part following row order
    arma::sp_mat::const_iterator it_end = adj_Matrix->end();
    for (arma::sp_mat::const_iterator it = adj_Matrix->begin(); it != it_end; ++it)
    {
      if (order[it.row()] < h_size)
      { // part 1
        (*part1)(order[it.row()], it.col()) = (*it);
      }
      else
      { // part 2
        (*part2)(order[it.row()] - h_size, it.col()) = (*it);
      }
    }
    t2 = std::chrono::high_resolution_clock::now();
    interval_split = t2 - t1;

    // rec_SVD
    arma::uvec *alpha2k = new arma::uvec(h_size);
    arma::uvec *alpha2k1 = new arma::uvec(n_rows - h_size);
    rec_SVD(part1, h_size, k, alpha2k);
    rec_SVD(part2, n_rows - h_size, k, alpha2k1);

    std::chrono::duration<double, std::milli> interval_cpyres;
    t1 = std::chrono::high_resolution_clock::now();
    // cpy res (alpha2k id's == alphak SVD rank)
    for (arma::u64 row = 0; row < h_size; ++row)
    {
      (*alphak)[row] = rank[(*alpha2k)[row]];
    }
    for (arma::u64 row = h_size; row < n_rows; ++row)
    {
      (*alphak)[row] = rank[(*alpha2k1)[row - h_size] + h_size];
    }
    t2 = std::chrono::high_resolution_clock::now();
    interval_cpyres = t2 - t1;
    /* std::cout << n_rows <<" " << interval_SVD.count() / 1000 << " ";
    std::cout << interval_tri.count() / 1000 << " ";
    std::cout << interval_split.count() / 1000 << " ";
    std::cout << interval_cpyres.count() / 1000 << "\n"; */

    // Release memory
    delete alpha2k;
    delete alpha2k1;
    delete part1;
    delete part2;
  }
  else if (n_rows == 3)// part 2 has one element, he call SVD on part1 only.
  { 
    std::chrono::duration<double, std::milli> interval_tri;
    t1 = std::chrono::high_resolution_clock::now();
    // sort value : alphak[rank SVD] = ID
    arma::uvec rank = arma::sort_index(U.col(1).rows(0, n_rows - 1), "descend").as_col();
    t2 = std::chrono::high_resolution_clock::now();
    interval_tri = t2 - t1;

    std::chrono::duration<double, std::milli> interval_split;
    t1 = std::chrono::high_resolution_clock::now();
    // Part1 
    arma::u64 h_size = adj_Matrix->n_rows / 2;
    arma::sp_mat *part1 = new arma::sp_mat(h_size, adj_Matrix->n_cols);

    // Compute row order : order[id] = SVD rank
    arma::uvec order(n_rows);
    for (arma::u64 row = 0; row < n_rows; ++row)
    {
      order(rank(row)) = row;
    }

    // Fill part following row order
    arma::sp_mat::const_iterator it_end = adj_Matrix->end();
    for (arma::sp_mat::const_iterator it = adj_Matrix->begin(); it != it_end; ++it)
    {
      if (order[it.row()] < h_size)
      {
        (*part1)(order[it.row()], it.col()) = (*it);
      }
    }
    t2 = std::chrono::high_resolution_clock::now();
    interval_split = t2 - t1;

    // rec_SVD
    arma::uvec *alpha2k = new arma::uvec(h_size);
    rec_SVD(part1, h_size, k, alpha2k);

    std::chrono::duration<double, std::milli> interval_cpyres;
    t1 = std::chrono::high_resolution_clock::now();
    // cpy res
    for (arma::u64 row = 0; row < h_size; ++row)
    {
      (*alphak)[row] = rank[(*alpha2k)[row]];
    }
    for (arma::u64 row = h_size; row < n_rows; ++row)
    {
      (*alphak)[row] = rank[row];
    }
    t2 = std::chrono::high_resolution_clock::now();
    interval_cpyres = t2 - t1;
    /* std::cout << n_rows <<" " << interval_SVD.count() / 1000 << " ";
    std::cout << interval_tri.count() / 1000 << " ";
    std::cout << interval_split.count() / 1000 << " ";
    std::cout << interval_cpyres.count() / 1000 << "\n"; */

    // Release Memory
    delete alpha2k;
    delete part1;
  }
}

int mapping_SVD(graph *g, map_node *m)
{
  // Init temporary ID
  std::map<node, int> temp_id_map;
  std::map<node, int>::iterator item;
  int temp_id = 0; // nb of id

  // Give temporary id to each node in graph
  //  A changer en supposant que V = {0,id_max_fichier} (avec des nodes qui peuvent ne pas apparaître)
  std::chrono::duration<double, std::milli> interval_matrix;
  auto t1 = std::chrono::high_resolution_clock::now();
  for (graph::const_iterator e = g->begin(); e != g->end(); ++e)
  {
    item = temp_id_map.find(e->first);
    if (item == temp_id_map.end())
    {
      temp_id_map.insert({e->first, temp_id});
      temp_id++;
    }

    item = temp_id_map.find(e->second);
    if (item == temp_id_map.end())
    {
      temp_id_map.insert({e->second, temp_id});
      temp_id++;
    }
  }

  // Init adjMatrix
  size_t size = is_power_of_two(temp_id) ? temp_id : next_power_of_2(temp_id);
  arma::sp_mat *adjMatrix = new arma::sp_mat(size, size);

  // Fill adjMatrix
  for (graph::const_iterator e = g->begin(); e != g->end(); ++e)
  {
    (*adjMatrix)(temp_id_map[e->first], temp_id_map[e->second]) = 1;
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  interval_matrix = t2 - t1;
  std::cout << "Temps matrix :" << interval_matrix.count() / 1000 << "\n";
  // std::cout << "Non zeros :" << adjMatrix->n_nonzero << "\n"; // 0,133466721 % de 0
  // print_sp_mat(adjMatrix);

  // rec_SVD on adjMatrix
  std::chrono::duration<double, std::milli> interval_svd;
  t1 = std::chrono::high_resolution_clock::now();
  arma::uvec *alphak = new arma::uvec(temp_id); // alphak[SVD rank] == ID
  uint k = 2;
  rec_SVD(adjMatrix, temp_id, k, alphak);
  t2 = std::chrono::high_resolution_clock::now();
  interval_svd = t2 - t1;
  std::cout << "Temps rec SVD :" << interval_svd.count() / 1000 << std::endl;

  // cpy new ID (SVD rank) for each node in the graph
  for (map_node::const_iterator it = temp_id_map.begin(); it != temp_id_map.end(); ++it)
  {
    m->insert({it->first, (*alphak)[it->second]});
  }
  delete alphak;
  delete adjMatrix;
  return EXIT_SUCCESS;
}

int fct_map_edge(int i, int j)
{
  int max = (i > j) ? i : j;
  if (max == 1)
    return 1;

  // Previous largest power of two
  int k = is_power_of_two(max) ? max / 2 : next_power_of_2(max) / 2;

  // Apply the function
  int value = 0;
  if (i <= k && j > k)
  {
    value += pow(k, 2) + fct_map_edge(i, j - k);
  }
  else if (i > k && j <= k)
  {
    value += 2 * pow(k, 2) + fct_map_edge(i - k, j);
  }
  else if (i > k && j > k)
  {
    value += 3 * pow(k, 2) + fct_map_edge(i - k, j - k);
  }

  return value;
}