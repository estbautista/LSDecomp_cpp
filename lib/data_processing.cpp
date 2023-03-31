/* data_processing.cpp
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
#include "data_processing.hpp"
#include <chrono> // For test purpose

float_lien *float_lien_alloc()
{
    float_lien *FL = new float_lien();
    if (FL == NULL)
    {
        GSL_ERROR_VAL("failed to allocate space for float_lien class", GSL_ENOMEM, 0);
    }
    return FL;
}

void float_lien_delete(float_lien *FL)
{
    if (FL != NULL)
        delete FL;
}

int read_link_stream(std::string edgelist_path, float_lien *float_lien_to_fill, char delimiter)
{
    std::vector<std::string> lines;
    std::string line;

    if (float_lien_to_fill == NULL)
    {
        GSL_ERROR_VAL("Error float_lien_to_fill need to be allocated", GSL_EFAULT, 0);
    }

    // Open fill
    std::ifstream input_file(edgelist_path);
    if (!input_file.is_open())
    {
        std::cerr << "Could not open the file - '" << edgelist_path << "'\n";
        exit(EXIT_FAILURE);
    }

    // Read file
    while (getline(input_file, line))
    {
        lines.push_back(line);
    }
    if (lines.size() == 0)
    {
        std::cout << "Empty file - '" << edgelist_path << "'" << std::endl;
        input_file.close();
        exit(EXIT_FAILURE);
    }

    // Check weigthed graph format or not
    int count = 0;
    for (size_t i = 0; (i = lines[0].find(delimiter, i)) != std::string::npos; i++)
        count++;

    if (count == 2)
    { // no weigth col
        for (const auto &l : lines)
        {
            int first, second, t;// nodeA, nodeB, timeX
            std::sscanf(l.c_str(), "%d,%d,%d\n", &first, &second, &t);

            auto item = float_lien_to_fill->time_graph.find(t);
            if (item != float_lien_to_fill->time_graph.end())
            { // key time already exist in time_edgelist
                float_lien_to_fill->time_graph[t].push_back({{first, second}, 1.});
            }
            else
            {
                float_lien_to_fill->time_graph.insert({t, {{{first, second}, 1.}}});
            }
        }
    }
    else if (count == 3)
    { // weigth col
        for (const auto &l : lines)
        {
            int first, second, t;// nodeA, nodeB, timeX
            weigth w;// weigth
            std::sscanf(l.c_str(), "%d,%d,%d,%lf\n", &first, &second, &t, &w);

            auto item = float_lien_to_fill->time_graph.find(t);
            if (item != float_lien_to_fill->time_graph.end())
            { // key time already exist in time_edgelist
                float_lien_to_fill->time_graph[t].push_back({{first, second}, w});
            }
            else
            {
                float_lien_to_fill->time_graph.insert({t, {{{first, second}, w}}});
            }
        }
    }
    else
    {
        std::cout << "File format not supported- '" << edgelist_path << "'" << std::endl;
        input_file.close();
        exit(EXIT_FAILURE);
    }
    

    input_file.close();
    return (EXIT_SUCCESS);
}

int weigthedGraph_to_Graph(weigthed_graph *wg, graph *g)
{
    if (wg == NULL)
    {
        GSL_ERROR_VAL("weigthed_graph wg is NULL", GSL_EFAULT, 0);
    }
    if (g == NULL)
    {
        GSL_ERROR_VAL("Memory has to be allocated for graph g", GSL_EFAULT, 0);
    }

    for (weigthed_graph::const_iterator e_w = wg->begin(); e_w != wg->end(); ++e_w)
    {
        g->push_back(e_w->first);
    }

    return EXIT_SUCCESS;
}


// Print function for debugging purpose
void print_float_lien(float_lien *FL)
{
    std::cout << "===== Print time_edgelist =====\n";
    for (time_edgelist::const_iterator it = FL->time_graph.begin(); it != FL->time_graph.end(); ++it)
    {
        std::cout << it->first << " :";
        for (size_t j = 0; j < it->second.size(); j++)
        {
            std::cout << "\t(" << it->second[j].first.first << ";" << it->second[j].first.second << ")=" << it->second[j].second;
        }
    }
    std::cout << "\n===== End Print =====\n";
}
