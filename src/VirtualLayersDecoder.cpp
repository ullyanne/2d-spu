/*
 * VirtualLayersDecoder.cpp
 *
 *  Created on: Jan 14, 2011
 *      Author: rtoso
 */

#include "VirtualLayersDecoder.h"

#include <math.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <set>

#include "MTRand.h"
#include "Packing.h"
using namespace std;

VirtualLayersDecoder::~VirtualLayersDecoder() {}

// Runs in \Theta(n \log n):
double VirtualLayersDecoder::decode(const std::vector<double>& chromosome) const
{
  unsigned strip_height_plus_penalty = 0;
  std::vector<ranking> rank(chromosome.size());

  unsigned idx = 0;
  for (unsigned i = 0; i < chromosome.size(); i++, idx++) {
    rank[idx].chromosome = chromosome[idx];
    rank[idx].index = i;
  }

  std::sort(rank.begin(), rank.end(), sort_rank);


  strip_height_plus_penalty =
      pack_with_one_layer(rank, items, max_width, ub, false);

  // std::vector<ranking> best_solution = seq;

  // unsigned pieces = 15;
  // std::vector<ranking> ls_seq =
  //     sliceLayers(rank, virtual_layers[current_layer - 1], pieces);

  // std::vector<vector<ranking>> remnd(2);
  // if (virtual_layers[current_layer - 1].size() > pieces) {
  //   remnd[0].insert(remnd[0].end(), virtual_layers[current_layer -
  //   1].begin(),
  //                   virtual_layers[current_layer - 1].end() - pieces);
  // }

  // if (rank.size() > pieces) {
  //   remnd[1].insert(remnd[1].end(), rank.begin() + pieces, rank.end());
  // }

  // for (unsigned i = 0; i < remnd[0].size(); i++) {
  //   std::cout << remnd[0][i].index << " " << std::flush;
  // }

  // std::vector<ranking> base;

  // if (current_layer > 0) {
  //   for (int j = 0; j < current_layer - 1; j++) {
  //     base.insert(base.end(), virtual_layers[j].begin(),
  //                 virtual_layers[j].end());
  //   }
  // }

  // double T = 1000;
  // double current_fitness = best_fitness;
  // double best_fit = current_fitness;
  // double max_iterations = 500;
  // double T_min = 0.1;
  // double alpha = 0.95;
  // MTRand ran;

  // for (unsigned iter = 0; iter < max_iterations && T > T_min; ++iter) {
  //   std::vector<ranking> neighbor = ls_seq;
  //   unsigned idx = ran.randInt(neighbor.size() - 1);
  //   unsigned idx2 = ran.randInt(neighbor.size() - 1);
  //   move_element(neighbor, idx, idx2);

  //   std::vector<ranking> seq2;
  //   seq2.insert(seq2.end(), base.begin(), base.end());
  //   seq2.insert(seq2.end(), remnd[0].begin(), remnd[0].end());
  //   seq2.insert(seq2.end(), neighbor.begin(), neighbor.end());
  //   seq2.insert(seq2.end(), remnd[1].begin(), remnd[1].end());

  //   unsigned new_fitness =
  //       pack_with_one_layer(seq2, items, max_width, ub, copy_change_later,
  //                           num_pieces_per_layer, false, best_height);

  //   double delta = new_fitness - current_fitness;
  //   if (delta < 0 || exp(-delta / T) > static_cast<double>(ran()) / RAND_MAX)
  //   {
  //     seq = neighbor;
  //     current_fitness = new_fitness;

  //     if (new_fitness < best_fit) {
  //       best_fitness = new_fitness;
  //       best_solution = neighbor;
  //     }
  //   }
  // }

  // if (best_fitness < 200) {
  //   cout << best_fitness << "\n";
  // }

  return strip_height_plus_penalty;
}