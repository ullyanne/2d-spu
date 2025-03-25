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

#include "Packing.h"
using namespace std;

VirtualLayersDecoder::~VirtualLayersDecoder() {}

// Runs in \Theta(n \log n):
double VirtualLayersDecoder::decode(const std::vector<double>& chromosome) const
{
  unsigned strip_height_plus_penalty = 0;
  std::vector<ranking> rank(chromosome.size());

  unsigned idx = 0;
  for (unsigned i = 0; i < virtual_layers[current_layer].size(); i++, idx++) {
    rank[idx].chromosome = chromosome[idx];
    rank[idx].index = virtual_layers[current_layer][i].index;
  }

  // for (unsigned i = 0; i < virtual_layers[current_layer + 1].size();
  //      i++, idx++) {
  //   rank[idx].chromosome = chromosome[idx];
  //   rank[idx].index = virtual_layers[current_layer + 1][i].index;
  // }

  std::sort(rank.begin(), rank.end(), sort_rank);

  std::vector<ranking> seq;
  seq.reserve(items.size());

  for (unsigned i = 0; i < virtual_layers.size(); i++) {
    if (i == current_layer) {
      seq.insert(seq.end(), rank.begin(), rank.end());
    }
    // else if (i == current_layer + 1) {
    //   continue;
    // }
    else {
      seq.insert(seq.end(), virtual_layers[i].begin(), virtual_layers[i].end());
    }
  }

  unsigned best_height = 0;
  std::vector<std::vector<ranking>> copy_change_later;
  strip_height_plus_penalty =
      pack_with_one_layer(seq, items, max_width, ub, copy_change_later,
                          num_pieces_per_layer, false, best_height);

  return strip_height_plus_penalty;
}