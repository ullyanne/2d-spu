/*
 * SampleDecoder.cpp
 *
 *  Created on: Jan 14, 2011
 *      Author: rtoso
 */

#include "SampleDecoder.h"

#include <math.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <set>

#include "Packing.h"
using namespace std;

SampleDecoder::~SampleDecoder() {}

// Runs in \Theta(n \log n):
double SampleDecoder::decode(const std::vector<double> &chromosome) const
{
  std::vector<ranking> rank(chromosome.size());

  for (unsigned i = 0; i < chromosome.size(); ++i) {
    rank[i].chromosome = chromosome[i];
    rank[i].index = seq[i].index;
  }

  std::sort(rank.begin(), rank.end(), sort_rank);

  unsigned strip_height_plus_penalty =
      pack_with_one_layer(rank, items, max_width, ub, false);

  return strip_height_plus_penalty;
}