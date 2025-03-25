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
double SampleDecoder::decode(
    const std::vector<double> &chromosome,
    std::unordered_map<unsigned, std::vector<unsigned>> &clients_to_layers,
    std::unordered_map<unsigned, unsigned> &layers_to_index) const
{
  unsigned strip_height_plus_penalty = 0;
  std::vector<ranking> rank(chromosome.size());

  for (unsigned i = 0; i < chromosome.size(); ++i) {
    rank[i].chromosome = chromosome[i];
    rank[i].index = subchromosome[i];
  }

  std::sort(rank.begin(), rank.end(), sort_rank);

  unsigned num_layers = 0;

  std::vector<ranking> lns_seq_copy = lns_seq;

  std::vector<unsigned> subchromosome_copy = subchromosome;
  for (unsigned i = 0; i < chromosome.size(); i++) {
    subchromosome_copy[i] = rank[i].index;
  }

  rearrangeSeq(lns_seq_copy, subchromosome_copy);

  strip_height_plus_penalty =
      pack(lns_seq_copy, items, max_width, ub, clients_to_layers,
           layers_to_index, num_layers);

  return strip_height_plus_penalty;
}