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
  std::vector<ranking> rank(chromosome.size());

  for (unsigned i = 0; i < chromosome.size(); ++i) {
    rank[i].chromosome = chromosome[i];
    rank[i].index = subchromosome[i];
  }

  std::sort(rank.begin(), rank.end(), sort_rank);

  std::vector<ranking> lns_seq_copy = lns_seq;

  cout << "paçoca\n";
  for (unsigned i = 0; i < lns_seq_copy.size(); i++) {
    cout << lns_seq_copy[i].index << " ";
  }

  cout << "coco\n";
  std::vector<unsigned> subchromosome_copy = subchromosome;
  for (unsigned i = 0; i < chromosome.size(); i++) {
    subchromosome_copy[i] = rank[i].index;
    cout << subchromosome_copy[i] << " ";
  }

  rearrangeSeq(lns_seq_copy, subchromosome_copy);

  cout << "pudim\n";
  for (unsigned i = 0; i < lns_seq_copy.size(); i++) {
    unsigned a = lns_seq_copy[i].index;
    cout << lns_seq_copy[i].index << " ";
  }
  unsigned strip_height_plus_penalty = 0;

  unsigned num_layers = 0;
  return 0;

  if (pack_compressed) {
    std::cout << "firstye! ";
    for (unsigned i = 0; i < subchromosome_copy.size(); i++) {
      std::cout << subchromosome_copy[i] << " ";
    }

    std::cout << "secondye! ";
    for (unsigned i = 0; i < lns_seq_copy.size(); i++) {
      std::cout << lns_seq_copy[i].index << " ";
    }

    std::vector<std::vector<unsigned>> virtual_layers(5);
    unsigned best_height = 0;
    strip_height_plus_penalty = pack_with_one_layer(
        lns_seq_copy, items, max_width, ub, virtual_layers, best_height);
  }
  else {
    strip_height_plus_penalty =
        pack(lns_seq_copy, items, max_width, ub, clients_to_layers,
             layers_to_index, num_layers);
  }

  return strip_height_plus_penalty;
}