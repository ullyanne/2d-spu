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

  // for (unsigned i = 0; i < lns_seq.size(); i++) {
  //   std::cout << "seq " << lns_seq[i].index << "\n";
  // }

  std::vector<unsigned> subchromosome_copy = subchromosome;
  for (unsigned i = 0; i < chromosome.size(); i++) {
    // std::cout << lns_seq_copy[i].index << "\n";
    // std::cout << "bom dia!" << subchromosome[i] << "\n";
    // std::cout << "b!" << rank[i].index << "\n";
    // std::cout << "d!" << lns_seq[rank[i].index].index << "\n";
    subchromosome_copy[i] = rank[i].index;
    // lns_seq_copy[subchromosome[i]] = lns_seq[rank[i].index];
  }

  // lns_seq_copy[0] = lns_seq[8];
  // lns_seq_copy[1] = lns_seq[7];
  // lns_seq_copy[2] = lns_seq[1];
  // lns_seq_copy[3] = lns_seq[0];
  // lns_seq_copy[4] = lns_seq[5];
  // lns_seq_copy[5] = lns_seq[3];
  // lns_seq_copy[7] = lns_seq[4];
  // lns_seq_copy[8] = lns_seq[2];
  // lns_seq_copy[9] = lns_seq[9];

  rearrangeSeq(lns_seq_copy, subchromosome_copy);
  // for (unsigned i = 0; i < lns_seq_copy.size(); i++) {
  //   std::cout << "seq copy " << lns_seq_copy[i].index << "\n";
  // }

  // best sol!
  // lns_seq_copy[0] = lns_seq[7];
  // lns_seq_copy[1] = lns_seq[1];
  // lns_seq_copy[2] = lns_seq[4];
  // lns_seq_copy[3] = lns_seq[5];
  // lns_seq_copy[4] = lns_seq[3];
  // lns_seq_copy[5] = lns_seq[0];
  // lns_seq_copy[6] = lns_seq[2];
  // lns_seq_copy[7] = lns_seq[9];
  // lns_seq_copy[8] = lns_seq[9];
  // lns_seq_copy[9] = lns_seq[6];

  // for (int i = 0; i < lns_seq_copy.size(); i++) {
  //   // std::cout << "hell!\n";
  //   std::cout << lns_seq_copy[i].index << " ";
  // }

  unsigned strip_height_plus_penalty = 0;

  unsigned num_layers = 0;

  strip_height_plus_penalty =
      pack(lns_seq_copy, items, max_width, ub, clients_to_layers,
           layers_to_index, num_layers);

  // if (strip_height_plus_penalty == 40) {
  //   cout << "\nh: " << strip_height_plus_penalty << "\n";

  //   for (int i = 0; i < lns_seq_copy.size(); i++) {
  //     cout << "baba " << lns_seq_copy[i].index << " ";
  //   }
  // }

  return strip_height_plus_penalty;
}