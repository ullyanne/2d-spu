/*
 * Decoder.cpp
 *
 *  Created on: Jan 14, 2011
 *      Author: rtoso
 */

#include "Decoder.h"

#include <math.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <set>

#include "MTRand.h"
#include "Packing.h"
using namespace std;

Decoder::~Decoder() {}

unsigned Decoder::local_search(std::vector<ranking>& solution,
                               unsigned current_cost, MTRand& rng,
                               int max_no_improve, double window_ls) const
{
  const int n = solution.size();
  int count = 0;

  int delta = n * window_ls;
  if (delta == 0) {
    delta = 1;
  }

  while (count < max_no_improve) {
    int i = rng.randInt(n - 1);
    int min_idx = max(0, i - delta);
    int max_idx = min(n - 1, i + delta);

    int j = min_idx + rng.randInt(max_idx - min_idx);

    while (j == i) {
      j = min_idx + rng.randInt(max_idx - min_idx);
    }

    swap(solution[i], solution[j]);

    double new_cost = pack(solution, items, max_width, ub, false);

    if (new_cost <= current_cost) {
      current_cost = new_cost;
    }
    else {
      swap(solution[i], solution[j]);
    }

    count++;
  }

  return current_cost;
}

double Decoder::decode(std::vector<double>& chromosome) const
{
  unsigned strip_height_plus_penalty = 0;
  std::vector<ranking> rank(chromosome.size());

  for (unsigned i = 0; i < chromosome.size(); i++) {
    rank[i].chromosome = chromosome[i];
    rank[i].index = i;
  }

  std::sort(rank.begin(), rank.end(), sort_rank);

  strip_height_plus_penalty = pack(rank, items, max_width, ub, false);

  return strip_height_plus_penalty;
}