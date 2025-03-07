#ifndef PACKING_H
#define PACKING_H

#include <algorithm>
#include <iostream>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Item.h"

typedef struct {
  double chromosome;
  unsigned index;
  unsigned client;
} ranking;

typedef struct {
  std::pair<int, int> bottom_point;
  std::pair<int, int> top_point;
} ems_t;

struct bottom_left_cmp {
  bool operator()(ems_t a, ems_t b) const
  {
    if (a.bottom_point.second == b.bottom_point.second) {
      if (a.bottom_point.first != b.bottom_point.first) {
        return a.bottom_point.first < b.bottom_point.first;
      }
      return a.top_point.first > b.top_point.first;
    }
    return a.bottom_point.second < b.bottom_point.second;
  }
};

unsigned pack(
    const std::vector<ranking> &rank, std::vector<item> items,
    const unsigned max_width, const unsigned ub,
    std::unordered_map<unsigned, std::vector<unsigned>> &clients_to_layers,
    std::unordered_map<unsigned, unsigned> &layers_to_index,
    unsigned &num_layers, bool debug_sol = false,
    std::fstream *solfile = nullptr);

unsigned pack_with_one_layer(const std::vector<ranking> &rank,
                             std::vector<item> items, const unsigned max_width,
                             const unsigned ub,
                             std::vector<std::vector<unsigned>> &virtual_layers,
                             unsigned &best_height, bool debug_sol = false,
                             std::fstream *solfile = nullptr);

unsigned pack_compressed(const std::vector<ranking> &rank,
                         std::vector<item> items, const unsigned max_width,
                         const unsigned ub, bool debug_sol,
                         std::fstream *solfile);

bool sort_rank(const ranking &a, const ranking &b);

void construct_sol(std::vector<ranking> &full_solution,
                   std::vector<double> chromosome,
                   std::vector<unsigned> subchromosome,
                   std::vector<ranking> lns_seq);

void rearrangeSeq(std::vector<ranking> &arr,
                  const std::vector<unsigned> &indices_to_move);

#endif  // PACKING_H