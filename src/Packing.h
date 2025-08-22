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
  int index;
  int client;
} ranking;

inline bool operator==(const ranking &a, const ranking &b)
{
  return a.index == b.index;
}

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

int pack_with_one_layer(const std::vector<ranking> &rank,
                        const std::vector<item> &items, const int &max_width,
                        const int &ub, bool debug_sol = false,
                        std::fstream *solfile = nullptr);

unsigned pack(const std::vector<ranking> &rank, const std::vector<item> &items,
              const unsigned &max_width, const unsigned &ub,
              bool debug_sol = false, std::fstream *solfile = nullptr);

void encode(std::vector<ranking> &rank, const std::vector<ranking> &seq,
            unsigned n);

bool sort_rank(const ranking &a, const ranking &b);

void construct_vl_sol(std::vector<ranking> &sol, std::vector<double> chromosome,
                      std::vector<item> items);

void construct_final_sol(std::vector<ranking> &sol,
                         std::vector<double> chromosome,
                         std::vector<ranking> seq);

#endif  // PACKING_H