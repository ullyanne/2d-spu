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
#include <set>
using namespace std;

SampleDecoder::~SampleDecoder() {}
typedef struct {
  double chromosome;
  unsigned index;
  unsigned client;
} ranking;

typedef struct {
  pair<int, int> bottom_point;
  pair<int, int> top_point;
} ems;

struct bottom_left_cmp {
  bool operator()(ems a, ems b) const
  {
    if (a.bottom_point.second == b.bottom_point.second) {
      return a.bottom_point.first < b.bottom_point.first;
    }
    return a.bottom_point.second < b.bottom_point.second;
  }
};

bool sort_rank(const ranking &a, const ranking &b)
{
  if (a.client != b.client) {
    return a.client > b.client;
  }
  return a.chromosome < b.chromosome;
}

bool item_can_fit(item item, ems ems)
{
  int space_width = ems.top_point.first - ems.bottom_point.first;
  int space_height = ems.top_point.second - ems.bottom_point.second;

  if (item.width <= space_width && item.height <= space_height) {
    return true;
  }
  return false;
}

bool is_a_line(ems space)
{
  if (space.bottom_point.first == space.top_point.first or
      space.bottom_point.second == space.top_point.second) {
    return true;
  }
  return false;
}

bool is_contained(ems a, ems b)
{
  if (a.top_point.first <= b.top_point.first &&
      a.top_point.second <= b.top_point.second &&
      a.bottom_point.first >= b.bottom_point.first &&
      a.bottom_point.second >= b.bottom_point.second) {
    return true;
  }
  return false;
}

bool is_maximal(ems new_space, set<ems, bottom_left_cmp> &layer)
{
  for (const auto &space : layer) {
    if (is_contained(new_space, space)) {
      return false;
    }
    if (is_contained(space, new_space)) {
      return false;
    }
  }
  return true;
}

bool is_ems_valid(ems space, set<ems, bottom_left_cmp> &layer)
{
  if (is_a_line(space)) {
    return false;
  }

  if (!is_maximal(space, layer)) {
    return false;
  }

  return true;
}

void fit_item(item item, ems space,
              set<ems, bottom_left_cmp> &layer)
{
  int x1 = space.bottom_point.first;
  int y1 = space.bottom_point.second;
  int x2 = space.top_point.first;
  int y2 = space.top_point.second;
  int x4 = x1 + item.width;
  int y4 = y1 + item.height;

  vector<int> diff_proc = {x1, y1, x1, y2, x4, y1, x2, y2,
                           x1, y1, x2, y1, x1, y4, x2, y2};

  for (int i = 0; i <= 12; i += 4) {
    space.bottom_point = make_pair(diff_proc[i], diff_proc[i + 1]);
    space.top_point = make_pair(diff_proc[i + 2], diff_proc[i + 3]);

    if (is_ems_valid(space, layer)) {
      layer.insert(space);
    }
  }
}

void open_new_layer(vector<set<ems, bottom_left_cmp>> &layers,
                    int &num_layers, int &strip_height, int max_width,
                    int item_height)
{
  layers.push_back(set<ems, bottom_left_cmp>());
  num_layers++;

  ems space;
  space.bottom_point = make_pair(0, strip_height);
  space.top_point = make_pair(max_width, item_height);

  layers[num_layers - 1].insert(space);
  strip_height += item_height;
}

// Runs in \Theta(n \log n):
double SampleDecoder::decode(const std::vector<double> &chromosome) const
{
  std::vector<ranking> rank(chromosome.size());

  for (unsigned i = 0; i < chromosome.size(); ++i) {
    rank[i].chromosome = chromosome[i];
    rank[i].index = i;
    rank[i].client = items[i].client;
  }
  // rank[0].client = 2;
  // rank[0].index = 5;
  // rank[1].client = 2;
  // rank[1].index = 4;
  // rank[2].client = 2;
  // rank[2].index = 3;
  // rank[3].client = 2;
  // rank[3].index = 2;
  // rank[4].client = 2;
  // rank[4].index = 1;
  // rank[5].client = 2;
  // rank[5].index = 0;
  // rank[6].client = 1;
  // rank[6].index = 8;
  // rank[7].client = 1;
  // rank[7].index = 7;
  // rank[8].client = 1;
  // rank[8].index = 9;
  // rank[9].client = 1;
  // rank[9].index = 6;

  // // Here we sort 'permutation', which will then produce a permutation of [n]
  // in
  // // pair::second:
  std::sort(rank.begin(), rank.end(), sort_rank);

  // for (unsigned i = 0; i < chromosome.size(); ++i) {
  //   cout << rank[i].client << " " << rank[i].chromosome << " " <<
  //   rank[i].index
  //        << " " << endl;
  // }

  vector<set<ems, bottom_left_cmp>> layers;

  int item_index = rank[0].index;
  item item = items[item_index];
  int prev_client = item.client;
  int strip_height = 0;
  int current_layer = 0;
  int num_layers = 0;
  long unsigned int items_placed = 0;
  bool fit;

  open_new_layer(layers, num_layers, strip_height, max_width, item.height);
  // layers.push_back(set<ems, decltype(bottom_left_compare)>());
  // ems space;
  // space.bottom_point = make_pair(0, 0);
  // space.top_point = make_pair(max_width, item.height);
  // layers[0].insert(space);
  /* aqui ainda não coloquei a peça; só criei o espaço!*/

  while (items_placed < items.size()) {
    item_index = rank[items_placed].index;
    item = items[item_index];
    fit = false;

    if (item.client != prev_client) {
      current_layer = num_layers - 1;
    }

    for (int i = current_layer; i < num_layers && !fit; i++) {
      if (!layers[i].empty()) {
        for (auto ems = layers[i].begin(); ems != layers[i].end(); ems++) {
          if (item_can_fit(item, *ems)) {
            auto space = *ems;
            ems = layers[i].erase(ems);
            fit_item(item, space, layers[i]);
            fit = true;
            items_placed++;
            break;
          }
        }
      }

      if (!fit && i == num_layers - 1) {
        open_new_layer(layers, num_layers, strip_height, max_width,
                       item.height);
        auto new_space = *layers[num_layers - 1].begin();
        layers[num_layers - 1].erase(new_space);
        fit_item(item, new_space, layers[num_layers - 1]);
        fit = true;
        items_placed++;
      }
    }

    prev_client = item.client;
  }

  cout << strip_height << endl;
  return strip_height;

  // std::pair<int, int> piece;
  // std::vector<std::pair<int, int>> spaces;
  // unsigned index;

  // for(unsigned i = 0; i < chromosome.size(); i++){
  //     index = ranking[i].second;
  //     piece = pieces[index];

  // }

  // std::cout << pieces[0].first << " " << pieces[0].second << std::endl;
  /* Primeira peça a ser inserida*/
  // piece = pieces[ranking[pieceIndex].second];

  // return strip_height;
}