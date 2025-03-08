#include "Packing.h"

#include <algorithm>
#include <boost/container/flat_set.hpp>
#include <deque>
#include <fstream>

using namespace std;
using boost::container::flat_set;

void open_new_layer(vector<flat_set<ems_t, bottom_left_cmp>> &layers,
                    unsigned &num_layers, unsigned &strip_height, int max_width,
                    int item_height)
{
  layers.push_back(flat_set<ems_t, bottom_left_cmp>());
  num_layers++;

  ems_t space;
  space.bottom_point = make_pair(0, strip_height);
  space.top_point = make_pair(max_width, item_height + strip_height);

  layers[num_layers - 1].insert(space);
  strip_height += item_height;
}

bool item_can_fit(item item, ems_t ems_t)
{
  int space_width = ems_t.top_point.first - ems_t.bottom_point.first;
  int space_height = ems_t.top_point.second - ems_t.bottom_point.second;

  if (item.width <= space_width && item.height <= space_height) {
    return true;
  }
  return false;
}

bool is_a_line(ems_t space)
{
  if (space.bottom_point.first == space.top_point.first or
      space.bottom_point.second == space.top_point.second) {
    return true;
  }
  return false;
}

bool is_contained(ems_t a, ems_t b)
{
  if (a.top_point.first <= b.top_point.first &&
      a.top_point.second <= b.top_point.second &&
      a.bottom_point.first >= b.bottom_point.first &&
      a.bottom_point.second >= b.bottom_point.second) {
    return true;
  }
  return false;
}

bool is_maximal(ems_t new_space, flat_set<ems_t, bottom_left_cmp> &layer)
{
  for (const auto &space : layer) {
    if (is_contained(new_space, space) || is_contained(space, new_space)) {
      return false;
    }
  }
  return true;
}

bool is_ems_valid(ems_t space, flat_set<ems_t, bottom_left_cmp> &layer)
{
  if (is_a_line(space)) {
    return false;
  }

  if (!is_maximal(space, layer)) {
    return false;
  }

  if (space.bottom_point.second > space.top_point.second) {
    return false;
  }

  if (space.bottom_point.first > space.top_point.first) {
    return false;
  }

  return true;
}

void calc_diff_process(flat_set<ems_t, bottom_left_cmp> &layer,
                       std::deque<ems_t> &layer_it, ems_t space, int x3, int y3,
                       int x4, int y4)
{
  int x1 = space.bottom_point.first;
  int y1 = space.bottom_point.second;
  int x2 = space.top_point.first;
  int y2 = space.top_point.second;

  vector<int> diff_proc = {x1, y1, x3, y2, x4, y1, x2, y2,
                           x1, y1, x2, y3, x1, y4, x2, y2};

  for (int i = 0; i <= 12; i += 4) {
    space.bottom_point = make_pair(diff_proc[i], diff_proc[i + 1]);
    space.top_point = make_pair(diff_proc[i + 2], diff_proc[i + 3]);

    if (is_ems_valid(space, layer)) {
      bool intersects_x =
          !(space.top_point.first <= x3 || space.bottom_point.first >= x4);
      bool below_piece = space.bottom_point.second < y3;
      if (intersects_x && below_piece) {
        continue;
      }

      layer.insert(space);
      layer_it.push_back(space);
    }
  }
}

bool does_intersect(int x3, int y3, int x4, int y4, ems_t space)
{
  return (x3 < space.top_point.first && x4 > space.bottom_point.first &&
          y3 < space.top_point.second && y4 > space.bottom_point.second);
}

void calculate_item_coordinates(int &x3, int &y3, int &x4, int &y4, item item,
                                ems_t space)
{
  x3 = space.bottom_point.first;
  y3 = space.bottom_point.second;
  x4 = x3 + item.width;
  y4 = y3 + item.height;
}

bool will_violate_unloading_constraint(item item, ems_t space,
                                       int clients_verification[])
{
  int x3, y3, x4, y4;
  calculate_item_coordinates(x3, y3, x4, y4, item, space);

  for (int i = x3; i < x4; i++) {
    if (clients_verification[i] != -1 &&
        clients_verification[i] < item.client) {
      return true;
    }
  }
  return false;
}

void fit_item(item item, ems_t space, flat_set<ems_t, bottom_left_cmp> &layer,
              int clients_verification[], unsigned ub, unsigned *penalty,
              bool debug_sol = false, fstream *solfile = nullptr)
{
  int x3, y3, x4, y4;
  calculate_item_coordinates(x3, y3, x4, y4, item, space);

  for (int i = x3; i < x4; i++) {
    if (clients_verification[i] != -1 &&
        clients_verification[i] < item.client) {
      *penalty += item.client - clients_verification[i];
      break;
    }
  }

  for (int i = x3; i < x4; i++) {
    clients_verification[i] = item.client;
  }

  if (debug_sol) {
    *solfile << x3 << " " << y3 << "\n" << x4 << " " << y4 << "\n";
  }

  layer.erase(space);

  std::deque<ems_t> to_process(layer.begin(), layer.end());

  calc_diff_process(layer, to_process, space, x3, y3, x4, y4);

  while (!to_process.empty()) {
    ems_t ems = to_process.front();
    to_process.pop_front();

    if (does_intersect(x3, y3, x4, y4, ems)) {
      ems_t ems_tmp = ems;
      layer.erase(ems);
      calc_diff_process(layer, to_process, ems_tmp, x3, y3, x4, y4);
    }
    else {
      bool intersects_x =
          !(space.top_point.first <= x3 || space.bottom_point.first >= x4);
      bool below_piece = space.bottom_point.second < y3;
      if (intersects_x && below_piece) {
        layer.erase(ems);
      }
    }
  }
}

void rearrangeSeq(std::vector<ranking> &arr,
                  const std::vector<unsigned> &indices_to_move)
{
  std::vector<ranking> result;
  std::vector<bool> copied(arr.size(), 0);

  for (unsigned j = 0; j < indices_to_move.size(); j++) {
    if (!copied[indices_to_move[j]]) {
      result.push_back(arr[indices_to_move[j]]);
      copied[indices_to_move[j]] = 1;
    }
  }

  for (unsigned i = 0; i < arr.size(); i++) {
    if (!copied[i]) {
      result.push_back(arr[i]);
      copied[i] = 1;
    }
  }

  arr = result;
}

bool sort_rank(const ranking &a, const ranking &b)
{
  return a.chromosome < b.chromosome;
}

void construct_sol(std::vector<ranking> &full_solution,
                   std::vector<double> chromosome,
                   std::vector<unsigned> subchromosome,
                   std::vector<ranking> lns_seq)
{
  std::vector<ranking> rank(chromosome.size());

  for (unsigned i = 0; i < chromosome.size(); ++i) {
    rank[i].chromosome = chromosome[i];
    rank[i].index = subchromosome[i];
  }
  std::sort(rank.begin(), rank.end(), sort_rank);

  full_solution = lns_seq;

  std::vector<unsigned> subchromosome_copy = subchromosome;
  for (unsigned i = 0; i < chromosome.size(); i++) {
    subchromosome_copy[i] = rank[i].index;
  }

  rearrangeSeq(full_solution, subchromosome_copy);
}

unsigned pack(
    const std::vector<ranking> &rank, vector<item> items,
    const unsigned max_width, const unsigned ub,
    std::unordered_map<unsigned, std::vector<unsigned>> &clients_to_layers,
    std::unordered_map<unsigned, unsigned> &layers_to_index,
    unsigned &num_layers, bool debug_sol, std::fstream *solfile)
{
  // vector<flat_set<ems_t, bottom_left_cmp>> layers;
  vector<flat_set<ems_t, bottom_left_cmp>> layers;

  unsigned item_index = rank[0].index;
  item item = items[item_index];
  unsigned strip_height = 0;
  num_layers = 0;
  long unsigned int items_placed = 0;
  bool fit;
  unsigned penalty = 0;

  std::vector<unsigned> layers_client;
  unsigned current_client = item.client;

  int clients_verification[max_width];

  for (unsigned i = 0; i < max_width; i++) {
    clients_verification[i] = -1;
  }

  open_new_layer(layers, num_layers, strip_height, max_width, item.height);
  layers_to_index[num_layers - 1] = 0;
  layers_client.push_back(num_layers - 1);

  while (items_placed < items.size()) {
    item_index = rank[items_placed].index;
    item = items[item_index];
    fit = false;

    if (!layers[num_layers - 1].empty()) {
      for (auto ems_t = layers[num_layers - 1].begin();
           ems_t != layers[num_layers - 1].end();) {
        if (item_can_fit(item, *ems_t)) {
          fit_item(item, *ems_t, layers[num_layers - 1], clients_verification,
                   ub, &penalty, debug_sol, solfile);
          fit = true;
          items_placed++;

          if (item.client != current_client) {
            clients_to_layers[current_client] = layers_client;

            layers_client.clear();

            current_client = item.client;
            layers_client.push_back(num_layers - 1);
          }
          break;
        }
        else {
          ems_t++;
        }
      }
    }

    if (!fit) {
      open_new_layer(layers, num_layers, strip_height, max_width, item.height);
      if (item.client != current_client) {
        clients_to_layers[current_client] = layers_client;

        layers_client.clear();

        current_client = item.client;
        layers_client.push_back(num_layers - 1);
      }
      else {
        layers_client.push_back(num_layers - 1);
      }
      auto new_space = *layers[num_layers - 1].begin();
      fit_item(item, new_space, layers[num_layers - 1], clients_verification,
               ub, &penalty, debug_sol, solfile);

      fit = true;
      items_placed++;
      layers_to_index[num_layers - 1] = items_placed - 1;
    }

    if (debug_sol) {
      *solfile << item_index << "\n" << item.client << "\n";
    }
  }

  clients_to_layers[current_client] = layers_client;

  if (penalty) {
    penalty += ub;
  }

  return strip_height + penalty;
}

void calc_strip_height(unsigned &strip_height, ems_t space,
                       unsigned item_height)
{
  unsigned height = space.bottom_point.second + item_height;
  if (height > strip_height) {
    strip_height = height;
  }
}

unsigned pack_with_one_layer(const std::vector<ranking> &rank,
                             vector<item> items, const unsigned max_width,
                             const unsigned ub, bool debug_sol,
                             std::fstream *solfile)
{
  vector<flat_set<ems_t, bottom_left_cmp>> layers;

  unsigned item_index = rank[0].index;
  item item = items[item_index];
  unsigned strip_height = 0;
  long unsigned int items_placed = 0;
  bool fit;
  unsigned penalty = 0;

  unsigned current_client = item.client;

  int clients_verification[max_width];

  for (unsigned i = 0; i < max_width; i++) {
    clients_verification[i] = -1;
  }

  layers.push_back(flat_set<ems_t, bottom_left_cmp>());
  ems_t space;
  space.bottom_point = make_pair(0, 0);
  space.top_point = make_pair(max_width, ub);
  layers[0].insert(space);

  while (items_placed < items.size()) {
    item_index = rank[items_placed].index;
    item = items[item_index];
    fit = false;

    if (!layers[0].empty()) {
      for (auto ems_t = layers[0].begin(); ems_t != layers[0].end();) {
        if (item_can_fit(item, *ems_t)) {
          calc_strip_height(strip_height, *ems_t, item.height);
          fit_item(item, *ems_t, layers[0], clients_verification, ub, &penalty,
                   debug_sol, solfile);
          if (debug_sol) {
            *solfile << item_index << "\n" << item.client << "\n";
          }
          fit = true;
          items_placed++;
          break;
        }
        else {
          ems_t++;
        }
      }
    }
  }

  if (penalty) {
    penalty += ub;
  }

  return strip_height + penalty;
}

unsigned pack_compressed(const std::vector<ranking> &rank, vector<item> items,
                         const unsigned max_width, const unsigned ub,
                         bool debug_sol, std::fstream *solfile)
{
  // vector<flat_set<ems_t, bottom_left_cmp>> layers;
  vector<flat_set<ems_t, bottom_left_cmp>> layers;

  unsigned item_index = rank[0].index;
  item item = items[item_index];
  unsigned strip_height = 0;
  long unsigned int items_placed = 0;
  bool fit;
  unsigned penalty;

  unsigned current_client = item.client;

  int clients_verification[max_width];

  for (unsigned i = 0; i < max_width; i++) {
    clients_verification[i] = -1;
  }

  layers.push_back(flat_set<ems_t, bottom_left_cmp>());
  ems_t space;
  space.bottom_point = make_pair(0, 0);
  space.top_point = make_pair(max_width, ub);
  layers[0].insert(space);

  unsigned index = 0;
  std::vector<ranking> rank_copy = rank;

  while (items_placed < items.size()) {
    item_index = rank_copy[index].index;
    item = items[item_index];
    fit = false;

    if (!layers[0].empty()) {
      for (auto ems_t = layers[0].begin(); ems_t != layers[0].end();) {
        if (item_can_fit(item, *ems_t) &&
            !will_violate_unloading_constraint(item, *ems_t,
                                               clients_verification)) {
          calc_strip_height(strip_height, *ems_t, item.height);
          fit_item(item, *ems_t, layers[0], clients_verification, ub, &penalty,
                   debug_sol, solfile);
          if (debug_sol) {
            *solfile << item_index << "\n" << item.client << "\n";
          }
          fit = true;
          items_placed++;
          index = 0;
          rank_copy.erase(rank_copy.begin() + index);
          break;
        }
        else if (index == rank_copy.size() - 1) {
          ems_t = layers[0].erase(ems_t);
        }
        else {
          index++;
          // ems_t++;
        }
      }
    }
  }

  return strip_height;
}