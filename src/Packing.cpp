#include "Packing.h"

#include <algorithm>
#include <boost/container/flat_set.hpp>
#include <chrono>
#include <climits>
#include <cmath>
#include <deque>
#include <fstream>
#include <optional>

using namespace std;
using boost::container::flat_set;

struct FlatSegmentTree {
  int size;
  vector<int> min_val, lazy;
  vector<bool> pending;

  FlatSegmentTree(int n)
  {
    size = 1;
    while (size < n) size *= 2;
    min_val.assign(2 * size, INT_MAX);
    lazy.assign(2 * size, INT_MAX);
    pending.assign(2 * size, false);
  }

  void apply(int x, int l, int r, int val)
  {
    min_val[x] = val;
    if (r > l) {
      lazy[x] = val;
      pending[x] = true;
    }
  }

  void push(int x, int l, int r)
  {
    if (pending[x]) {
      int mid = (l + r) / 2;
      apply(2 * x, l, mid, lazy[x]);
      apply(2 * x + 1, mid + 1, r, lazy[x]);
      pending[x] = false;
    }
  }

  void update(int x, int l, int r, int ql, int qr, int val)
  {
    if (qr < l || r < ql) return;
    if (ql <= l && r <= qr) {
      apply(x, l, r, val);
      return;
    }
    push(x, l, r);
    int mid = (l + r) / 2;
    update(2 * x, l, mid, ql, qr, val);
    update(2 * x + 1, mid + 1, r, ql, qr, val);
    min_val[x] = min(min_val[2 * x], min_val[2 * x + 1]);
  }

  int query(int x, int l, int r, int ql, int qr)
  {
    if (qr < l || r < ql) return INT_MAX;
    if (ql <= l && r <= qr) return min_val[x];
    push(x, l, r);
    int mid = (l + r) / 2;
    return min(query(2 * x, l, mid, ql, qr),
               query(2 * x + 1, mid + 1, r, ql, qr));
  }

  // Interface externa
  void update(int l, int r, int val) { update(1, 0, size - 1, l, r, val); }

  int query(int l, int r) { return query(1, 0, size - 1, l, r); }
};

inline int width(const ems_t &e)
{
  return e.top_point.first - e.bottom_point.first;
}
inline int height(const ems_t &e)
{
  return e.top_point.second - e.bottom_point.second;
}

void open_new_layer(vector<flat_set<ems_t, bottom_left_cmp>> &layers,
                    unsigned &num_layers, unsigned &strip_height,
                    const int &max_width, const int &item_height)
{
  layers.emplace_back();
  ++num_layers;

  int new_strip_height = strip_height + item_height;

  // ems_t space;
  // space.bottom_point = make_pair(0, strip_height);
  // space.top_point = make_pair(max_width, item_height + strip_height);

  // layers[num_layers - 1].insert(space);

  layers.back().insert({{0, strip_height}, {max_width, new_strip_height}});

  strip_height = new_strip_height;
}

inline bool item_can_fit(const item &item, const ems_t &ems)
{
  // int space_width = ems_t.top_point.first - ems_t.bottom_point.first;
  // int space_height = ems_t.top_point.second - ems_t.bottom_point.second;

  // if (item.width <= space_width && item.height <= space_height) {
  //   return true;
  // }
  // return false;
  return item.width <= (width(ems)) && item.height <= (height(ems));
}

inline bool is_a_line(const ems_t &space)
{
  return (space.bottom_point.first == space.top_point.first ||
          space.bottom_point.second == space.top_point.second);
}

inline bool is_contained(const ems_t &a, const ems_t &b)
{
  return (a.top_point.first <= b.top_point.first &&
          a.top_point.second <= b.top_point.second &&
          a.bottom_point.first >= b.bottom_point.first &&
          a.bottom_point.second >= b.bottom_point.second);
}

bool is_maximal(const ems_t &new_space, flat_set<ems_t, bottom_left_cmp> &layer)
{
for (auto space = layer.begin(); space != layer.end(); ) {
    if (space->bottom_point.second > new_space.top_point.second) {
      break;
    }

    if (is_contained(new_space, *space)) {
      return false;
    }

    if (is_contained(*space, new_space)) {
      space = layer.erase(space); 
    } else {
      ++space;
    }
  }

  return true;

}

bool is_ems_valid(const ems_t &space, flat_set<ems_t, bottom_left_cmp> &layer)
{
  if (is_a_line(space)) {
    return false;
  }

  if (space.bottom_point.second > space.top_point.second) {
    return false;
  }

  if (space.bottom_point.first > space.top_point.first) {
    return false;
  }

  if (!is_maximal(space, layer)) {
    return false;
  }

  return true;
}

void calc_diff_process(flat_set<ems_t, bottom_left_cmp> &layer,
                       std::vector<ems_t> &layer_it, ems_t space, const int &x3,
                       const int &y3, const int &x4, const int &y4)
{
  int x1 = space.bottom_point.first;
  int y1 = space.bottom_point.second;
  int x2 = space.top_point.first;
  int y2 = space.top_point.second;

  // vector<int> diff_proc = {x1, y1, x3, y2, x4, y1, x2, y2,
  //                          x1, y1, x2, y3, x1, y4, x2, y2};
  for (int i = 0; i <= 12; i += 4) {
    switch (i) {
      case 0:
        space = {{x1, y1}, {x3, y2}};
        break;
      case 4:
        space = {{x4, y1}, {x2, y2}};
        break;
      case 8:
        space = {{x1, y1}, {x2, y3}};
        break;
      case 12:
        space = {{x1, y4}, {x2, y2}};
        break;
    }

    // space.bottom_point = make_pair(diff_proc[i], diff_proc[i + 1]);
    // space.top_point = make_pair(diff_proc[i + 2], diff_proc[i + 3]);

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

inline bool does_intersect(const int &x3, const int &y3, const int &x4,
                           const int &y4, const ems_t &space)
{
  return (x3 < space.top_point.first && x4 > space.bottom_point.first &&
          y3 < space.top_point.second && y4 > space.bottom_point.second);
}

void calculate_item_coordinates(int &x3, int &y3, int &x4, int &y4,
                                const item &item, const ems_t &space)
{
  x3 = space.bottom_point.first;
  y3 = space.bottom_point.second;
  x4 = x3 + item.width;
  y4 = y3 + item.height;
}

void fit_item(const item &item, ems_t space,
              flat_set<ems_t, bottom_left_cmp> &layer, const unsigned &ub,
              unsigned *penalty, FlatSegmentTree &seg_tree,
              bool debug_sol = false, fstream *solfile = nullptr)
{
  // auto start = std::chrono::high_resolution_clock::now();
  int x3, y3, x4, y4;
  calculate_item_coordinates(x3, y3, x4, y4, item, space);

  int min_client = seg_tree.query(x3, x4 - 1);
  if (min_client < item.client) {
    *penalty += item.client - min_client;
  }
  seg_tree.update(x3, x4 - 1, item.client);

  // for (int i = x3; i < x4; i++) {
  //   if (clients_verification[i] != -1 &&
  //       clients_verification[i] < item.client) {
  //     *penalty += item.client - clients_verification[i];
  //     break;
  //   }
  // }

  // std::fill(clients_verification + x3, clients_verification + x4,
  // item.client);

  if (debug_sol) {
    *solfile << x3 << " " << y3 << "\n" << x4 << " " << y4 << "\n";
  }

  layer.erase(space);

  // std::deque<ems_t> to_process(layer.begin(), layer.end());
  std::vector<ems_t> to_process(layer.begin(), layer.end());

  // auto end = std::chrono::high_resolution_clock::now();
  // auto duration =
  //     std::chrono::duration_cast<std::chrono::microseconds>(end - start)
  //         .count();

  // cout << "Time taken to calculate diff process: " << duration
  //      << " microseconds\n";

  // start = std::chrono::high_resolution_clock::now();
  calc_diff_process(layer, to_process, space, x3, y3, x4, y4);

  // end = std::chrono::high_resolution_clock::now();
  // duration = std::chrono::duration_cast<std::chrono::microseconds>(end -
  // start)
  //                .count();
  // cout << "Time taken to calculate diff process: " << duration
  //      << " microseconds\n";

  // start = std::chrono::high_resolution_clock::now();
  while (!to_process.empty()) {
    // ems_t ems = to_process.front();
    // to_process.pop_front();

    ems_t ems = to_process.back();
    to_process.pop_back();

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
  // end = std::chrono::high_resolution_clock::now();
  // duration = std::chrono::duration_cast<std::chrono::microseconds>(end -
  // start)
  //                .count();
  // cout << "Time taken to process deque " << duration << " microseconds\n";
}

bool sort_rank(const ranking &a, const ranking &b)
{
  return a.chromosome < b.chromosome;
}

void construct_vl_sol(std::vector<ranking> &sol, std::vector<double> chromosome,
                      std::vector<item> items)
{
  std::vector<ranking> rank(chromosome.size());

  unsigned idx = 0;
  for (unsigned i = 0; i < chromosome.size(); i++, idx++) {
    rank[idx].chromosome = chromosome[idx];
    rank[idx].index = i;
    rank[idx].client = items[i].client;
  }

  std::sort(rank.begin(), rank.end(), sort_rank);

  sol = rank;
}

void construct_final_sol(std::vector<ranking> &sol,
                         std::vector<double> chromosome,
                         std::vector<ranking> seq)
{
  std::vector<ranking> rank(chromosome.size());

  unsigned idx = 0;
  for (unsigned i = 0; i < chromosome.size(); i++, idx++) {
    rank[idx].chromosome = chromosome[idx];
    rank[idx].index = seq[idx].index;
    rank[idx].client = seq[idx].client;
  }

  std::sort(rank.begin(), rank.end(), sort_rank);

  sol = rank;
}

void calc_strip_height(unsigned &strip_height, ems_t space,
                       const unsigned &item_height)
{
  unsigned height = space.bottom_point.second + item_height;
  if (height > strip_height) {
    strip_height = height;
  }
}

unsigned pack_with_one_layer(const std::vector<ranking> &rank,
                             const vector<item> &items,
                             const unsigned &max_width, const unsigned &ub,
                             bool debug_sol, std::fstream *solfile)
{
  vector<flat_set<ems_t, bottom_left_cmp>> layers;

  unsigned item_index = rank[0].index;
  item item = items[item_index];
  unsigned strip_height = 0;
  long unsigned int items_placed = 0;
  bool fit;
  unsigned penalty = 0;

  unsigned current_client = item.client;

  FlatSegmentTree seg_tree(max_width);

  // int clients_verification[max_width];
  // std::fill(clients_verification, clients_verification + max_width, -1);

  layers.push_back(flat_set<ems_t, bottom_left_cmp>());
  layers[0].insert({{0, 0}, {max_width, ub}});

  unsigned current_virtual_layer = 0;

  while (items_placed < rank.size()) {
    item_index = rank[items_placed].index;

    item = items[item_index];
    fit = false;

    if (!layers[0].empty()) {
      for (auto ems_t = layers[0].begin(); ems_t != layers[0].end();) {
        if (item_can_fit(item, *ems_t)) {
          calc_strip_height(strip_height, *ems_t, item.height);

          // auto start = std::chrono::high_resolution_clock::now();
          fit_item(item, *ems_t, layers[0], ub, &penalty, seg_tree, debug_sol,
                   solfile);

          // auto end = std::chrono::high_resolution_clock::now();

          // auto duration =
          //     std::chrono::duration_cast<std::chrono::microseconds>(end -
          //     start)
          //         .count();

          // cout << "Time taken to fit item " << item_index << ": " << duration
          //      << " microseconds" << endl;

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

void encode(std::vector<ranking> &rank, const std::vector<ranking> &seq,
            unsigned n)
{
  unsigned idx = 0;
  for (const auto &item : seq) {
    rank[idx].index = item.index;
    rank[idx].chromosome = static_cast<double>(idx) / (n * 10);
    idx++;
  }
}
