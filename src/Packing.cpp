#include "Packing.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <deque>
#include <fstream>
#include <optional>
#include <set>

using namespace std;

typedef struct {
  std::pair<int, int> bottom_point;
  int width;
} mos_t;

struct bottom_left_cmp_open_space {
  bool operator()(const mos_t &a, const mos_t &b) const
  {
    // Primeiro compara por y (bottom_point.second)
    if (a.bottom_point.second != b.bottom_point.second) {
      return a.bottom_point.second < b.bottom_point.second;
    }

    return a.bottom_point.first < b.bottom_point.first;
  }
};

struct dominated_cmp {
  bool operator()(const mos_t &a, const mos_t &b) const
  {
    if (a.bottom_point.first != b.bottom_point.first) {
      return a.bottom_point.first < b.bottom_point.first;  // Ordena por x
    }
    return a.bottom_point.second < b.bottom_point.second;  // Depois por y
  }
};

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

inline bool can_item_fit(const item &item, const mos_t &mos)
{
  return item.width <= mos.width;
}

bool intersects(const mos_t &mos, int space_bp_x, int space_bp_y,
                int space_width, const int &item_width, const int &item_height)
{
  int item_x1 = space_bp_x;
  int item_x2 = space_bp_x + (int)item_width;
  int item_y1 = space_bp_y;
  int item_y2 = space_bp_y + (int)item_height;

  int mos_x1 = mos.bottom_point.first;
  int mos_x2 = mos_x1 + mos.width;
  int mos_y1 = mos.bottom_point.second;
  int mos_y2 = std::numeric_limits<int>::max();

  return (item_x1 < mos_x2 && item_x2 > mos_x1 && item_y1 < mos_y2 &&
          item_y2 > mos_y1);
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

int find_mos(const item &item, std::set<mos_t> &slist,
             FlatSegmentTree &seg_tree)
{
  for (size_t i = 0; i < slist.size(); i++) {
    if (can_item_fit(item, slist[i])) {
      return i;
    }
  }
  return -1;
}

void update_strip_height(int &strip_height, int &mos_index,
                         const int &item_height, const std::set<mos_t> &slist)
{
  mos_t mos = slist[mos_index];
  int height = mos.bottom_point.second + item_height;
  if (height > strip_height) {
    strip_height = height;
  }
}

bool is_dominated(const mos_t &p, const std::vector<mos_t> &layer)
{
  // Get the range of MOS with same x
  auto range = std::equal_range(
      layer.begin(), layer.end(), p, [](const mos_t &a, const mos_t &b) {
        return a.bottom_point.first < b.bottom_point.first;  // compare only x
      });

  for (auto it = range.first; it != range.second; ++it) {
    if (it->bottom_point.second <= p.bottom_point.second &&
        it->width >= p.width) {
      return true;  // p is dominated
    }
  }
  return false;
}

void place_item(const item &item, int &mos_index, std::vector<mos_t> &slist,
                FlatSegmentTree &seg_tree, int *penalty, bool debug_sol = false,
                fstream *solfile = nullptr)
{
  mos_t space = slist[mos_index];
  int space_width = space.width;
  int item_tp_x = space.bottom_point.first + item.width;
  int space_bp_x = space.bottom_point.first;
  int space_bp_y = space.bottom_point.second;

  int min_client = seg_tree.query(space_bp_x, item_tp_x - 1);

  if (min_client < item.client) {
    *penalty += item.client - min_client;
  }

  seg_tree.update(space_bp_x, item_tp_x - 1, item.client);

  if (debug_sol) {
    *solfile << space_bp_x << " " << space.bottom_point.second << "\n"
             << item_tp_x << " " << space.bottom_point.second + item.height
             << "\n";
  }

  slist.erase(slist.begin() + mos_index);

  std::vector<mos_t> to_process;

  mos_t new_space1 = {{space_bp_x, space_bp_y + item.height}, space_width};
  to_process.push_back(new_space1);

  if (space_width - item.width > 0) {
    mos_t new_space2 = {{space_bp_x + item.width, space_bp_y},
                        space_width - item.width};
    to_process.push_back(new_space2);
  }

  for (auto mos = slist.begin(); mos != slist.end();) {
    if (intersects(*mos, space_bp_x, space_bp_y, space_width, item.width,
                   item.height)) {
      int mos_x = mos->bottom_point.first;
      int mos_y = mos->bottom_point.second;
      int mos_width = mos->width;

      mos_t old_space = *mos;
      mos = slist.erase(mos);

      if (mos_x < space_bp_x) {
        mos_t new_space3 = {{mos_x, mos_y}, space_bp_x - mos_x};
        to_process.push_back(new_space3);
      }

      if (mos_x + mos_width > item_tp_x) {
        mos_t new_space4 = {{item_tp_x, mos_y}, mos_x + mos_width - item_tp_x};
        to_process.push_back(new_space4);
      }

      mos_t new_space5 = {{mos_x, space_bp_y + item.height}, mos_width};
      to_process.push_back(new_space5);
    }
    else {
      mos++;
    }
  }

  std::vector<mos_t> slistByX = slist;
  std::sort(slistByX.begin(), slistByX.end(), dominated_cmp());

  for (auto new_space : to_process) {
    if (!is_dominated(new_space, slistByX)) {
      slist.push_back(new_space);
      slistByX.push_back(new_space);
      std::sort(slistByX.begin(), slistByX.end(), dominated_cmp());
    }
  }

  std::sort(slist.begin(), slist.end(), bottom_left_cmp_open_space());
}

int pack_with_one_layer(const std::vector<ranking> &rank,
                        const vector<item> &items, const int &max_width,
                        const int &ub, bool debug_sol, std::fstream *solfile)
{
  set<mos_t> slist = {{{0, 0}, max_width}};
  FlatSegmentTree seg_tree(max_width);

  int item_index = rank[0].index;
  item item = items[item_index];

  int strip_height = 0;
  long int items_placed = 0;
  int penalty = 0;

  int current_client = item.client;

  while (items_placed < rank.size()) {
    item_index = rank[items_placed].index;

    item = items[item_index];

    if (!slist.empty()) {
      int mos_index = find_mos(item, slist, seg_tree);

      if (mos_index == -1) {
        return ub + strip_height;
      }

      update_strip_height(strip_height, mos_index, item.height, slist);

      place_item(item, mos_index, slist, seg_tree, &penalty, debug_sol,
                 solfile);

      items_placed++;
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
