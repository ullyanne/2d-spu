#include "Packing.h"

#include <algorithm>

using namespace std;

void open_new_layer(vector<set<ems, bottom_left_cmp>> &layers,
                    unsigned &num_layers, unsigned &strip_height, int max_width,
                    int item_height)
{
  layers.push_back(set<ems, bottom_left_cmp>());
  num_layers++;

  ems space;
  space.bottom_point = make_pair(0, strip_height);
  space.top_point = make_pair(max_width, item_height + strip_height);

  layers[num_layers - 1].insert(space);
  strip_height += item_height;
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
    if (is_contained(new_space, space) || is_contained(space, new_space)) {
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

  if (space.bottom_point.second > space.top_point.second) {
    return false;
  }

  if (space.bottom_point.first > space.top_point.first) {
    return false;
  }

  return true;
}

void calc_diff_process(set<ems, bottom_left_cmp> &layer, ems space, int x3,
                       int y3, int x4, int y4)
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
      layer.insert(space);
    }
  }
}

bool does_intersect(int x3, int y3, int x4, int y4, ems space)
{
  return (x3 < space.top_point.first && x4 > space.bottom_point.first &&
          y3 < space.top_point.second && y4 > space.bottom_point.second);
}

void fit_item(item item, ems space, set<ems, bottom_left_cmp> &layer,
              int clients_verification[], unsigned ub, unsigned *penalty)
{
  layer.erase(space);

  int x3 = space.bottom_point.first;
  int y3 = space.bottom_point.second;
  int x4 = x3 + item.width;
  int y4 = y3 + item.height;

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

  calc_diff_process(layer, space, x3, y3, x4, y4);

  for (auto ems = layer.begin(); ems != layer.end();) {
    if (does_intersect(x3, y3, x4, y4, *ems)) {
      auto ems_tmp = *ems;
      ems = layer.erase(ems);
      calc_diff_process(layer, ems_tmp, x3, y3, x4, y4);
    }
    else {
      ems++;
    }
  }
}

void rearrangeSeq(std::vector<ranking> &arr,
                  const std::vector<unsigned> &indices_to_move)
{
  int n = arr.size();
  std::vector<ranking> result(n);

  // Inserir elementos especificados na nova ordem
  for (size_t i = 0; i < indices_to_move.size(); ++i) {
    result[i] = arr[indices_to_move[i]];
  }

  // Inserir os elementos restantes
  int insert_index = indices_to_move.size();
  for (int i = 0; i < n; ++i) {
    if (std::find(indices_to_move.begin(), indices_to_move.end(), i) ==
        indices_to_move.end()) {
      result[insert_index++] = arr[i];
    }
  }

  // Substituir o array original pelo reorganizado
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

  for (int i = 0; i < chromosome.size(); ++i) {
    std::cout << "rank: " << rank[i].index << "\n";
  }

  full_solution = lns_seq;

  std::vector<unsigned> subchromosome_copy = subchromosome;
  for (unsigned i = 0; i < chromosome.size(); i++) {
    // std::cout << lns_seq_copy[i].index << "\n";
    // std::cout << "bom dia!" << subchromosome[i] << "\n";
    // std::cout << "b!" << rank[i].index << "\n";
    // std::cout << "d!" << lns_seq[rank[i].index].index << "\n";
    subchromosome_copy[i] = rank[i].index;
    // lns_seq_copy[subchromosome[i]] = lns_seq[rank[i].index];
  }

  rearrangeSeq(full_solution, subchromosome_copy);

  for (int i = 0; i < full_solution.size(); i++) {
    // std::cout << "hi!\n";
    std::cout << full_solution[i].index << " ";
  }
}

unsigned pack(
    const std::vector<ranking> &rank, vector<item> items,
    const unsigned max_width, const unsigned ub,
    std::unordered_map<unsigned, std::vector<unsigned>> &clients_to_layers,
    std::unordered_map<unsigned, unsigned> &layers_to_index,
    unsigned &num_layers)
{
  vector<set<ems, bottom_left_cmp>> layers;

  unsigned item_index = rank[0].index;
  item item = items[item_index];
  unsigned strip_height = 0;
  num_layers = 0;
  long unsigned int items_placed = 0;
  bool fit;
  unsigned penalty = 0;

  std::vector<unsigned> layers_client;
  unsigned current_client = 1;

  int clients_verification[max_width];

  for (unsigned i = 0; i < max_width; i++) {
    clients_verification[i] = -1;
  }

  // for (int i = 0; i < rank.size(); i++) {
  //   std::cout << "rank: " << rank[i].index << "\n";
  // }

  open_new_layer(layers, num_layers, strip_height, max_width, item.height);
  layers_to_index[num_layers - 1] = 0;
  layers_client.push_back(num_layers - 1);

  while (items_placed < items.size()) {
    item_index = rank[items_placed].index;
    item = items[item_index];
    fit = false;

    if (!layers[num_layers - 1].empty()) {
      for (auto ems = layers[num_layers - 1].begin();
           ems != layers[num_layers - 1].end();) {
        if (item_can_fit(item, *ems)) {
          fit_item(item, *ems, layers[num_layers - 1], clients_verification, ub,
                   &penalty);
          // cout << "woops " << item_index << "\n";
          fit = true;
          items_placed++;

          if (item.client != current_client) {
            clients_to_layers[current_client] = layers_client;

            // for (unsigned lay : layers_client) {
            //   std::cout << "lay changed!: " << lay << "\n";
            // }

            layers_client.clear();

            // for (unsigned lay : layers_client) {
            //   std::cout << "lay cleared!: " << lay << "\n";
            // }

            // for (unsigned lay : clients_to_layers[current_client]) {
            //   std::cout << "new vec!: " << lay << "\n";
            // }
            current_client = item.client;
            layers_client.push_back(num_layers - 1);
          }
          break;
        }
        else {
          ems++;
        }
      }
    }

    if (!fit) {
      if (item.client != current_client) {
        clients_to_layers[current_client] = layers_client;

        // for (unsigned lay : layers_client) {
        //   std::cout << "lay changed!: " << lay << "\n";
        // }

        layers_client.clear();

        // for (unsigned lay : layers_client) {
        //   std::cout << "lay cleared!: " << lay << "\n";
        // }

        // for (unsigned lay : clients_to_layers[current_client]) {
        //   std::cout << "new vec!: " << lay << "\n";
        // }
        current_client = item.client;
        layers_client.push_back(num_layers - 1);
      }
      else {
        layers_client.push_back(num_layers - 1);
      }
      open_new_layer(layers, num_layers, strip_height, max_width, item.height);
      auto new_space = *layers[num_layers - 1].begin();
      fit_item(item, new_space, layers[num_layers - 1], clients_verification,
               ub, &penalty);
      // cout << "woops1 " << item_index << "\n";
      fit = true;
      items_placed++;
      layers_to_index[num_layers - 1] = items_placed - 1;
    }
  }

  clients_to_layers[current_client] = layers_client;

  if (penalty) {
    penalty += ub;
  }

  // for (int i = 0; i < rank.size(); i++) {
  //   std::cout << "rank: " << rank[i].index << "\n";
  // }

  // std::cout << strip_height + penalty << std::endl;
  return strip_height + penalty;
}