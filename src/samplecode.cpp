#include <unistd.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <unordered_map>

#include "BRKGA.h"
#include "MTRand.h"
#include "Packing.h"
#include "SampleDecoder.h"
#include "VLBRKGA.h"
#include "VirtualLayersDecoder.h"

using namespace std;
using std::log10;

bool sort_lns(const ranking& a, const ranking& b,
              const std::vector<item>& items)
{
  if (items[a.index].client != items[b.index].client) {
    return items[a.index].client > items[b.index].client;
  }
  return (items[a.index].height * items[a.index].width) >
         (items[b.index].height * items[b.index].width);
}

void choose_client(MTRand& random, unsigned& client, unsigned& next_client,
                   unsigned num_clients)
{
  client = random.randInt(num_clients - 1) + 1;

  if (client == 1 && num_clients == 1) {
    next_client = 1;
  }
  else if (client == 1) {
    next_client = client + 1;
  }
  else if (client == num_clients) {
    next_client = client - 1;
  }
  else {
    unsigned ant_or_suc = random.randInt(1);
    if (ant_or_suc == 0) {
      next_client = client - 1;
    }
    else {
      next_client = client + 1;
    }
  }
}

int main(int argc, char* argv[])
{
  unsigned num_items;  // size of chromosomes
  unsigned p = 50;     // size of population
  double pe = 0.3;     // fraction of population to be the elite-set
  double pm = 0.7;     // fraction of population to be replaced by mutants
  double rhoe =
      0.80;  // probability that offspring inherit an allele from elite parent
  unsigned K = 4;          // number of independent populations
  unsigned MAXT = 14;      // number of threads for parallel decoding
  unsigned X_INTVL = 600;  // exchange best individuals at every 100 generations
  unsigned MAX_GENS = 100;  // run for 1000 gens
  unsigned lns_max_iter = 1000;
  bool debug_sol = false;

  int best_fitness = numeric_limits<int>::max();
  int opt;
  string input_filename;

  while ((opt = getopt(argc, argv, "p:e:m:o:k:t:f:b:g:l:d")) != -1) {
    switch (opt) {
      case 'f':
        input_filename = optarg;
        break;

      case 'p':
        p = std::stoi(optarg);
        break;

      case 'e':
        pe = std::strtod(optarg, nullptr);
        break;

      case 'm':
        pm = std::strtod(optarg, nullptr);
        break;

      case 'o':
        rhoe = std::strtod(optarg, nullptr);
        break;

      case 'k':
        K = std::stoi(optarg);
        break;

      case 't':
        MAXT = std::stoi(optarg);
        break;

      case 'b':
        X_INTVL = std::stoi(optarg);
        break;

      case 'g':
        MAX_GENS = std::stoi(optarg);
        break;

      case 'l':
        lns_max_iter = std::stoi(optarg);
        break;

      case 'd':
        debug_sol = true;
        break;
    }
  }

  ifstream input_file(input_filename);

  unsigned num_clients, max_width, num_items_client, client, width, height;

  input_file >> num_clients >> num_items >> max_width;

  vector<item> items(num_items);
  vector<ranking> lns_seq(num_items);

  int index = 0;

  for (unsigned i = 0; i < num_clients; i++) {
    input_file >> client;
    input_file >> num_items_client;

    for (unsigned j = 0; j < num_items_client; j++) {
      input_file >> height >> width;
      items[index].height = height;
      items[index].width = width;
      items[index].client = client;

      lns_seq[index].index = index;
      lns_seq[index].client = client;

      index++;
    }
  }

  int val[num_items] = {0};
  int area = 0;
  int lb2 = 0;
  int lb1, lb;
  unsigned ub = 0;

  for (unsigned i = 0; i < num_items; i++) {
    ub += items[i].height;
    area += items[i].width * items[i].height;
    val[i] += items[i].height;
    lb2 = max(lb2, val[i]);

    for (unsigned j = i + 1; j < num_items; j++) {
      if (items[j].client != items[i].client &&
          items[i].width + items[j].width > max_width && val[i] > val[j]) {
        val[j] = val[i];
      }
    }
  }

  lb1 = area / max_width;
  lb = max(lb1, lb2);

  const std::string logs_folder = "logs";
  const std::string sol_folder = "sol";

  if (!std::filesystem::exists(logs_folder)) {
    std::filesystem::create_directory(logs_folder);
  }

  if (!std::filesystem::exists(sol_folder)) {
    std::filesystem::create_directory(sol_folder);
  }

  std::filesystem::path filePath(input_filename);
  std::filesystem::path solPath(input_filename);
  std::string parentPathStr = filePath.parent_path().string();
  std::size_t pos = parentPathStr.find("instances/");
  std::string relativePath =
      parentPathStr.substr(pos + std::string("instances/").length());
  std::filesystem::path logFilePath = std::filesystem::path(logs_folder) /
                                      relativePath /
                                      (filePath.stem().string() + ".log");

  std::filesystem::path solFilePath = std::filesystem::path(sol_folder) /
                                      relativePath /
                                      (filePath.stem().string() + "Sol.txt");

  if (!std::filesystem::exists(logFilePath.parent_path())) {
    std::filesystem::create_directories(logFilePath.parent_path());
  }

  if (!std::filesystem::exists(solFilePath.parent_path())) {
    std::filesystem::create_directories(solFilePath.parent_path());
  }

  std::fstream logfile(logFilePath, std::ios::app);
  std::fstream solfile(solFilePath, std::ios::out);

  std::stringstream params;

  std::cout << filePath.stem().string() << "\n";

  auto start = std::chrono::high_resolution_clock::now();

  std::sort(lns_seq.begin(), lns_seq.end(),
            [&items](const ranking& a, const ranking& b) {
              return sort_lns(a, b, items);
            });

  unsigned num_layers = 1;
  unsigned count_layer_items = 0;
  unsigned group_size_per_layer = 20;
  unsigned current_client = lns_seq[0].client;
  // std::vector<std::vector<ranking>> v_layers;
  // v_layers.emplace_back();

  // if (num_items <= 40) {
  //   v_layers[num_layers - 1] = lns_seq;
  // }
  // else {
  //   for (ranking item : lns_seq) {
  //     if (item.client != current_client) {
  //       current_client = item.client;

  //       if (count_layer_items >= group_size_per_layer) {
  //         num_layers++;
  //         count_layer_items = 0;
  //         v_layers.emplace_back();
  //       }
  //     }

  //     v_layers[num_layers - 1].push_back(item);
  //     count_layer_items++;
  //   }
  // }

  unsigned classes_per_layer = num_clients;
  num_layers = ceil(static_cast<double>(num_clients) / classes_per_layer);
  int layer_lim = num_clients - classes_per_layer;

  std::vector<std::vector<ranking>> v_layers(num_layers);
  unsigned idx = 0;
  for (ranking item : lns_seq) {
    if (item.client <= layer_lim) {
      layer_lim -= classes_per_layer;
      if (layer_lim < 0) {
        layer_lim = 0;
      }
      idx++;
    }
    v_layers[idx].push_back(item);
  }

  unordered_map<unsigned, vector<unsigned>> clients_to_layers;
  unordered_map<unsigned, unsigned> layers_to_index;
  unsigned best_height;

  std::vector<ranking> best_sol = lns_seq;

  unsigned n;
  best_height = pack(lns_seq, items, max_width, ub, clients_to_layers,
                     layers_to_index, n);

  cout << best_height << "\n";

  // unsigned num_pieces_per_layer = num_items / log(num_items);

  // if (num_items <= 40) {
  //   num_pieces_per_layer = num_items;
  // }

  // if (num_items >= 200) {
  //   num_pieces_per_layer = num_items / log2(num_items);
  // }

  // unsigned num_virtual_layers =
  //     ceil(static_cast<double>(num_items) / num_pieces_per_layer);

  unsigned bb = best_height;
  std::vector<std::vector<ranking>> tmp(2);
  height =
      pack_with_one_layer(best_sol, items, max_width, ub, tmp, 5, false, bb);

  cout << "h1: " << height << "\n";

  if (height < best_height) {
    best_height = height;
  }
  cout << "h2: " << best_height << "\n";

  unsigned p2 = 300;
  double pe2 = 0.2;
  double pm2 = 0.15;
  double rhoe2 = 0.7;
  unsigned MAX_GENS2 = 250;
  unsigned X_INTVL2 = 40;
  unsigned K2 = 4;
  unsigned max_improvements = 150;

  int index_plot[num_layers];

  params << "Tamanho da população: " << p2 << "\n"
         << "Fração elite: " << pe2 << "\n"
         << "Fração mutante: " << pm2 << "\n"
         << "Probabilidade herdar alelo elite: " << rhoe2 << "\n"
         << "Número de populações independentes: " << K2 << "\n"
         << "Número de threads: " << MAXT << "\n"
         << "Trocar os melhores indivíduos a cada " << X_INTVL2 << " gerações\n"
         << "Número máximo de gerações: " << MAX_GENS2 << "\n"
         << "Número máximo de iterações sem melhora: " << max_improvements
         << "\n";

  logfile << params.str();

  auto startpart2 = std::chrono::high_resolution_clock::now();
  std::chrono::time_point<std::chrono::_V2::system_clock> startiteration;
  std::chrono::time_point<std::chrono::_V2::system_clock> enditeration;

  std::chrono::time_point<std::chrono::_V2::system_clock> startevolution;
  std::chrono::time_point<std::chrono::_V2::system_clock> endevolution;

  std::chrono::time_point<std::chrono::_V2::system_clock> startconstruct;
  std::chrono::time_point<std::chrono::_V2::system_clock> endconstruct;

  std::chrono::time_point<std::chrono::_V2::system_clock> startconstructsol;
  std::chrono::time_point<std::chrono::_V2::system_clock> endconstructsol;

  std::chrono::time_point<std::chrono::_V2::system_clock> startiniti;
  std::chrono::time_point<std::chrono::_V2::system_clock> endiniti;

  for (unsigned i = 0; i < num_layers; i++) {
    unsigned no_improvement = 0;
    if (i == 0) {
      startiteration = std::chrono::high_resolution_clock::now();
    }

    if (v_layers[i].size() == 0) {
      continue;
    }

    cout << "size! " << v_layers[i].size() << "\n";

    VirtualLayersDecoder decoder(v_layers, v_layers[i].size(), i, items,
                                 max_width,
                                 ub);  // initialize the decoder

    const long unsigned rngSeed = 10;  // seed to the random number generator
    MTRand rng2;                       // initialize the random number generator
    startiniti = std::chrono::high_resolution_clock::now();
    // initialize the BRKGA-based heuristic
    VLBRKGA<VirtualLayersDecoder, MTRand> algorithm(
        v_layers[i].size(), p2, pe2, pm2, rhoe2, decoder, rng2, K2, MAXT);

    endiniti = std::chrono::high_resolution_clock::now();
    unsigned generation = 0;      // current generation
    const unsigned X_NUMBER = 2;  // exchange top 2 best

    best_fitness = best_height;
    do {
      startevolution = std::chrono::high_resolution_clock::now();
      algorithm.evolve();  // evolve the population for one generation
      int cost = algorithm.getBestFitness();

      if (cost < best_fitness) {
        best_fitness = cost;
        no_improvement = 0;

        logfile << "Melhor até agora: " << best_fitness << " "
                << "Geração: " << generation << "\n"
                << "Iteração: " << i << "\n";
      }
      else {
        no_improvement++;
      }
      // if (cost < best_fitness) {
      //   best_fitness = cost;
      //   logfile << "Melhor até agora: " << best_fitness << " "
      //           << "Geração: " << generation << "\n"
      //           << "Iteração: " << i << "\n";
      // }

      if ((++generation) % X_INTVL2 == 0) {
        algorithm.exchangeElite(X_NUMBER);  // exchange top individuals
      }
      endevolution = std::chrono::high_resolution_clock::now();
    } while (generation < MAX_GENS2 && no_improvement < max_improvements);

    startconstructsol = std::chrono::high_resolution_clock::now();
    std::vector<ranking> sol;
    construct_vl_sol(sol, algorithm.getBestChromosome(), i, v_layers);
    endconstructsol = std::chrono::high_resolution_clock::now();
    best_sol = sol;
    best_fitness = algorithm.getBestFitness();
    best_height = best_fitness;
    cout << "kk1 " << best_height << "\n" << flush;

    if (i == 0) {
      enditeration = std::chrono::high_resolution_clock::now();
    }

    if (i != 0) {
      unsigned ls_pieces = 70;
      std::vector<ranking> ls_seq =
          slice_layers(v_layers[i], v_layers[i - 1], ls_pieces);
      std::vector<ranking> best_ls_sol = ls_seq;

      std::vector<std::vector<ranking>> remnd(2);
      if (v_layers[i - 1].size() > ls_pieces) {
        remnd[0].insert(remnd[0].end(), v_layers[i - 1].begin(),
                        v_layers[i - 1].end() - ls_pieces);
      }

      if (v_layers[i].size() > ls_pieces) {
        remnd[1].insert(remnd[1].end(), v_layers[i].begin() + ls_pieces,
                        v_layers[i].end());
      }

      std::vector<ranking> base;
      for (int j = 0; j < i - 1; j++) {
        base.insert(base.end(), v_layers[j].begin(), v_layers[j].end());
      }

      double T = 5000;
      vector<ranking> best_full_sol = best_sol;
      double current_fitness = best_fitness;
      double best_fit = current_fitness;
      double max_iterations = 500;
      double T_min = 0.01;
      double alpha = 0.995;
      MTRand ran;
      MTRand sa;

      for (; T > T_min;) {
        std::vector<ranking> neighbor = ls_seq;
        int idx = ran.randInt(neighbor.size() - 1);
        int idx2 = ran.randInt(neighbor.size() - 1);
        while (idx2 == idx) {
          idx2 = ran.randInt(neighbor.size() - 1);
        }
        swap(neighbor[idx], neighbor[idx2]);
        // move_element(neighbor, idx, idx2);

        std::vector<ranking> seq2;
        seq2.insert(seq2.end(), base.begin(), base.end());
        seq2.insert(seq2.end(), remnd[0].begin(), remnd[0].end());
        seq2.insert(seq2.end(), neighbor.begin(), neighbor.end());
        seq2.insert(seq2.end(), remnd[1].begin(), remnd[1].end());

        unsigned new_fitness = ls_pack(seq2, items, max_width, ub, tmp, 5,
                                       false, best_height, false);

        int delta = new_fitness - current_fitness;
        if (delta <= 0 || exp(-delta / T) > sa.rand()) {
          ls_seq = neighbor;
          current_fitness = new_fitness;

          if (new_fitness <= best_fit) {
            best_fit = new_fitness;
            best_ls_sol = neighbor;
            best_full_sol = seq2;
          }
        }

        T *= alpha;
      }

      if (best_fit <= best_fitness) {
        unsigned s1 =
            min(ls_pieces, static_cast<unsigned>(v_layers[i - 1].size()));
        unsigned s2 = min(ls_pieces, static_cast<unsigned>(v_layers[i].size()));

        std::copy_n(best_ls_sol.begin(), s1, v_layers[i - 1].end() - s1);
        std::copy_n(best_ls_sol.begin() + s1, s2, v_layers[i].begin());
        best_fitness = best_fit;
        best_sol = best_full_sol;
        best_height = best_fit;
      }
      index_plot[i - 1] = v_layers[i - 1][v_layers[i - 1].size() - 1].index;
    }

    // if (i != num_layers - 1) {
    //   unsigned n1 = min(50, static_cast<int>(v_layers[i].size()));
    //   unsigned n2 = min(20, static_cast<int>(v_layers[i + 1].size()));

    //   std::vector<ranking> seq;
    //   seq.insert(seq.end(), v_layers[i].end() - n1, v_layers[i].end());
    //   seq.insert(seq.end(), v_layers[i + 1].begin(),
    //              v_layers[i + 1].begin() + n2);
    //   // std::sort(seq.begin(), seq.end(),
    //   //           [&items](const ranking& a, const ranking& b) {
    //   //             return sort_lns(a, b, items);
    //   //           });

    //   cout << "he! \n";
    //   // for (unsigned i = 0; i < seq.size(); i++) {
    //   //   cout << seq[i].client << " " << flush;
    //   //   cout << seq[i].index << "\n" << flush;
    //   // }

    //   std::vector<ranking> l1 =
    //       std::vector<ranking>(v_layers[i].begin(), v_layers[i].end() - n1);

    //   unsigned best = best_height;

    //   std::vector<ranking> global_best_solution;
    //   unsigned global_best_fit = std::numeric_limits<unsigned>::max();

    //   std::vector<ranking> origin;
    //   std::vector<ranking> best_solution = seq;

    //   for (unsigned j = 0; j < 500; j++) {
    //     double T = 1000;
    //     unsigned current_fitness = best;
    //     unsigned best_fit = current_fitness;
    //     unsigned max_iterations = 1000;
    //     double T_min = 0.1;
    //     double alpha = 0.9;
    //     MTRand ran;

    //     for (int iter = 0; iter < max_iterations && T > T_min; ++iter) {
    //       std::vector<ranking> neighbor = seq;
    //       int idx = ran.randInt(neighbor.size() - 1);
    //       int idx2 = ran.randInt(neighbor.size() - 1);
    //       move_element(neighbor, idx, idx2);

    //       std::vector<ranking> seq2;
    //       seq2.insert(seq2.end(), l1.begin(), l1.end());
    //       seq2.insert(seq2.end(), neighbor.begin(), neighbor.end());
    //       seq2.insert(seq2.end(), v_layers[i + 1].begin() + n2,
    //                   v_layers[i + 1].end());

    //       unsigned new_fitness = pack_with_one_layer(
    //           seq2, items, max_width, ub, tmp, 5, false, best_height);

    //       int delta = new_fitness - current_fitness;
    //       if (delta < 0 ||
    //           (delta < 30) &&
    //               exp(-delta / T) > static_cast<double>(ran()) / RAND_MAX) {
    //         seq = neighbor;
    //         current_fitness = new_fitness;

    //         if (new_fitness < best_fit) {
    //           best_fit = new_fitness;
    //           best_solution = neighbor;
    //         }
    //       }
    //     }

    //     if (best_fit < global_best_fit) {
    //       cout << global_best_fit << " " << best_fit << "\n";
    //       global_best_fit = best_fit;
    //       global_best_solution = best_solution;
    //     }
    //   }
    //   best_sol = global_best_solution;
    //   break;
    //   int j = 0;
    //   for (auto it = v_layers[i].end() - n1; it != v_layers[i].end(); ++it) {
    //     *it = global_best_solution[j];
    //     j++;
    //   }
    //   v_layers[i] = global_best_solution;
    // }
  }

  for (unsigned i = 0; i < num_layers; i++) {
    cout << "i l - " << index_plot[i] << "\n";
  }

  auto endpart2 = std::chrono::high_resolution_clock::now();
  auto durationpart2 = std::chrono::duration_cast<std::chrono::milliseconds>(
                           endpart2 - startpart2)
                           .count();

  auto durationiteration =
      std::chrono::duration_cast<std::chrono::milliseconds>(enditeration -
                                                            startiteration)
          .count();

  auto durationevolution =
      std::chrono::duration_cast<std::chrono::milliseconds>(endevolution -
                                                            startevolution)
          .count();

  auto durationconstruct =
      std::chrono::duration_cast<std::chrono::milliseconds>(endconstruct -
                                                            startconstruct)
          .count();

  cout << "endpart2 " << durationpart2 << "\n";
  cout << "enditeration " << durationiteration << "\n";
  cout << "endevolution " << durationevolution << "\n";
  cout << "endconstruct " << durationconstruct << "\n";
  cout << "endconstructsol "
       << std::chrono::duration_cast<std::chrono::milliseconds>(
              endconstructsol - startconstructsol)
              .count()
       << "\n";
  cout << "endiniti "
       << std::chrono::duration_cast<std::chrono::milliseconds>(endiniti -
                                                                startiniti)
              .count();

  cout << "kk " << best_height << "\n";

  // for (int i = 0; i < best_sol.size(); i++) {
  //   cout << best_sol[i].index << " ";
  // }

  // double T = 1000;
  // double current_fitness = best_height;

  // std::vector<ranking> chromosome = best_sol;
  // std::vector<ranking> best_solution = chromosome;
  // double best_fit = current_fitness;
  // unsigned max_iterations = 1000;
  // double T_min = 0.1;
  // double alpha = 0.95;

  // for (int iter = 0; iter < max_iterations && T > T_min; ++iter) {
  //   // Escolhe um tipo de perturbação aleatório
  //   std::vector<ranking> neighbor = chromosome;

  //   // Avalia a nova solução
  //   double new_fitness =
  //       pack_with_one_layer(neighbor, items, max_width, ub, virtual_layers,
  //                           num_pieces_per_layer, false, best_height);

  //   // Critério de aceitação
  //   double delta = new_fitness - current_fitness;
  //   if (delta < 0 || exp(-delta / T) > static_cast<double>(rand()) /
  //   RAND_MAX) {
  //     chromosome = neighbor;
  //     current_fitness = new_fitness;

  //     if (new_fitness < best_fit) {
  //       best_fit = new_fitness;
  //       best_solution = neighbor;
  //     }
  //   }

  //   // Resfriamento mais suave
  //   T *= alpha;
  // }

  // chromosome = best_solution;

  // cout << best_fit << "\n";

  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();

  debug_sol = true;
  if (debug_sol) {
    solfile << max_width << "\n";

    cout << best_height << "\n";

    std::vector<std::vector<ranking>> vl;
    unsigned bh = 0;
    unsigned dummy = pack_with_one_layer(best_sol, items, max_width, ub, vl, 2,
                                         false, bh, debug_sol, &solfile);

    // else {
    //   unsigned dummy = pack(best_sol, items, max_width, ub,
    //   clients_to_layers,
    //                         layers_to_index, num_layers, debug_sol,
    //                         &solfile);
    // }
  }

  std::cout << "\nMelhor altura = " << best_height << "\n";
  logfile << "\nGrupo de itens por camada = " << group_size_per_layer << "\n";
  logfile << "Largura da faixa = " << max_width << "\n";
  logfile << "Melhor altura = " << best_height << "\n\n";
  logfile << "Lower Bound 1 = " << lb1 << "\n";
  logfile << "Lower Bound 2 = " << lb2 << "\n";
  logfile << "Avaliação = " << (double)best_height / lb << "\n\n";

  logfile << "Tempo de execução: " << duration / 1000 << "s" << "\n\n"
          << "---------------------------" << "\n\n";

  cout << best_height << "\n";
  return 0;
}