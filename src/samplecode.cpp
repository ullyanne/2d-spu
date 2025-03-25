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
  if (a.client != b.client) {
    return a.client > b.client;
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
  params << "Número de cromossomos: " << num_items << "\n"
         << "Tamanho da população: " << p << "\n"
         << "Fração elite: " << pe << "\n"
         << "Fração mutante: " << pm << "\n"
         << "Probabilidade herdar alelo elite: " << rhoe << "\n"
         << "Número de populações independentes: " << K << "\n"
         << "Número de threads: " << MAXT << "\n"
         << "Trocar os melhores indivíduos a cada " << X_INTVL << " gerações\n"
         << "Número máximo de gerações: " << MAX_GENS << "\n"
         << "\n";

  logfile << params.str();

  std::cout << filePath.stem().string() << "\n";

  auto start = std::chrono::high_resolution_clock::now();

  std::sort(lns_seq.begin(), lns_seq.end(),
            [&items](const ranking& a, const ranking& b) {
              return sort_lns(a, b, items);
            });

  // for (const auto& item : lns_seq) {
  //   std::cout << item.client << " " << item.index << " "
  //             << items[item.index].width << " " << items[item.index].height
  //             << "\n";
  // }

  unordered_map<unsigned, vector<unsigned>> clients_to_layers;
  unordered_map<unsigned, unsigned> layers_to_index;
  unsigned best_height;

  std::vector<ranking> best_sol = lns_seq;
  unsigned num_layers;
  best_height = pack(lns_seq, items, max_width, ub, clients_to_layers,
                     layers_to_index, num_layers);

  unsigned swap_num = 0;
  // cout << "first! " << best_height << "\n";
  // MTRand random;
  // MTRand random_lns;
  // MTRand choose_two;

  // unsigned swap_num = 10;
  // unsigned max_swaps = 2;
  // unsigned attempts = 0;

  // /* choose layers of adjacent clients to swap */

  // for (unsigned i = 0; i < swap_num; ++i) {
  //   if (attempts >= max_swaps) {
  //     break;
  //   }

  //   height = best_height;
  //   lns_seq = best_sol;
  //   unsigned client;
  //   unsigned next_client;
  //   choose_client(random, client, next_client, num_clients);

  //   unsigned c_size = clients_to_layers[client].size();

  //   unsigned layer_1 =
  //       clients_to_layers[client][random_lns.randInt(c_size - 1)];

  //   unsigned cn_size = clients_to_layers[next_client].size();

  //   unsigned layer_2 =
  //       clients_to_layers[next_client][random_lns.randInt(cn_size - 1)];

  //   vector<unsigned> subchromosome;

  //   unsigned index_end = lns_seq.size();
  //   num_layers = layers_to_index.size();

  //   if (layer_1 < num_layers - 1) {
  //     index_end = layers_to_index[layer_1 + 1];
  //   }
  //   for (unsigned i = layers_to_index[layer_1]; i < index_end; ++i) {
  //     subchromosome.push_back(i);
  //   }

  //   if (layer_1 != layer_2) {
  //     index_end = lns_seq.size();
  //     if (layer_2 < num_layers - 1) {
  //       index_end = layers_to_index[layer_2 + 1];
  //     }
  //     for (unsigned i = layers_to_index[layer_2]; i < index_end; ++i) {
  //       subchromosome.push_back(i);
  //     }
  //   }

  //   // std::cout << "l! " << layer_1 << " " << layer_2 << "\n";
  //   unsigned subchromosome_size = subchromosome.size();
  //   if (subchromosome_size == 0) {
  //     continue;
  //   }

  //   SampleDecoder decoder(lns_seq, subchromosome, items, max_width,
  //                         ub);  // initialize the decoder

  //   const long unsigned rngSeed = 0;  // seed to the random number generator
  //   MTRand rng;                       // initialize the random number
  //   generator

  //   // initialize the BRKGA-based heuristic
  //   BRKGA<SampleDecoder, MTRand> algorithm(subchromosome_size, p, pe, pm,
  //   rhoe,
  //                                          decoder, rng, K, MAXT);

  //   unsigned generation = 0;      // current generation
  //   const unsigned X_NUMBER = 2;  // exchange top 2 best

  //   best_fitness = best_height;
  //   do {
  //     algorithm.evolve();  // evolve the population for one generation
  //     int cost = algorithm.getBestFitness();

  //     if (cost < best_fitness) {
  //       best_fitness = cost;
  //       // printf("%d %d\n", best, generation);
  //       logfile << "Melhor até agora: " << best_fitness << " "
  //               << "Geração: " << generation << "\n"
  //               << "Iteração: " << i << "\n";
  //     }

  //     if ((++generation) % X_INTVL == 0) {
  //       algorithm.exchangeElite(X_NUMBER);  // exchange top individuals
  //     }
  //   } while (generation < MAX_GENS);

  //   if (best_fitness < best_height) {
  //     best_height = best_fitness;
  //     /* update clients to index, layers to index etc... */
  //     /* remember subchromosome is just a part of the sol. construct sol
  //      * merges lns_seq with subchromosome*/
  //     layers_to_index = algorithm.getBestLayersToIndex();
  //     clients_to_layers = algorithm.getBestClientsToLayers();
  //     construct_sol(best_sol, algorithm.getBestChromosome(), subchromosome,
  //                   lns_seq);
  //     attempts = 0;
  //   }
  //   else {
  //     attempts++;
  //   }
  // }

  cout << best_height << "\n";

  unsigned best_swap;
  std::vector<ranking> best_swap_sol = best_sol;

  unsigned num_pieces_per_layer = num_items / log(num_items);

  if (num_items <= 40) {
    num_pieces_per_layer = num_items;
  }

  if (num_items >= 200) {
    num_pieces_per_layer = num_items / log2(num_items);
  }

  unsigned num_virtual_layers =
      ceil(static_cast<double>(num_items) / num_pieces_per_layer);

  unsigned bb = best_height;
  std::vector<std::vector<ranking>> virtual_layers(num_virtual_layers);
  height = pack_with_one_layer(best_sol, items, max_width, ub, virtual_layers,
                               num_pieces_per_layer, true, bb);

  cout << "h1: " << height << "\n";

  if (height < best_height) {
    best_height = height;
  }
  cout << "h2: " << best_height << "\n";

  cout << virtual_layers[0].size() << "\n";

  unsigned p2 = 250;
  double pe2 = 0.15;
  double pm2 = 0.2;
  double rhoe2 = 0.55;
  unsigned MAX_GENS2 = 100;
  unsigned X_INTVL2 = 40;
  unsigned K2 = 3;
  unsigned max_improvements = 50;

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
  for (unsigned i = 0; i < num_virtual_layers; i++) {
    unsigned no_improvement = 0;
    if (i == 0) {
      startiteration = std::chrono::high_resolution_clock::now();
    }

    if (virtual_layers[i].size() == 0) {
      continue;
    }

    // cout << "first layer " << "\n";
    // for (unsigned k = 0; k < virtual_layers[i].size(); k++) {
    //   cout << virtual_layers[i][k].index << " ";
    // }

    // cout << "second layer " << "\n";
    // for (unsigned k = 0; k < virtual_layers[i + 1].size(); k++) {
    //   cout << virtual_layers[i + 1][k].index << " ";
    // }

    VirtualLayersDecoder decoder(virtual_layers, num_pieces_per_layer, i, items,
                                 max_width,
                                 ub);  // initialize the decoder

    const long unsigned rngSeed = 10;  // seed to the random number generator
    MTRand rng2;                       // initialize the random number generator
    startiniti = std::chrono::high_resolution_clock::now();
    // initialize the BRKGA-based heuristic
    VLBRKGA<VirtualLayersDecoder, MTRand> algorithm(
        virtual_layers[i].size(), p2, pe2, pm2, rhoe2, decoder, rng2, K2, MAXT);

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

    best_fitness = algorithm.getBestFitness();
    if (best_fitness < best_height) {
      /* update clients to index, layers to index etc... */
      /* remember subchromosome is just a part of the sol. construct sol
       * merges lns_seq with subchromosome*/
      startconstructsol = std::chrono::high_resolution_clock::now();
      std::vector<ranking> sol;
      construct_vl_sol(sol, algorithm.getBestChromosome(), i, virtual_layers);

      endconstructsol = std::chrono::high_resolution_clock::now();

      if (sol != best_sol) {
        startconstruct = std::chrono::high_resolution_clock::now();
        best_height = best_fitness;
        cout << "kk1 " << best_height << "\n";
        best_sol = sol;
        std::vector<std::vector<ranking>> copy_virtual_layers(
            num_virtual_layers);
        height = pack_with_one_layer(best_sol, items, max_width, ub,
                                     copy_virtual_layers, num_pieces_per_layer,
                                     true, best_height);
        virtual_layers = copy_virtual_layers;
        endconstruct = std::chrono::high_resolution_clock::now();
      }
    }

    if (i == 0) {
      enditeration = std::chrono::high_resolution_clock::now();
    }
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

  // unsigned no_improvement = 0;

  // while (true) {
  //   cout << no_improvement << "@\n";
  //   cout << best_height << "\n";

  //   if (no_improvement >= 20) {
  //     break;
  //   }

  //   if (best_swap < best_height) {
  //     best_height = best_swap;
  //   }

  //   for (unsigned i = 0; i < num_virtual_layers; i++) {
  //     if (virtual_layers[i].size() == 0) {
  //       continue;
  //     }

  //     VirtualLayersDecoder decoder(virtual_layers, i, items,
  //     max_width,
  //                                  ub);  // initialize the decoder

  //     const long unsigned rngSeed = 0;  // seed to the random number
  //     generator MTRand rng(rngSeed);  // initialize the random number
  //     generator

  //     // initialize the BRKGA-based heuristic
  //     VLBRKGA<VirtualLayersDecoder, MTRand> algorithm(
  //         virtual_layers[i].size(), 200, pe, pm, rhoe, decoder, rng,
  //         4, MAXT);

  //     unsigned generation = 0;      // current generation
  //     const unsigned X_NUMBER = 2;  // exchange top 2 best

  //     best_fitness = best_height;
  //     do {
  //       algorithm.evolve();  // evolve the population for one
  //       generation int cost = algorithm.getBestFitness();

  //       if (cost < best_fitness) {
  //         best_fitness = cost;
  //         // printf("%d %d\n", best, generation);
  //         logfile << "Melhor até agora: " << best_fitness << " "
  //                 << "Geração: " << generation << "\n"
  //                 << "Iteração: " << i << "\n";
  //       }

  //       if ((++generation) % X_INTVL == 0) {
  //         algorithm.exchangeElite(X_NUMBER);  // exchange top
  //         individuals
  //       }
  //     } while (generation < 200);

  //     if (best_fitness < best_height) {
  //       no_improvement = 0;
  //       best_height = best_fitness;
  //       best_swap = best_height;
  //       cout << "kk1 " << best_height << "\n";
  //       /* update clients to index, layers to index etc... */
  //       /* remember subchromosome is just a part of the sol.
  //       construct sol
  //        * merges lns_seq with subchromosome*/
  //       std::vector<ranking> sol;
  //       construct_vl_sol(sol, algorithm.getBestChromosome(), i,
  //       virtual_layers); best_sol = sol;
  //     }
  //     else {
  //       no_improvement++;
  //     }
  //   }
  // }

  cout << "kk " << best_height << "\n";

  // for (int i = 0; i < best_sol.size(); i++) {
  //   cout << best_sol[i].index << " ";
  // }

  // cout << best_height << "\n";
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
  //   int i = rand() % neighbor.size();

  //   int perturbation_type = rand() % 3;
  //   switch (perturbation_type) {
  //     case 0:  // Swap aleatório
  //       std::swap(neighbor[i], neighbor[rand() % neighbor.size()]);
  //       break;
  //     case 1:  // Swap com vizinho
  //       std::swap(neighbor[i], neighbor[(i + 1) % neighbor.size()]);
  //       break;
  //     case 2:  // Inverter uma sublista aleatória
  //     {
  //       int j = rand() % neighbor.size();
  //       if (i > j) std::swap(i, j);
  //       std::reverse(neighbor.begin() + i, neighbor.begin() + j);
  //       break;
  //     }
  //     case 3:  // Trocar um valor por um completamente aleatório
  //       neighbor[i].chromosome = static_cast<double>(rand()) / RAND_MAX;
  //       break;
  //   }

  //   // Avalia a nova solução
  //   double new_fitness =
  //       pack_with_one_layer(neighbor, items, max_width, ub, virtual_layers,
  //                           num_pieces_per_layer, false, best_height);

  //   // Critério de aceitação
  //   double delta = new_fitness - current_fitness;
  //   if (delta < 0 || exp(-delta / T) > static_cast<double>(rand()) / RAND_MAX) {
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
    unsigned dummy = pack_with_one_layer(best_sol, items, max_width, ub, vl,
                                         num_pieces_per_layer, false, bh,
                                         debug_sol, &solfile);

    // else {
    //   unsigned dummy = pack(best_sol, items, max_width, ub,
    //   clients_to_layers,
    //                         layers_to_index, num_layers, debug_sol,
    //                         &solfile);
    // }
  }

  std::cout << "Melhor altura = " << best_height << "\n";

  logfile << "Número de swaps entre camadas = " << swap_num << "\n";
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