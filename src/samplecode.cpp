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

bool sort_by_client_and_area_desc(const ranking& a, const ranking& b,
                                  const std::vector<item>& items)
{
  if (items[a.index].client != items[b.index].client) {
    return items[a.index].client > items[b.index].client;
  }
  return (items[a.index].height * items[a.index].width) >
         (items[b.index].height * items[b.index].width);
}

bool sort_by_client_and_width_desc(const ranking& a, const ranking& b,
                                   const std::vector<item>& items)
{
  if (items[a.index].client != items[b.index].client) {
    return items[a.index].client > items[b.index].client;
  }
  return (items[a.index].width) > (items[b.index].width);
}

bool sort_by_client_and_height_desc(const ranking& a, const ranking& b,
                                    const std::vector<item>& items)
{
  if (items[a.index].client != items[b.index].client) {
    return items[a.index].client > items[b.index].client;
  }
  return (items[a.index].height) > (items[b.index].height);
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
  vector<ranking> seq(num_items);

  int index = 0;

  for (unsigned i = 0; i < num_clients; i++) {
    input_file >> client;
    input_file >> num_items_client;

    for (unsigned j = 0; j < num_items_client; j++) {
      input_file >> height >> width;
      items[index].height = height;
      items[index].width = width;
      items[index].client = client;

      seq[index].index = index;
      seq[index].client = client;

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

  vector<vector<ranking>> initial_seqs(3);

  std::sort(seq.begin(), seq.end(),
            [&items](const ranking& a, const ranking& b) {
              return sort_by_client_and_area_desc(a, b, items);
            });

  initial_seqs[0] = seq;

  std::sort(seq.begin(), seq.end(),
            [&items](const ranking& a, const ranking& b) {
              return sort_by_client_and_width_desc(a, b, items);
            });

  initial_seqs[1] = seq;

  std::sort(seq.begin(), seq.end(),
            [&items](const ranking& a, const ranking& b) {
              return sort_by_client_and_height_desc(a, b, items);
            });

  initial_seqs[2] = seq;

  unordered_map<unsigned, vector<unsigned>> clients_to_layers;
  unordered_map<unsigned, unsigned> layers_to_index;
  unsigned best_height;

  std::vector<ranking> best_sol = seq;

  unsigned n;
  best_height =
      pack(seq, items, max_width, ub, clients_to_layers, layers_to_index, n);

  cout << best_height << "\n";

  height = pack_with_one_layer(best_sol, items, max_width, ub, false);

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
  unsigned K2 = 3;
  unsigned max_improvements = 90;
  unsigned no_improvement = 0;

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

  MTRand ran;

  VirtualLayersDecoder decoder(initial_seqs, items, max_width, ub,
                               ran);  // initialize the decoder

  const long unsigned rngSeed = 10;  // seed to the random number generator
  MTRand rng2;                       // initialize the random number generator

  // initialize the BRKGA-based heuristic
  VLBRKGA<VirtualLayersDecoder, MTRand> algorithm(
      num_items, p2, pe2, pm2, rhoe2, decoder, rng2, K2, MAXT);

  unsigned generation = 0;      // current generation
  const unsigned X_NUMBER = 2;  // exchange top 2 best

  best_fitness = best_height;
  do {
    algorithm.evolve();  // evolve the population for one generation
    int cost = algorithm.getBestFitness();

    if (cost < best_fitness) {
      best_fitness = cost;
      no_improvement = 0;

      logfile << "Melhor até agora: " << best_fitness << " "
              << "Geração: " << generation << "\n";
    }
    else {
      no_improvement++;
    }

    if ((++generation) % X_INTVL2 == 0) {
      algorithm.exchangeElite(X_NUMBER);  // exchange top individuals
    }
  } while (generation < MAX_GENS2 && no_improvement < max_improvements);

  std::vector<ranking> sol;
  construct_vl_sol(sol, algorithm.getBestChromosome(), items);
  best_sol = sol;
  best_fitness = algorithm.getBestFitness();
  best_height = best_fitness;
  std::cout << "kk1 " << best_height << "\n" << flush;

  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();

  debug_sol = true;
  if (debug_sol) {
    solfile << max_width << "\n";

    cout << best_height << "\n";

    unsigned dummy = pack_with_one_layer(best_sol, items, max_width, ub,
                                         debug_sol, &solfile);
  }

  std::cout << "Melhor altura = " << best_height << "\n";
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