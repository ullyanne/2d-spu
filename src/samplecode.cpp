#include <unistd.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

#include "BRKGA.h"
#include "MTRand.h"
#include "SampleDecoder.h"

using namespace std;

int main(int argc, char *argv[])
{
  unsigned num_items;  // size of chromosomes
  unsigned p = 1000;   // size of population
  double pe = 0.30;    // fraction of population to be the elite-set
  double pm = 0.4;     // fraction of population to be replaced by mutants
  double rhoe =
      0.80;  // probability that offspring inherit an allele from elite parent
  unsigned K = 3;          // number of independent populations
  unsigned MAXT = 2;       // number of threads for parallel decoding
  unsigned X_INTVL = 100;  // exchange best individuals at every 100 generations
  unsigned MAX_GENS = 110;  // run for 1000 gens

  int best = numeric_limits<int>::max();
  int opt;
  string input_filename;

  while ((opt = getopt(argc, argv, "p:e:m:o:k:t:f:b:g:")) != -1) {
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
    }
  }

  ifstream input_file(input_filename);

  int num_clients, max_width, num_items_client, client, width, height;

  input_file >> num_clients >> num_items >> max_width;

  vector<item> items(num_items);

  int index = 0;

  for (int i = 0; i < num_clients; i++) {
    input_file >> client;
    input_file >> num_items_client;

    for (int j = 0; j < num_items_client; j++) {
      input_file >> height >> width;
      items[index].height = height;
      items[index].width = width;
      items[index].client = client;
      index++;
    }
  }

  int val[num_items] = {0};
  int area = 0;
  int lb2 = 0;
  int lb1, lb;

  for (unsigned i = 0; i < num_items; i++) {
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

  if (!std::filesystem::exists(logs_folder)) {
    std::filesystem::create_directory(logs_folder);
  }

  std::filesystem::path filePath(input_filename);
  std::string parentPathStr = filePath.parent_path().string();
  std::size_t pos = parentPathStr.find("instances/");
  std::string relativePath =
      parentPathStr.substr(pos + std::string("instances/").length());
  std::filesystem::path logFilePath = std::filesystem::path(logs_folder) /
                                      relativePath /
                                      (filePath.stem().string() + ".log");

  if (!std::filesystem::exists(logFilePath.parent_path())) {
    std::filesystem::create_directories(logFilePath.parent_path());
  }

  std::fstream logfile(logFilePath, std::ios::app);

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

  SampleDecoder decoder(items, max_width);  // initialize the decoder

  const long unsigned rngSeed = 0;  // seed to the random number generator
  MTRand rng(rngSeed);              // initialize the random number generator

  // initialize the BRKGA-based heuristic
  BRKGA<SampleDecoder, MTRand> algorithm(num_items, p, pe, pm, rhoe, decoder,
                                         rng, K, MAXT);

  unsigned generation = 0;      // current generation
  const unsigned X_NUMBER = 2;  // exchange top 2 best

  auto start = std::chrono::high_resolution_clock::now();

  do {
    algorithm.evolve();  // evolve the population for one generation
    int cost = algorithm.getBestFitness();

    if (best > cost) {
      best = cost;
      // printf("%d %d\n", best, generation);
      logfile << "Melhor até agora: " << best << " "
              << "Geração: " << generation << "\n";
    }

    if ((++generation) % X_INTVL == 0) {
      algorithm.exchangeElite(X_NUMBER);  // exchange top individuals
    }
  } while (generation < MAX_GENS);

  best = algorithm.getBestFitness();

  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();

  std::cout << "Melhor altura = " << best << "\n";

  logfile << "Largura da faixa = " << max_width << "\n";
  logfile << "Melhor altura = " << best << "\n\n";
  logfile << "Lower Bound 1 = " << lb1 << "\n";
  logfile << "Lower Bound 2 = " << lb2 << "\n";
  logfile << "Avaliação = " << (double)best / lb << "\n\n";

  logfile << "Tempo de execução: " << duration << "ms" << "\n\n"
          << "---------------------------" << "\n\n";
  return 0;
}