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

int main(int argc, char *argv[])
{
  auto start = std::chrono::high_resolution_clock::now();
  unsigned n;        // size of chromosomes
  unsigned p = 50;   // size of population
  double pe = 0.20;  // fraction of population to be the elite-set
  double pm = 0.1;   // fraction of population to be replaced by mutants
  double rhoe =
      0.70;  // probability that offspring inherit an allele from elite parent
  unsigned K = 3;     // number of independent populations
  unsigned MAXT = 1;  // number of threads for parallel decoding

  int best = std::numeric_limits<int>::max();
  int opt;
  int num_pieces;
  int width, height;
  std::string input_filename;

  while ((opt = getopt(argc, argv, "p:e:m:o:k:t:f:")) != -1) {
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
    }
  }

  std::ifstream input_file(input_filename);

  input_file >> n;

  std::vector<std::pair<int, int>> pieces;

  for (int i = 0; i < n; i++) {
    input_file >> height >> width;

    pieces.push_back({height, width});

    // std::cout << height << " " << width << std::endl;
  }

  const std::string logs_folder = "logs";

  if (!std::filesystem::exists(logs_folder)) {
    std::filesystem::create_directory(logs_folder);
  }

  std::filesystem::path filePath(input_filename);
  std::string logfile_name = filePath.stem().string();
  std::fstream logfile(logs_folder + "/" + logfile_name + ".log",
                       std::ios::app);

  std::stringstream params;
  params << "Número de cromossomos: " << n << "\n"
         << "Tamanho da população: " << p << "\n"
         << "Fração elite: " << pe << "\n"
         << "Fração mutante: " << pm << "\n"
         << "Probabilidade herdar alelo elite: " << rhoe << "\n"
         << "Número de populações independentes: " << K << "\n"
         << "Número de threads: " << MAXT << "\n";

  logfile << params.str();

  SampleDecoder decoder(pieces);  // initialize the decoder

  const long unsigned rngSeed = 0;  // seed to the random number generator
  MTRand rng(rngSeed);              // initialize the random number generator

  // initialize the BRKGA-based heuristic
  BRKGA<SampleDecoder, MTRand> algorithm(n, p, pe, pm, rhoe, decoder, rng, K,
                                         MAXT);

  unsigned generation = 0;  // current generation
  const unsigned X_INTVL =
      10;  // exchange best individuals at every 100 generations
  const unsigned X_NUMBER = 2;    // exchange top 2 best
  const unsigned MAX_GENS = 100;  // run for 1000 gens
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

  logfile << "Melhor altura = " << best << "\n\n";

  logfile << "Tempo de execução: " << duration << "ms" << "\n\n"
          << "---------------------------" << "\n\n";
  return 0;
}