#include <unistd.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <unordered_map>

#include "BRKGA.h"
#include "Decoder.h"
#include "MTRand.h"
#include "Packing.h"
#include "random_manager.h"

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

void show_help(const char* progname)
{
  cout << "Usage: " << progname << " [options]\n"
       << "Options:\n"
       << "  -f <file>      Specifies the path of the input instance\n"
       << "  -p <number>    Number of individuals in each population\n"
       << "  -g <number>    Frequency (in generations) of elite migration "
          "between populations\n"
       << "  -n <number>    Number of elite individuals transferred between "
          "populations\n"
       << "  -e <float>     Fraction of the population formed by elite\n"
       << "  -m <float>     Fraction of the population formed by mutants\n"
       << "  -o <float>     Probability of inheriting an allele from the elite "
          "parent (rho_e)\n"
       << "  -k <number>    Number of independent populations (K)\n"
       << "  -t <number>    Sets the time limit for execution (defaults to "
          "60s)\n"
       << "  -d             Activates debug mode, the output is a file you can "
          "use as input in debug.py script\n"
       << "  -x <float>     Fraction of the initial population from naive "
          "solutions (gamma_1)\n"
       << "  -y <float>     Fraction of the initial population from naive "
          "solutions (gamma_2)\n"
       << "  -z <float>     Fraction of the initial population from naive "
          "solutions (gamma_3)\n"
       << "  -i <number>    Swap window size used in initial solution "
          "generation (w_i)\n"
       << "  -l <number>    Swap window size used in RLS (w_RLS)\n"
       << "  -v <float>     Probability of applying RLS to an elite individual "
          "(rho_RLS)\n"
       << "  -s <number>    Sets the seed for generating random numbers\n"
       << "  -a <number>    Number of RLS iterations (i_RLS)\n"
       << "  -h             Shows this help menu\n";
}

int main(int argc, char* argv[])
{
  unsigned num_items;  // size of chromosomes
  unsigned p = 300;    // size of population
  double pe = 0.2;     // fraction of population to be the elite-set
  double pm = 0.15;    // fraction of population to be replaced by mutants
  double rhoe =
      0.7;  // probability that offspring inherit an allele from elite parent
  unsigned K = 3;         // number of independent populations
  unsigned MAXT = 1;      // number of threads for parallel decoding
  unsigned X_INTVL = 30;  // exchange best individuals at every 100 generations
                          // // run for 1000 gens
  bool debug_sol = false;

  int opt;
  string input_filename;
  unsigned global_seed;
  unsigned attempts;
  double i1, i2, i3;
  double window_init, window_ls;
  double top_elite;
  unsigned max_improvements;
  unsigned time_limit = 60;
  unsigned X_NUMBER;

  while ((opt = getopt(argc, argv, "p:g:n:e:m:o:k:t:f:s:a:x:y:z:i:l:v:dh")) !=
         -1) {
    switch (opt) {
      case 'f':
        input_filename = optarg;
        break;

      case 'p':
        p = std::stoi(optarg);
        break;

      case 'g':
        X_INTVL = std::stoi(optarg);
        break;

      case 'n':
        X_NUMBER = std::stoi(optarg);
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
        time_limit = std::stoi(optarg);
        break;

      case 'd':
        debug_sol = true;
        break;

      case 'x':
        i1 = std::strtod(optarg, nullptr);
        break;

      case 'y':
        i2 = std::strtod(optarg, nullptr);
        break;

      case 'z':
        i3 = std::strtod(optarg, nullptr);
        break;

      case 'i':
        window_init = std::strtod(optarg, nullptr);
        break;

      case 'l':
        window_ls = std::strtod(optarg, nullptr);
        break;

      case 'v':
        top_elite = std::strtod(optarg, nullptr);
        break;

      case 's':
        global_seed = std::stoi(optarg);
        set_global_seed(global_seed);
        break;

      case 'a':
        attempts = std::stoi(optarg);
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

  if (!std::filesystem::exists(logs_folder)) {
    std::filesystem::create_directory(logs_folder);
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

  if (!std::filesystem::exists(logFilePath.parent_path())) {
    std::filesystem::create_directories(logFilePath.parent_path());
  }

  std::fstream logfile(logFilePath, std::ios::app);

  std::stringstream params;

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

  std::vector<std::vector<ranking>> rank_groups(
      3, std::vector<ranking>(num_items));

  for (unsigned j = 0; j < 3; j++) {
    encode(rank_groups[j], initial_seqs[j], num_items);
    std::sort(
        rank_groups[j].begin(), rank_groups[j].end(),
        [](const ranking& a, const ranking& b) { return a.index < b.index; });
  }

  unsigned best_height = numeric_limits<int>::max();
  std::vector<ranking> best_sol = seq;

  Decoder decoder(rank_groups, items, max_width,
                  ub);  // initialize the decoder

  MTRand rng(get_seed(0));  // initialize the random number generator
  MTRand r_ls(get_seed(1));
  MTRand r_init(get_seed(2));
  MTRand r_elite(get_seed(3));

  // initialize the BRKGA-based heuristic
  BRKGA<Decoder, MTRand> algorithm(num_items, p, pe, pm, rhoe, decoder, rng,
                                   r_ls, r_init, r_elite, attempts, i1, i2, i3,
                                   window_init, window_ls, top_elite, K, MAXT);

  unsigned generation = 0;  // current generation

  const auto tl = std::chrono::seconds(time_limit);

  do {
    algorithm.evolve();  // evolve the population for one generation

    if ((++generation) % X_INTVL == 0) {
      algorithm.exchangeElite(X_NUMBER);  // exchange top individuals
    }

    ++generation;

    auto now = std::chrono::high_resolution_clock::now();
    if (now - start >= tl) {
      break;
    }

  } while (true);

  best_height = algorithm.getBestFitness();

  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start)
          .count();

  if (debug_sol) {
    const std::string sol_folder = "sol";

    if (!std::filesystem::exists(sol_folder)) {
      std::filesystem::create_directory(sol_folder);
    }

    std::filesystem::path solFilePath = std::filesystem::path(sol_folder) /
                                        relativePath /
                                        (filePath.stem().string() + "Sol.txt");

    if (!std::filesystem::exists(solFilePath.parent_path())) {
      std::filesystem::create_directories(solFilePath.parent_path());
    }

    std::fstream solfile(solFilePath, std::ios::out);

    solfile << max_width << "\n";
    solfile << best_height << "\n";

    std::vector<ranking> sol;
    construct_sol(sol, algorithm.getBestChromosome(), items);
    best_sol = sol;
    unsigned dummy = pack(best_sol, items, max_width, ub, debug_sol, &solfile);
  }

  params << "Tamanho da população: " << p << "\n"
         << "Fração elite: " << pe << "\n"
         << "Fração mutante: " << pm << "\n"
         << "Probabilidade herdar alelo elite: " << rhoe << "\n"
         << "Número de populações independentes: " << K << "\n"
         << "Troca de informação de " << X_NUMBER << " individuos a cada "
         << X_INTVL << " gerações\n"
         << "Número de threads: " << MAXT << "\n"
         << "\n";

  logfile << params.str();

  cout << best_height << "\n";
  cout << "Tempo de execução: " << (double)duration / 1000000.0 << "s\n";

  logfile << "Largura da faixa = " << max_width << "\n";
  logfile << "Melhor altura = " << best_height << "\n\n";
  logfile << "Lower Bound 1 = " << lb1 << "\n";
  logfile << "Lower Bound 2 = " << lb2 << "\n";
  logfile << "Avaliação = " << (double)best_height / lb << "\n\n";

  logfile << "Tempo de execução: " << (double)duration / 1000000.0 << "s"
          << "\n\n"
          << "---------------------------" << "\n\n";

  return 0;
}
