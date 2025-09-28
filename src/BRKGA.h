/*
 * BRKGA.h
 *
 * This class encapsulates a Biased Random-key Genetic Algorithm (for
 * minimization problems) with K independent Populations stored in two vectors
 * of Population, current and previous. It supports multi-threading via OpenMP,
 * and implements the following key methods:
 *
 * - BRKGA() constructor: initializes the populations with parameters
 * described below.
 * - evolve() operator: evolve each Population following the BRKGA
 * methodology. This method supports OpenMP to evolve up to K independent
 * Populations in parallel. Please note that double Decoder::decode(...) MUST be
 * thread-safe.
 *
 * Required hyperparameters:
 * - n: number of genes in each chromosome
 * - p: number of elements in each population
 * - pe: pct of elite items into each population
 * - pm: pct of mutants introduced at each generation into the population
 * - rhoe: probability that an offspring inherits the allele of its elite parent
 *
 * Optional parameters:
 * - K: number of independent Populations
 * - MAX_THREADS: number of threads to perform parallel decoding -- WARNING:
 * Decoder::decode() MUST be thread-safe!
 *
 * Required templates are:
 * RNG: random number generator that implements the methods below.
 *     - RNG(unsigned long seed) to initialize a new RNG with 'seed'
 *     - double rand() to return a double precision random deviate in range
 * [0,1)
 *     - unsigned long randInt() to return a >=32-bit unsigned random deviate in
 * range [0,2^32-1)
 *     - unsigned long randInt(N) to return a unsigned random deviate in range
 * [0, N] with N < 2^32
 *
 * Decoder: problem-specific decoder that implements any of the decode methods
 * outlined below. When compiling and linking BRKGA with -fopenmp (i.e., with
 * multithreading support via OpenMP), the method must be thread-safe.
 *     - double decode(const vector< double >& chromosome) const, if you don't
 * want to change chromosomes inside the framework, or
 *     - double decode(vector< double >& chromosome) const, if you'd like to
 * update a chromosome
 *
 *  Created on : Jun 22, 2010 by rtoso
 *  Last update: Sep 28, 2010 by rtoso
 *      Authors: Rodrigo Franco Toso <rtoso@cs.rutgers.edu>
 */

#ifndef BRKGA_H
#define BRKGA_H

#include <omp.h>

#include <algorithm>
// #include <deque>
#include <exception>
#include <limits>
#include <stdexcept>
#include <typeinfo>

#include "Packing.h"
#include "Population.h"
#include "random_manager.h"

template <class Decoder, class RNG>
class BRKGA {
 public:
  /*
   * Default constructor
   * Required hyperparameters:
   * - n: number of genes in each chromosome
   * - p: number of elements in each population
   * - pe: pct of elite items into each population
   * - pm: pct of mutants introduced at each generation into the population
   * - rhoe: probability that an offspring inherits the allele of its elite
   * parent
   *
   * Optional parameters:
   * - K: number of independent Populations
   * - MAX_THREADS: number of threads to perform parallel decoding
   *                WARNING: Decoder::decode() MUST be thread-safe; safe if
   * implemented as
   *                + double Decoder::decode(std::vector< double >& chromosome)
   * const
   */
  BRKGA(unsigned n, unsigned p, double pe, double pm, double rhoe,
        const Decoder& refDecoder, RNG& refRNG, RNG& r_ls, RNG& r_init,
        RNG& r_elite, unsigned max_ls_attempts, double i1, double i2, double i3,
        double window_init, double window_ls, double top_elite, unsigned K = 1,
        unsigned MAX_THREADS = 1);

  /**
   * Destructor
   */
  ~BRKGA();

  /**
   * Resets all populations with brand new keys
   */
  void reset();

  void partial_reset();

  /**
   * Evolve the current populations following the guidelines of BRKGAs
   * @param generations number of generations (must be even and nonzero)
   * @param J interval to exchange elite chromosomes (must be even; 0 ==> no
   * synchronization)
   * @param M number of elite chromosomes to select from each population in
   * order to exchange
   */
  void evolve(unsigned generations = 1);

  /**
   * Exchange elite-solutions between the populations
   * @param M number of elite chromosomes to select from each population
   */
  void exchangeElite(unsigned M);

  /**
   * Returns the current population
   */
  const Population& getPopulation(unsigned k = 0) const;

  /**
   * Returns the chromosome with best fitness so far among all populations
   */
  const std::vector<double>& getBestChromosome() const;

  /**
   * Returns the best fitness found so far among all populations
   */
  double getBestFitness() const;

  // Return copies to the internal parameters:
  unsigned getN() const;
  unsigned getP() const;
  unsigned getPe() const;
  unsigned getPm() const;
  unsigned getPo() const;
  double getRhoe() const;
  unsigned getK() const;
  unsigned getMAX_THREADS() const;

 private:
  // Hyperparameters:
  const unsigned n;   // number of genes in the chromosome
  const unsigned p;   // number of elements in the population
  const unsigned pe;  // number of elite items in the population
  const unsigned pm;  // number of mutants introduced at each generation into
                      // the population
  const double rhoe;  // probability that an offspring inherits the allele of
                      // its elite parent

  // Templates:
  RNG& refRNG;  // reference to the random number generator
  RNG& r_ls;
  RNG& r_init;
  RNG& r_elite;

  double i1;
  double i2;
  double i3;
  double window_init;
  double window_ls;
  double top_elite;

  unsigned max_ls_attempts;
  const Decoder& refDecoder;  // reference to the problem-dependent Decoder

  // Parallel populations parameters:
  const unsigned K;            // number of independent parallel populations
  const unsigned MAX_THREADS;  // number of threads for parallel decoding

  // Data:
  std::vector<Population*> previous;  // previous populations
  std::vector<Population*> current;   // current populations

  // Local operations:
  void initialize(const unsigned i,
                  const unsigned start = 0);  // initialize current population
                                              // 'i' with random keys
  void evolution(Population& curr, Population& next);
  bool isRepeated(const std::vector<double>& chrA,
                  const std::vector<double>& chrB) const;
};

template <class Decoder, class RNG>
BRKGA<Decoder, RNG>::BRKGA(unsigned _n, unsigned _p, double _pe, double _pm,
                           double _rhoe, const Decoder& decoder, RNG& rng,
                           RNG& r_ls, RNG& r_init, RNG& r_elite,
                           unsigned max_ls_attempts, double i1, double i2,
                           double i3, double window_init, double window_ls,
                           double top_elite, unsigned _K, unsigned MAX)
    : n(_n),
      p(_p),
      pe(unsigned(_pe * p)),
      pm(unsigned(_pm * p)),
      rhoe(_rhoe),
      refRNG(rng),
      r_ls(r_ls),
      r_init(r_init),
      r_elite(r_elite),
      max_ls_attempts(max_ls_attempts),
      refDecoder(decoder),
      i1(i1),
      i2(i2),
      i3(i3),
      window_init(window_init),
      window_ls(window_ls),
      top_elite(top_elite),
      K(_K),
      MAX_THREADS(MAX),
      previous(K, 0),
      current(K, 0)
{
  // Error check:
  using std::range_error;
  if (n == 0) {
    throw range_error("Chromosome size equals zero.");
  }
  if (p == 0) {
    throw range_error("Population size equals zero.");
  }
  if (pe == 0) {
    throw range_error("Elite-set size equals zero.");
  }
  if (pe > p) {
    throw range_error("Elite-set size greater than population size (pe > p).");
  }
  if (pm > p) {
    throw range_error("Mutant-set size (pm) greater than population size (p).");
  }
  if (pe + pm > p) {
    throw range_error("elite + mutant sets greater than population size (p).");
  }
  if (K == 0) {
    throw range_error("Number of parallel populations cannot be zero.");
  }

  // Initialize and decode each chromosome of the current population, then copy
  // to previous:
  for (unsigned i = 0; i < K; ++i) {
    // Allocate:
    current[i] = new Population(n, p);

    // Initialize:
    initialize(i);

    // Then just copy to previous:
    previous[i] = new Population(*current[i]);
  }
}

template <class Decoder, class RNG>
BRKGA<Decoder, RNG>::~BRKGA()
{
  for (unsigned i = 0; i < K; ++i) {
    delete current[i];
    delete previous[i];
  }
}

template <class Decoder, class RNG>
const Population& BRKGA<Decoder, RNG>::getPopulation(unsigned k) const
{
  return (*current[k]);
}

template <class Decoder, class RNG>
double BRKGA<Decoder, RNG>::getBestFitness() const
{
  double best = current[0]->chromosome_packing_info[0].fitness;
  for (unsigned i = 1; i < K; ++i) {
    if (current[i]->chromosome_packing_info[0].fitness < best) {
      best = current[i]->chromosome_packing_info[0].fitness;
    }
  }
  return best;
}

template <class Decoder, class RNG>
const std::vector<double>& BRKGA<Decoder, RNG>::getBestChromosome() const
{
  unsigned bestK = 0;
  for (unsigned i = 1; i < K; ++i) {
    if (current[i]->getBestFitness() < current[bestK]->getBestFitness()) {
      bestK = i;
    }
  }

  return current[bestK]->getChromosome(0);  // The top one :-)
}

template <class Decoder, class RNG>
void BRKGA<Decoder, RNG>::reset()
{
  for (unsigned i = 0; i < K; ++i) {
    initialize(i);
  }
}

template <class Decoder, class RNG>
void BRKGA<Decoder, RNG>::partial_reset()
{
  for (unsigned i = 0; i < K; ++i) {
    initialize(i, pe);
  }
}

template <class Decoder, class RNG>
void BRKGA<Decoder, RNG>::evolve(unsigned generations)
{
  if (generations == 0) {
    throw std::range_error("Cannot evolve for 0 generations.");
  }

  for (unsigned i = 0; i < generations; ++i) {
    for (unsigned j = 0; j < K; ++j) {
      evolution(*current[j],
                *previous[j]);  // First evolve the population (curr, next)
      std::swap(current[j],
                previous[j]);  // Update (prev = curr; curr = prev == next)
    }
  }
}

template <class Decoder, class RNG>
void BRKGA<Decoder, RNG>::exchangeElite(unsigned M)
{
  if (M == 0 || M >= p) {
    throw std::range_error("M cannot be zero or >= p.");
  }

  for (unsigned i = 0; i < K; ++i) {
    // Population i will receive some elite members from each Population j
    // below:
    unsigned dest = p - 1;  // Last chromosome of i (will be updated below)
    for (unsigned j = 0; j < K; ++j) {
      if (j == i) {
        continue;
      }

      // Copy the M best of Population j into Population i:
      for (unsigned m = 0; m < M; ++m) {
        // Copy the m-th best of Population j into the 'dest'-th position of
        // Population i:
        const std::vector<double>& bestOfJ = current[j]->getChromosome(m);

        std::copy(bestOfJ.begin(), bestOfJ.end(),
                  current[i]->getChromosome(dest).begin());

        current[i]->chromosome_packing_info[dest].fitness =
            current[j]->chromosome_packing_info[m].fitness;

        --dest;
      }
    }
  }

  for (int j = 0; j < int(K); ++j) {
    current[j]->sortFitness();
  }
}

template <class Decoder, class RNG>
inline void BRKGA<Decoder, RNG>::initialize(const unsigned i,
                                            const unsigned start)
{
  std::vector<std::vector<double>> chromosome_groups(3, std::vector<double>(n));

  int lims[3];
  lims[0] = i1 * p - 1;
  lims[1] = lims[0] + (i2 * p) - 1;
  lims[2] = lims[1] + (i3 * p) - 1;

  for (unsigned j = 0; j < 3; j++) {
    for (unsigned k = 0; k < n; ++k) {
      chromosome_groups[j][k] = refDecoder.rank_groups[j][k].chromosome;
    }
  }

  (*current[i])(0) = chromosome_groups[0];
  (*current[i])(1) = chromosome_groups[1];
  (*current[i])(2) = chromosome_groups[2];

  int window = window_init * n;
  if (window == 0) {
    window = 1;
  }

  for (unsigned j = start + 3; j < p; ++j) {
    int first_piece = r_init.randInt(n - 1);

    int s_min = std::max(0, first_piece - window);
    int s_max = std::min(int(n - 1), first_piece + window);

    int second_piece = s_min + r_init.randInt(s_max - s_min);

    while (second_piece == first_piece) {
      second_piece = s_min + r_init.randInt(s_max - s_min);
    }

    if (j < lims[0]) {
      std::swap(chromosome_groups[0][first_piece],
                chromosome_groups[0][second_piece]);
      (*current[i])(j) = chromosome_groups[0];
    }
    else if (j >= lims[0] && j < lims[1]) {
      std::swap(chromosome_groups[1][first_piece],
                chromosome_groups[1][second_piece]);
      (*current[i])(j) = chromosome_groups[1];
    }
    else if (j >= lims[1]) {
      std::swap(chromosome_groups[2][first_piece],
                chromosome_groups[2][second_piece]);
      (*current[i])(j) = chromosome_groups[2];
    }
  }

// Decode:
#ifdef _OPENMP
#pragma omp parallel for num_threads(MAX_THREADS)
#endif
  for (int j = start; j < int(p); ++j) {
    double fitness = refDecoder.decode((*current[i])(j));

    current[i]->setFitness(j, fitness);
  }

  // Sort:

  current[i]->sortFitness();
}

template <class Decoder, class RNG>
inline void BRKGA<Decoder, RNG>::evolution(Population& curr, Population& next)
{
  // We now will set every chromosome of 'current', iterating with 'i':
  unsigned i = 0;  // Iterate chromosome by chromosome
  unsigned j = 0;  // Iterate allele by allele

  // 2. The 'pe' best chromosomes are maintained, so we just copy these into
  // 'current':

  for (unsigned k = 0; k < pe; ++k) {
    for (j = 0; j < n; ++j) {
      next(k, j) = curr(curr.chromosome_packing_info[k].chromosome, j);
    }

    unsigned fitness = curr.chromosome_packing_info[k].fitness;

    double rand_val = r_elite.randDblExc();
    if (rand_val < top_elite) {
      std::vector<double>& chr = next(k);

      std::vector<ranking> rank(chr.size());
      for (unsigned i = 0; i < chr.size(); i++) {
        rank[i].chromosome = chr[i];
        rank[i].index = i;
      }
      std::sort(rank.begin(), rank.end(), sort_rank);

      unsigned new_fitness =
          refDecoder.local_search(rank, curr.chromosome_packing_info[k].fitness,
                                  r_ls, max_ls_attempts, window_ls);

      if (new_fitness < fitness) {
        fitness = new_fitness;

        std::vector<ranking> new_sol(chr.size());
        encode(new_sol, rank, chr.size());
        std::sort(new_sol.begin(), new_sol.end(),
                  [](const ranking& a, const ranking& b) {
                    return a.index < b.index;
                  });

        for (unsigned i = 0; i < new_sol.size(); i++) {
          chr[i] = new_sol[i].chromosome;
        }
      }
    }

    next.chromosome_packing_info[k].fitness = fitness;
    next.chromosome_packing_info[k].chromosome = k;
  }

  i = pe;
  // 3. We'll mate 'p - pe - pm' pairs; initially, i = pe, so we need to iterate
  // until i < p - pm:
  while (i < p - pm) {
    // Select an elite parent:
    const unsigned eliteParent = (refRNG.randInt(pe - 1));

    // Select a non-elite parent:
    const unsigned noneliteParent = pe + (refRNG.randInt(p - pe - 1));

    // Mate:
    for (j = 0; j < n; ++j) {
      const unsigned sourceParent =
          ((refRNG.rand() < rhoe) ? eliteParent : noneliteParent);

      next(i, j) =
          curr(curr.chromosome_packing_info[sourceParent].chromosome, j);
    }

    ++i;
  }

  // We'll introduce 'pm' mutants:
  while (i < p) {
    for (j = 0; j < n; ++j) {
      next(i, j) = refRNG.rand();
    }
    ++i;
  }

// Time to compute fitness, in parallel:
#ifdef _OPENMP
#pragma omp parallel for num_threads(MAX_THREADS)
#endif
  for (int i = int(pe); i < int(p); ++i) {
    next.setFitness(i, refDecoder.decode(next.population[i]));
  }

  // Now we must sort 'current' by fitness, since things might have changed:
  next.sortFitness();
}

template <class Decoder, class RNG>
unsigned BRKGA<Decoder, RNG>::getN() const
{
  return n;
}

template <class Decoder, class RNG>
unsigned BRKGA<Decoder, RNG>::getP() const
{
  return p;
}

template <class Decoder, class RNG>
unsigned BRKGA<Decoder, RNG>::getPe() const
{
  return pe;
}

template <class Decoder, class RNG>
unsigned BRKGA<Decoder, RNG>::getPm() const
{
  return pm;
}

template <class Decoder, class RNG>
unsigned BRKGA<Decoder, RNG>::getPo() const
{
  return p - pe - pm;
}

template <class Decoder, class RNG>
double BRKGA<Decoder, RNG>::getRhoe() const
{
  return rhoe;
}

template <class Decoder, class RNG>
unsigned BRKGA<Decoder, RNG>::getK() const
{
  return K;
}

template <class Decoder, class RNG>
unsigned BRKGA<Decoder, RNG>::getMAX_THREADS() const
{
  return MAX_THREADS;
}

#endif
