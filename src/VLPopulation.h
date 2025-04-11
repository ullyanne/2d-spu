/*
 * Population.h
 *
 * Encapsulates a population of chromosomes represented by a vector of doubles.
 * We don't decode nor deal with random numbers here; instead, we provide
 * private support methods to set the fitness of a specific chromosome as well
 * as access methods to each allele. Note that the BRKGA class must have access
 * to such methods and thus is a friend.
 *
 *  Created on : Jun 21, 2010 by rtoso
 *  Last update: Nov 15, 2010 by rtoso
 *      Authors: Rodrigo Franco Toso <rtoso@cs.rutgers.edu>
 */

#ifndef VLPOPULATION_H
#define VLPOPULATION_H

#include <algorithm>
#include <exception>
#include <stdexcept>
#include <vector>

typedef struct {
  double fitness;
  unsigned chromosome;
} vl_packing_info;

class VLPopulation {
  template <class Decoder, class RNG>
  friend class VLBRKGA;

 public:
  unsigned getN() const;  // Size of each chromosome
  unsigned getP() const;  // Size of population

  // double operator()(unsigned i, unsigned j) const;	// Direct access to
  // allele j of chromosome i

  // These methods REQUIRE fitness to be sorted, and thus a call to
  // sortFitness() beforehand (this is done by BRKGA, so rest assured:
  // everything will work just fine with BRKGA).
  double getBestFitness() const;  // Returns the best fitness in this population
  double getFitness(unsigned i) const;  // Returns the fitness of chromosome i
  const std::vector<double>& getChromosome(
      unsigned i) const;  // Returns i-th best chromosome

 private:
  VLPopulation(const VLPopulation& other);
  VLPopulation(unsigned n, unsigned p);
  ~VLPopulation();

  std::vector<std::vector<double>>
      population;  // Population as vectors of prob.
  std ::vector<vl_packing_info> chromosome_packing_info;

  void sortFitness();  // Sorts 'fitness' by its first parameter
  void setFitness(unsigned i, double f);  // Sets the fitness of chromosome i
  std::vector<double>& getChromosome(unsigned i);  // Returns a chromosome

  double& operator()(unsigned i,
                     unsigned j);  // Direct access to allele j of chromosome i
  std::vector<double>& operator()(unsigned i);  // Direct access to chromosome i
};

VLPopulation::VLPopulation(const VLPopulation& pop)
    : population(pop.population),
      chromosome_packing_info(pop.chromosome_packing_info)
{
}

VLPopulation::VLPopulation(const unsigned n, const unsigned p)
    : population(p, std::vector<double>(n, 0.0)), chromosome_packing_info(p)
{
  if (p == 0) {
    throw std::range_error("Population size p cannot be zero.");
  }
  if (n == 0) {
    throw std::range_error("Chromosome size n cannot be zero.");
  }
}

VLPopulation::~VLPopulation() {}

unsigned VLPopulation::getN() const { return population[0].size(); }

unsigned VLPopulation::getP() const { return population.size(); }

double VLPopulation::getBestFitness() const { return getFitness(0); }

double VLPopulation::getFitness(unsigned i) const
{
  return chromosome_packing_info[i].fitness;
}

const std::vector<double>& VLPopulation::getChromosome(unsigned i) const
{
  return population[chromosome_packing_info[i].chromosome];
}

std::vector<double>& VLPopulation::getChromosome(unsigned i)
{
  return population[chromosome_packing_info[i].chromosome];
}

void VLPopulation::setFitness(unsigned i, double f)
{
  chromosome_packing_info[i].fitness = f;
  chromosome_packing_info[i].chromosome = i;
}

void VLPopulation::sortFitness()
{
  std::sort(chromosome_packing_info.begin(), chromosome_packing_info.end(),
            [](const vl_packing_info& a, const vl_packing_info& b) {
              return a.fitness < b.fitness;
            });
}

// double Population::operator()(unsigned chromosome, unsigned allele) const {
//	return population[chromosome][allele];
// }

double& VLPopulation::operator()(unsigned chromosome, unsigned allele)
{
  return population[chromosome][allele];
}

std::vector<double>& VLPopulation::operator()(unsigned chromosome)
{
  return population[chromosome];
}

#endif
