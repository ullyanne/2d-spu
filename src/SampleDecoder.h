/*
 * SampleDecoder.h
 *
 * Any decoder must have the format below, i.e., implement the method
 * decode(std::vector< double >&) returning a double corresponding to the
 * fitness of that vector. If parallel decoding is to be used in the BRKGA
 * framework, then the decode() method _must_ be thread-safe; the best way to
 * guarantee this is by adding 'const' to the end of decode() so that the
 * property will be checked at compile time.
 *
 * The chromosome inside the BRKGA framework can be changed if desired. To do
 * so, just use the first signature of decode() which allows for modification.
 * Please use double values in the interval [0,1) when updating, thus obeying
 * the BRKGA guidelines.
 *
 *  Created on: Jan 14, 2011
 *      Author: rtoso
 */

#ifndef SAMPLEDECODER_H
#define SAMPLEDECODER_H

#include <algorithm>
#include <limits>
#include <list>
#include <vector>

#include "Item.h"
#include "Packing.h"

using namespace std;

class SampleDecoder {
 public:
  SampleDecoder(const vector<ranking> lns_seq, const vector<item> items,
                const int max_width, const int ub)
      : lns_seq(lns_seq), items(items), max_width(max_width), ub(ub)
  {
  }
  ~SampleDecoder();
  double decode(
      const vector<double> &chromosome,
      std::unordered_map<unsigned, std::vector<unsigned>> &clients_to_layers,
      std::unordered_map<unsigned, unsigned> &layers_to_index) const;

 private:
  const vector<ranking> lns_seq;
  const vector<item> items;
  const int max_width;
  const int ub;
};

#endif