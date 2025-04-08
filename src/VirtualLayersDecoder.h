/*
 * VirtualLayersDecoder.h
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

#ifndef VIRTUALLAYERSDECODER_H
#define VIRTUALLAYERSDECODER_H

#include <algorithm>
#include <limits>
#include <list>
#include <vector>

#include "Item.h"
#include "Packing.h"

using namespace std;

class VirtualLayersDecoder {
 public:
  VirtualLayersDecoder(const vector<vector<ranking>> initial_seqs,
                       const vector<item> items, const int max_width,
                       const int ub)
      : initial_seqs(initial_seqs), items(items), max_width(max_width), ub(ub)
  {
  }
  ~VirtualLayersDecoder();
  double decode(const vector<double> &chromosome) const;

  const vector<vector<ranking>> initial_seqs;
  const vector<item> items;
  const int max_width;
  const int ub;

 private:
};

#endif