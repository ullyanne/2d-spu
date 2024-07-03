/*
 * SampleDecoder.cpp
 *
 *  Created on: Jan 14, 2011
 *      Author: rtoso
 */

#include "SampleDecoder.h"
#include <math.h>
#include <iostream>
#include <fstream>

SampleDecoder::~SampleDecoder() {}

// Runs in \Theta(n \log n):
double SampleDecoder::decode(const std::vector<double> &chromosome) const
{
    std::vector<std::pair<double, unsigned>> ranking(chromosome.size());
    return 0;
    for (unsigned i = 0; i < chromosome.size(); ++i)
    {
        ranking[i] = std::pair<double, unsigned>(chromosome[i], i);
    }

    // Here we sort 'permutation', which will then produce a permutation of [n] in pair::second:
    std::sort(ranking.begin(), ranking.end());


    std::pair<int, int> piece;
    std::vector<std::pair<int, int>> spaces;
    unsigned index;

    spaces.push_back({0, 0});

    for(unsigned i = 0; i < chromosome.size(); i++){
        index = ranking[i].second;
        piece = pieces[index];

    }


    std::cout << pieces[0].first << " " << pieces[0].second << std::endl;
    /* Primeira peça a ser inserida*/
    // piece = pieces[ranking[pieceIndex].second];


    return 0;
    // return strip_height;
}