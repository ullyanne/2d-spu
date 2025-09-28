#ifndef RANDOM_MANAGER_H
#define RANDOM_MANAGER_H

#include "MTRand.h"

extern unsigned int global_seed;

void set_global_seed(unsigned int seed);
unsigned int get_seed(unsigned i);

#endif