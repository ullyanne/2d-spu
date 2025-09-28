#include "random_manager.h"

unsigned int global_seed;

void set_global_seed(unsigned int seed) { global_seed = seed; }

unsigned int get_seed(unsigned i) { return global_seed + i; }