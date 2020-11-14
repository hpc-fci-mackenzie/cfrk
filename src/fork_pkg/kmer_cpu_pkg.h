#ifndef _kmer_cuh
#define _kmer_cuh

#include "tipos_cpu_pkg.h"

void kmer_main(struct read *rd, lint nN, lint nS, int k);

void ComputeFreq(struct read *rd, int k, int idx, lint nN);

#endif
