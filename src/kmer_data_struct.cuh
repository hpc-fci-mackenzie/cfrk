#ifndef _kmer_cuh
#define _kmer_cuh

#include "tipos_data_struct.h"

void kmer_main(struct chunk *rd, lint nN, lint nS, int k, ushort device);

__global__ void ComputeFrequence(char *Seq, struct counter *d_counter, lint *d_start, int *d_length, const int k, lint nN, ushort offset, int n_sequence, int n_combination);


#endif
