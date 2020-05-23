#ifndef _kmer_cuh
#define _kmer_cuh

#include "tipos_data_struct.h"

void kmer_main(struct chunk *rd, lint nN, lint nS, int k, ushort device);

__global__ void SetMatrix(struct counter *Mat, ushort offset, int nF, int id_sequence);

__global__ void ComputeFrequency(char *Seq, struct counter *d_counter, lint *d_start, int *d_length, const int k, lint nN, ushort offset, lint n_sequence, int n_combination);


#endif
