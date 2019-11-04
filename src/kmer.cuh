#ifndef _kmer_cuh
#define _kmer_cuh

#include "tipos.h"

void kmer_main(struct read *rd, lint nN, lint nS, int k, ushort device, cudaStream_t *stream);

__global__ void SetMatrix(float *Mat, ushort offset, int val, int nF);

__global__ void ComputeIndex(char *Seq, int *Index, const int k, lint nN, ushort offset);

__global__ void ComputeFreq(int *Index, int *Freq, lint *start, int *length, ushort offset, int fourk, lint nS, lint nN);

__global__ void ComputeFreqNew(int *Index, int *Freq, lint *start, int *length, ushort offset, int fourk, lint nS);

#endif
