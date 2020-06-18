#ifndef _kmer_cuh
#define _kmer_cuh

#include "tipos_cpu.h"

void kmer_main(struct read *rd, lint nN, lint nS, int k, ushort device);

#ifndef CPU
__global__ void SetMatrix(int *Mat, ushort offset, int val, int nF);
__global__ void ComputeIndex(char *Seq, int *Index, const int k, lint nN, ushort offset);
__global__ void ComputeFreq(int *Index, int *Freq, lint *start, int *length, ushort offset, int fourk, lint nS, lint nN);
__global__ void ComputeFreqNew(int *Index, int *Freq, lint *start, int *length, ushort offset, int fourk, lint nS);
#else
void ComputeFreq(struct read *rd, int k, int idx, lint nN);
#endif

#endif
