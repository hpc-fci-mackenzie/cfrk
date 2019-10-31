#ifndef _kmer_cuh
#define _kmer_cuh

#include "tipos.h"

void kmer_main(struct read *rd, lint nN, lint nS, int k, ushort device);

#ifndef CPU
__global__ void SetMatrix(int *Mat, ushort offset, int val, int nF);
__global__ void ComputeIndex(char *Seq, int *Index, const int k, lint nN, ushort offset);
__global__ void ComputeFreq(int *Index, int *Freq, lint *start, int *length, ushort offset, int fourk, lint nS, lint nN);
__global__ void ComputeFreqNew(int *Index, int *Freq, lint *start, int *length, ushort offset, int fourk, lint nS);
#else

/*
	The execution configuration is specified by inserting an expression of the form

		<<< Dg, Db, Ns, S >>> between the function name and the parenthesized argument list, where:

	- Dg is of type dim3 (see Section B.3.2) and specifies the dimension and size of the grid, such that Dg.x * Dg.y * Dg.z equals the number of blocks being launched; Dg.z must be equal to 1 for devices of compute capability 1.x;
	- Db is of type dim3 (see Section B.3.2) and specifies the dimension and size of each block, such that Db.x * Db.y * Db.z equals the number of threads per block;
	- Ns is of type size_t and specifies the number of bytes in shared memory that is dynamically allocated per block for this call in addition to the statically allocated memory; this dynamically allocated memory is used by any of the variables declared as an external array as mentioned in Section B.2.3; Ns is an optional argument which defaults to 0;
	- S is of type cudaStream_t and specifies the associated stream; S is an optional argument which defaults to 0.
*/

void SetMatrix(int *Mat, ushort offset, int val, int nF);
void ComputeIndex(char *Seq, int *Index, const int k, lint nN, ushort offset);
void ComputeFreq(int *Index, int *Freq, lint *start, int *length, ushort offset, int fourk, lint nS, lint nN);
void ComputeFreqNew(int *Index, int *Freq, lint *start, int *length, ushort offset, int fourk, lint nS);
#endif

#endif
