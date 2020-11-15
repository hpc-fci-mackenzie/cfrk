#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tipos_cpu.h"
#include "kmer_cpu.h"

void kmer_main(struct read *rd, lint nN, lint nS, int k)
{
    #pragma omp parallel for
    for (int i = 0; i < nS; i++)
    {
        ComputeFreq(rd, k, i, nN);

    }
}