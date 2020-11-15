#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tipos_cpu_pkg.h"
#include "kmer_cpu_pkg.h"

void kmer_main(struct read *rd, lint nN, lint nS, int k)
{
    for (int i = 0; i < nS; i++)
    {
        ComputeFreq(rd, k, i, nN);

    }
}