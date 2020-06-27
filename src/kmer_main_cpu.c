#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tipos_cpu.h"
#include "kmer_cpu.h"

void kmer_main(struct read *rd, lint nN, lint nS, int k)
{
    rd->counter = (struct counter *) malloc(sizeof(struct counter) * nS);
    #pragma omp parallel for
    for (int i = 0; i < nS; i++)
    {
        int n_combination = rd->length[i] - k + 1;
        rd->counter[i].index = (char **) malloc(sizeof(char *) * n_combination);
        rd->counter[i].frequency = (int *) malloc(sizeof(int) * n_combination);
        #pragma omp parallel for
        for (int j = 0; j < n_combination; j++)
        {
            rd->counter[i].index[j] = (char *) malloc(sizeof(char) * k);
            rd->counter[i].index[j][0] = -1;
            rd->counter[i].frequency[j] = 0;
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < nS; i++)
    {
        ComputeFreq(rd, k, i, nN);
    }
}