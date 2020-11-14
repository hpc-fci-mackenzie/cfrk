#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tipos_cpu_pkg.h"
#include "kmer_cpu_pkg.h"

void kmer_main(struct read *rd, lint nN, lint nS, int k)
{
    rd->counter = (struct counter *) malloc(sizeof(struct counter) * nS);
    for (int i = 0; i < nS; i++)
    {
        int n_combination = rd->length[i] - k + 1;
        rd->counter[i].kmer = (lint *) malloc(sizeof(lint) * n_combination);
        rd->counter[i].frequency = (char *) malloc(sizeof(char) * n_combination);
        for (int j = 0; j < n_combination; j++)
        {
            rd->counter[i].kmer[j] = -1;
            rd->counter[i].frequency[j] = 0;
        }
    }
    for (int i = 0; i < nS; i++)
    {
        ComputeFreq(rd, k, i, nN);

    }
}