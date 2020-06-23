#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tipos_cpu.h"
#include "kmer_cpu.h"

/*
#include <cuda.h>
#include <math.h>
#include <stdint.h>
*/

void kmer_main(struct read *rd, lint nN, lint nS, int k, ushort device)
{
   rd->counter = (struct counter *) malloc(sizeof(struct counter)*nS);
   for (int i = 0; i < nS; i++)
   {
       int n_combination = rd->length[i] - k + 1;
       rd->counter[i].index = (char **) malloc(sizeof(char *)*n_combination);
       rd->counter[i].frequency = (int *) malloc(sizeof(int)*n_combination);

       for (int j = 0; j < n_combination; j++)
       {
           rd->counter[i].index[j] = (char *) malloc(sizeof(char)*k);
           rd->counter[i].index[j][0] = -1;
           rd->counter[i].frequency[j] = 0;
       }
   }
    for (int i = 0; i < nS; i++) {
        ComputeFreq(rd, k, i, nN);
    }
}