#include <stdio.h>
// #include <cuda.h>
#include <math.h>
#include "tipos_cpu_pkg.h"
#include <string.h>
#include <stdlib.h>


/*
   Compute k-mer frequency
   Na segunda fase é calculada a repetição dos valores convertidos na fase
   anterior. O kernel ComputeFreq é responsável pela contabilização da
   repetição dos k-mers. O resultado do cálculo é armazenado em um vetor
   denominado vetor de frequência.
*/
void ComputeFreq(struct read *rd, int k, int idx, lint nN)
{
    int end = rd->start[idx] + (rd->length[idx] + 1);


    lint index = -1;
    for (int i = rd->start[idx]; i < end; i++)
    {
        index = 0;
        if (i < nN)
        {
            for (lint c = 0; c < k; c++)
            {
                lint nuc = rd->data[c + i];

                if (nuc != -1) //Verifica se ha alguem que nao deve ser processado
                {
                    index += nuc * pow(4, (k - 1) - c);
                }
                else
                {
                    index = -1;
                    break;
                }
            }
            if (index > -1)
            {
                for (int j = 0; j < (rd->length[idx] - k + 1); j++)
                {
                    if (rd->counter[idx].kmer[j] == -1)
                    {
                        rd->counter[idx].kmer[j] = index;
                        rd->counter[idx].frequency[j]++;
                        break;
                    }
                    else if (rd->counter[idx].kmer[j] == index)
                    {
                        rd->counter[idx].frequency[j]++;
                        break;
                    }
                }
            }
        }
    }
}

