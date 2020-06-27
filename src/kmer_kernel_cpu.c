#include <stdio.h>
// #include <cuda.h>
#include <math.h>
#include "tipos_cpu.h"
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
    int *index = (char *) malloc(sizeof(char) * k);

    for (lint i = rd->start[idx]; i < end; i++)
    {
        index[0] = -1;
        if (i < nN)
        {
            for (lint c = 0; c < k; c++)
            {
                int nuc = rd->data[c + i];

                if (nuc != -1) //Verifica se ha alguem que nao deve ser processado
                {
                    //printf("%d", nuc);
                    index[c] = nuc;
                }
                else
                {
                    index[0] = -1;
                    break;
                }
            }
            if (index[0] != -1)
            {
                for (int j = 0; j < (rd->length[idx] - k + 1); j++)
                {
                    if (rd->counter[idx].index[j][0] == -1)
                    {
                        for (int w = 0; w < k; w++)
                        {
                            rd->counter[idx].index[j][w] = index[w];
                        }
                        rd->counter[idx].frequency[j]++;
                        break;
                    }
                    else
                    {
                        int is_equal = 1;
                        for (int w = 0; w < k; w++)
                        {
                            if (rd->counter[idx].index[j][w] != index[w])
                            {
                                is_equal = 0;
                            }
                        }
                        if (is_equal == 1)
                        {
                            rd->counter[idx].frequency[j]++;
                            break;
                        }
                    }
                }
            }
        }
    }
}

