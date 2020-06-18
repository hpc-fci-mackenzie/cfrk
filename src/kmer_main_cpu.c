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
//   int *d_Index; // Index vector
//   char *d_Seq; // Seq matrix
//   int *d_Freq;// Frequence vector
//   int fourk; // 4 power k
//   lint *d_start;
//   int *d_length; // The beggining and the length of each sequence
//   size_t size[5];

   // printf("@deb | kmer_main | nN: %ld | nS: %ld | k: %d | device: %d\n", nN, nS, k, device);
//   fourk = POW(k);
   // printf("> fourk: %d\n", fourk);
   
   //---------------------------------------------------------------------------
//   size[0] = nN * sizeof(char); // d_Seq and Seq size
//   size[1] = nN * sizeof(int); // d_Index and Index size
//   size[2] = nS * sizeof(int); // d_length
//   size[3] = nS * fourk * sizeof(int); // Freq and d_Freq
//   size[4] = nS * sizeof(lint); // d_start
   // printf("> size[0]: %ld\n", size[0]);
   // printf("> size[1]: %ld\n", size[1]);
   // printf("> size[2]: %ld\n", size[2]);
   // printf("> size[3]: %ld\n", size[3]);
   // printf("> size[4]: %ld\n", size[4]);

   //---------------------------------------------------------------------------

   rd->counter = (struct counter *) malloc(sizeof(struct counter)*nS);
//    #pragma omp parallel for
   for (int i = 0; i < nS; i++)
   {
       int n_combination = rd->length[i] - k + 1;
       rd->counter[i].index = (int **) malloc(sizeof(int *)*n_combination);
       rd->counter[i].frequency = (int *) malloc(sizeof(int)*n_combination);

       for (int j = 0; j < n_combination; j++)
       {
           rd->counter[i].index[j] = (int *) malloc(sizeof(int)*k);
           rd->counter[i].index[j][0] = -1;
           rd->counter[i].frequency[j] = 0;
       }
   }

   //************************************************

    for (int i = 0; i < nS; i++) {
        ComputeFreq(rd, k, i, nN);
    }

   //---------------------------------------------------------------------------
//    for (int i = 0; i < nS; i++)
//    {
//        int n_combination = rd->length[i] - k + 1;
//
//        for (int j = 0; j < n_combination; j++)
//        {
//            if (rd->counter[i].frequency[j] > 0)
//                printf("Sequence: %d Index: %d %d %d %d Frequency: %d\n",
//                       i,
//                       (int) rd->counter[i].index[j][0],
//                        (int) rd->counter[i].index[j][1],
//                        (int) rd->counter[i].index[j][2],
//                        (int) rd->counter[i].index[j][3],
//                        rd->counter[i].frequency[j]);
//        }
//    }
}