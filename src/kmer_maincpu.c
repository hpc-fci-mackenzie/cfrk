#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tipos.h"
#include "kmer.cuh"

/*
#include <cuda.h>
#include <math.h>
#include <stdint.h>
*/

void kmer_main(struct read *rd, lint nN, lint nS, int k, ushort device)
{
   int *d_Index; // Index vector
   char *d_Seq; // Seq matrix
   int *d_Freq;// Frequence vector
   int fourk; // 4 power k
   lint *d_start;
   int *d_length; // The beggining and the length of each sequence
   size_t size[5];

   // printf("@deb | kmer_main | nN: %ld | nS: %ld | k: %d | device: %d\n", nN, nS, k, device);
   fourk = POW(k);
   // printf("> fourk: %d\n", fourk);
   
   //---------------------------------------------------------------------------
   size[0] = nN * sizeof(char); // d_Seq and Seq size
   size[1] = nN * sizeof(int); // d_Index and Index size
   size[2] = nS * sizeof(int); // d_length
   size[3] = nS * fourk * sizeof(int); // Freq and d_Freq
   size[4] = nS * sizeof(lint); // d_start
   // printf("> size[0]: %ld\n", size[0]);
   // printf("> size[1]: %ld\n", size[1]);
   // printf("> size[2]: %ld\n", size[2]);
   // printf("> size[3]: %ld\n", size[3]);
   // printf("> size[4]: %ld\n", size[4]);

   //---------------------------------------------------------------------------
   d_Seq    = (char*)calloc(size[0], sizeof(char));
   d_Index  = ( int*)calloc(size[1], sizeof(int));
   d_start  = (lint*)calloc(size[4], sizeof(lint));
   d_length = ( int*)calloc(size[2], sizeof(int));
   d_Freq   = ( int*)calloc(size[3], sizeof(int));

   // Copies data between host and device.
   memcpy(d_Seq, rd->data, size[0]);
   memcpy(d_start, rd->start, size[4]);
   memcpy(d_length, rd->length, size[2]);

   //************************************************
   ComputeIndex(d_Seq, d_Index, k, nN, size[1]);
   ComputeFreq(d_Index, d_Freq, d_start, d_length, size[4], fourk, nS, nN);

   rd->Freq = (int*)malloc(size[3]);
   memcpy(rd->Freq, d_Freq, size[3]);
   //************************************************
   free(d_Seq);
   free(d_Freq);
   free(d_Index);
   free(d_start);
   free(d_length);
   //---------------------------------------------------------------------------
}