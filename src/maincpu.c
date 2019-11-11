/*
   CFRK - Contabilizador da Frequência de Repetição de K-mer (Versão CPU)
   Instiuição: Universidade Presbiteriana Mackenzie
   
   Este código foi baseado no código original CFRK-MT desenvolvido por
   Fabrício Gomes Vilasboas (LNCC | Laboratório Nacional de Computação Científica)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
// #include <omp.h>
#include "tipos.h"
#include "fastaIO.h"
#include "kmer.cuh"

#define FILENAME_LENGTH 512

//*** Global Variables ***
struct read *chunk;
lint *nS, *nN;
int offset;
int device;
int k;
char file_out[FILENAME_LENGTH];

void PrintFreq(struct seq *seq, struct read *pchunk, int nChunk, int chunkSize)
{
   FILE *out;
   int cont = 0;
   int cont_seq = 0;
   char str[32];
   lint fourk = pow(4, k);

   out = fopen(file_out, "w");

   for (int j = 0; j < nChunk; j++)
   {
      cont = 0;
      cont_seq = 0;
      for (int i = 0; i < (chunkSize*fourk); i++)
      {
         if (i % fourk == 0)
         {
            cont = 0;
            sprintf(str, "\n");
            // printf("%s", str);
            fwrite(str, sizeof(char), sizeof(str), out);
            // printf("%s", seq[cont_seq].header);
            cont_seq++;
         }
         if (pchunk[j].Freq[i] != 0)
         {
            sprintf(str, "%d ", pchunk[j].Freq[i]);
            printf("\tFreq[%d]: %s\n", i, str);
            fwrite(str, sizeof(char), sizeof(str), out);
         }
         cont++;
      }
   }
   fclose(out);
}

struct read* SelectChunkRemain(struct read *rd, ushort chunkSize, ushort it, lint max, lint gnS, lint *nS, lint gnN, lint *nN, int nt)
{
   struct read *chunk = 0;
   lint i;
   lint j;
   lint length = 0;

   // Size to be allocated
   for (i = 0; i < max; i++)
   {
      lint id = chunkSize*it + i;
      if (id > gnS-1)
      {
         break;
      }
      length += rd->length[id]+1;
   }

   printf("@deb | SelectChunkRemain | length: %ld\n", length);

   chunk = (struct read *)malloc(sizeof(struct read));
   chunk->data = (char*)malloc(sizeof(char) * length);
   chunk->length = (int*)malloc(sizeof(int) * chunkSize);
   chunk->start = (lint*)malloc(sizeof(lint) * chunkSize);

   // Copy rd->data to chunk->data
   lint start = rd->start[chunkSize*it];
   lint end = start + (lint)length;

   // Substituir esta parte por pthreads!
   // #pragma omp parallel for num_threads(nt)
   for (j = start; j < end; j++)
   {
      chunk->data[j-start] = rd->data[j];
   }

   chunk->length[0] = rd->length[chunkSize*it];
   chunk->start[0] = 0;

   // Copy start and length
   for (i = 1; i < max; i++)
   {
      lint id = chunkSize*it + i;
      chunk->length[i] = rd->length[id];
      chunk->start[i] = chunk->start[i-1]+(chunk->length[i-1]+1);
   }

   *nN = length;
   *nS = max;

   return chunk;
}

void SelectChunk(struct read *chunk, const int nChunk, struct read *rd, ushort chunkSize, lint max, lint gnS, lint *nS, lint gnN, lint *nN, int nt)
{
   lint i, j, it;

   // printf("@deb | Parameters\n\tnChunk: %d\n\tchunkSize: %d\n\tmax: %ld\n\tgnS: %ld\n\tnS: %ld\n\tgnN: %ld\n\tnN: %ld\n\tnt: %d\n", nChunk, chunkSize, max, gnS, *nS, gnN, *nN, nt);

   for (it = 0; it < nChunk; it++)
   {
      lint length = 0;

      // Size to be allocated
      for (i = 0; i < max; i++)
      {
         lint id = chunkSize*it + i;
         
         if (id > gnS-1)
         {
            break;
         }

         length += rd->length[id]+1;
      }

      // cudaMallocHost((void**)&chunk[it].data, sizeof(char)*length);
      // cudaMallocHost((void**)&chunk[it].length, sizeof(int)*chunkSize);
      // cudaMallocHost((void**)&chunk[it].start, sizeof(lint)*chunkSize);
      chunk[it].data = (char*)calloc(length, sizeof(char));
      chunk[it].length = (int*)calloc(chunkSize, sizeof(int));
      chunk[it].start = (lint*)calloc(chunkSize, sizeof(lint));

      // printf("- @deb | chunk[it].data size: %ld (%ld)\n", length, sizeof(char));
      // printf("- @deb | chunk[it].length size: %d (%ld)\n", chunkSize, sizeof(int));
      // printf("- @deb | chunk[it].start size: %d (%ld)\n", chunkSize, sizeof(lint));

      // Copy rd->data to chunk->data
      lint start = rd->start[chunkSize*it];
      lint end = start + (lint)length;

      // printf("@deb (%ld)| length: %ld - start: %ld - end: %ld\n", it, length, start, end);

      // #pragma omp parallel for num_threads(nt)
      for(j = start; j < end; j++) {
         chunk[it].data[j-start] = rd->data[j];
         // printf("chunk[%ld].data[%ld]: %d\n", it, j-start, chunk[it].data[j-start]);
      }

      chunk[it].length[0] = rd->length[chunkSize*it];
      chunk[it].start[0] = 0;

      // printf("@deb | chunk[%ld].length[0]: %d\n", it, chunk[it].length[0]);

      // Copy start and length
      for (i = 1; i < max; i++)
      {
         lint id = chunkSize*it + i;
         chunk[it].length[i] = rd->length[id];
         chunk[it].start[i] = chunk[it].start[i-1]+(chunk[it].length[i-1]+1);
      }

      nN[it] = length;
      nS[it] = max;
   }
}

void *LaunchKmer(void* threadId)
{
   int start = 0;
   int end = offset;
   // printf("\t>>> ttid: %ld - start: %d - end: %d - offset: %d\n", *tid, start, end, offset);

   for (int i = start; i < end; i++)
   {
      printf("\tkmer_main - chunk(%d)\n", i);
      kmer_main(&chunk[i], nN[i], nS[i], k, device);
      // Waits for stream tasks to complete.
      free(chunk[i].data);
      free(chunk[i].length);
      free(chunk[i].start);
   }

   return NULL;
}

int main(int argc, char* argv[])
{
   lint gnN = 0, gnS = 0, chunkSize = 8192;
   int devCount = 1;
   int nt = 12;
   char dataset[FILENAME_LENGTH] = {0};

   /*
      (!) CONFIGURAÇÃO PADRÃO
         dataset: teste.fasta
         outfile: resultado.cfrk
         número de threads: 12
         tamanho do chunksize: 8192
         k: 2
         nt: 12
   */

   memset(file_out, 0, FILENAME_LENGTH);

   if(argc > 1) {
      if(argv[1]) {
         strcpy(dataset, argv[1]);
      }

      if(argv[2]) {
         strcpy(file_out, argv[2]);
      }

      if(argv[3]) {
         k = atoi(argv[3]);
      }

      if(argv[4]) {
         nt = atoi(argv[4]);
      }

      if(argv[5]) {
         chunkSize = atoi(argv[5]);
      }
   } else {
      strcpy(dataset, "teste.fasta");
      strcpy(file_out, "resultado.cfrk");
      k = 2;
   }

   printf("Usage:\n> dataset: %s\n> file_out: %s\n> k: %d\n> chunkSize: %ld\n> nt: %d\n", dataset, file_out, k, chunkSize, nt);

   /*
      (!) LEITURA DE SEQUÊNCIAS (FASTA)
   */
   struct read *rd;
   lint st = time(NULL);
   puts("\t... Reading sequences ...");
   rd = (struct read*)calloc(1, sizeof(struct read));
   struct seq *seq = ReadFASTASequences(argv[1], &gnN, &gnS, rd, 1);
   // ReadFASTASequences(argv[1], &gnN, &gnS, rd, 1);

   lint et = time(NULL);
   printf("> Reading time: %ld\n", (et - st));
   printf("> nS: %ld, nN: %ld\n", gnS, gnN);

   /*
      (!) DIVIDINDO OS DADOS EM CHUNKS
   */
   int nChunk = floor(gnS / chunkSize);
   printf("> nChunk: %d\n", nChunk);
   int i = 0;
   
   chunk = (struct read *)calloc(nChunk, sizeof(struct read));
   nS = (lint*)calloc(nChunk, sizeof(lint));
   nN = (lint*)calloc(nChunk, sizeof(lint));

   SelectChunk(chunk, nChunk, rd, chunkSize, chunkSize, gnS, nS, gnN, nN, nt);

   offset = floor(nChunk/devCount);
   printf("> offset: %d\n", offset);

   pthread_t threads[devCount];

   for (i = 0; i < devCount; i++)
   {
      // problema de cast (void*)i
      pthread_create(&threads[i], NULL, LaunchKmer, (void*)&i);
   }

   for (i = 0; i < devCount; i++)
   {
      pthread_join(threads[i], NULL);
   }

   int threadRemain = nChunk - (offset * devCount);
   printf("> threadRemain: %d\n", threadRemain);

   if (threadRemain > 0)
   {
      kmer_main(&chunk[nChunk-1], nN[nChunk-1], nS[nChunk-1], k, device);
   }

   int chunkRemain = abs(gnS - (nChunk*chunkSize));
   lint rnS, rnN;
   printf("> chunkRemain: %d\n", chunkRemain);
   struct read *chunk_remain;

   if(chunkRemain) {
      chunk_remain = SelectChunkRemain(rd, chunkSize, nChunk, chunkRemain, gnS, &rnS, gnN, &rnN, nt);
      kmer_main(chunk_remain, rnN, rnS, k, device);
   }

   st = time(NULL);
   PrintFreq(seq, chunk, nChunk, chunkSize);
   et = time(NULL);

   if(chunkRemain) {
      PrintFreq(seq, chunk_remain, 1, rnS);
   }
   
   printf("> Writing time: %ld\n", (et-st));
   return 0;
}
