#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include "kmer.cuh"
#include "tipos.h"
#include "fastaIO.h"

int SelectDevice(int devCount)
{

   int i, device = 0;
   cudaDeviceProp prop[devCount];

   if (devCount > 0)
   {
      for (i = 0; i < devCount; i++)
      {
         cudaGetDeviceProperties(&prop[i], i);
      }

      for (i = 0; i < devCount; i++)
      {
         if (prop[i].totalGlobalMem > prop[device].totalGlobalMem)
         {
            device = i;
         }
      }
   }
   else
      return 0;

return device;
}

struct read* SelectChunk(struct read *rd, ushort chunkSize, ushort it, lint gnS, lint *nS, lint gnN, lint *nN)
{
   struct read *chunk;
   int i;
   lint length = 0;

   // Size to be allocated
   for (i = 0; i < chunkSize; i++)
   {
      int id = chunkSize*it + i;
      if (id > gnS-1)
      {
         break;
      }
      length += rd->length[id]+1;
   }

   cudaMallocHost((void**)&chunk, sizeof(struct read));
   cudaMallocHost((void**)&chunk->data, sizeof(short)*length);
   cudaMallocHost((void**)&chunk->length, sizeof(int)*chunkSize);
   cudaMallocHost((void**)&chunk->start, sizeof(int)*chunkSize);

   // Copy rd->data to chunk->data
   lint start = rd->start[chunkSize*it];
   lint end = start + length;
   for (i = start; i < end; i++)
   {
      chunk->data[i-start] = rd->data[i];
   }

   // Copy start and length
   for (i = 0; i < chunkSize; i++)
   {
      int id = chunkSize*it + i;
      chunk->length[i] = rd->length[id];
      chunk->start[i] = rd->start[id];
   }

   *nN = length;
   *nS = chunkSize;
return chunk;
}

int main(int argc, char* argv[])
{

   ushort k;
   int device;
   lint gnN, gnS, nN, nS, chunkSize = 10000;
   int devCount;

   if ( argc < 3)
   {
      printf("Usage: ./kmer <dataset.fasta> <k> ");
      return 1;
   }

   cudaDeviceReset();
   
   k = atoi(argv[2]);

   cudaGetDeviceCount(&devCount);
   device = SelectDevice(devCount);

   printf("dataset: %s, k: %d", argv[1], k);

   lint st = time(NULL);
   puts("\n\n\t\tReading seqs!!!");
   struct read *rd;
   cudaMallocHost((void**)&rd, sizeof(struct read));
   ReadFASTASequences(argv[1], &gnN, &gnS, rd, 1);
   printf("\nnS: %ld, nN: %ld\n", gnS, gnN);
   lint et = time(NULL);

   printf("\n\t\tReading time: %ld\n", (et-st));

   struct read *chunk[(gnS/chunkSize)+1];
   for (int i = 0; i < (gnS/chunkSize)+1; i++)
   {
      chunk[i] = SelectChunk(rd, chunkSize, i, gnS, &nS, gnN, &nN);
      kmer_main(chunk[i], nN, nS, k, device);
      //cudaDeviceReset();
   }

return 0;
}
