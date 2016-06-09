#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <stdint.h>
#include "kmer.cuh"
#include "tipos.h"
#include "fastaIO.h"

void DeviceInfo(uint8_t device)
{
   cudaDeviceProp prop;

   cudaGetDeviceProperties(&prop, device);

   printf("\n\n***** Device information *****\n\n");

   printf("\tId: %d\n", device);
   printf("\tName: %s\n", prop.name);
   printf("\tTotal global memory: %ld\n", prop.totalGlobalMem);
   printf("\tMax grid size: %d, %d, %d\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
   printf("\tMax thread dim: %d, %d, %d\n", prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
   printf("\tWarp size: %d\n", prop.warpSize);
   printf("\tMax threads per multiprocessor: %d\n", prop.maxThreadsPerMultiProcessor);

   printf("\n************************************\n\n");
}

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

struct read* SelectChunk(struct read *rd, ushort chunkSize, ushort it, lint max, lint gnS, lint *nS, lint gnN, lint *nN)
{
   struct read *chunk;
   int i;
   lint length = 0;

   // Size to be allocated
   for (i = 0; i < max; i++)
   {
      int id = chunkSize*it + i;
      if (id > gnS-1)
      {
         break;
      }
      //printf("id: %d\n", id);
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

   chunk->length[0] = rd->length[chunkSize*it];
   chunk->start[0] = 0;
   // Copy start and length
   for (i = 1; i < max; i++)
   {
      int id = chunkSize*it + i;
      chunk->length[i] = rd->length[id];
      chunk->start[i] = chunk->start[i-1]+(chunk->length[i-1]+1);
   }

   *nN = length;
   *nS = max;
return chunk;
}

int main(int argc, char* argv[])
{

   int k;
   int device;
   lint gnN, gnS, nN, nS, chunkSize = 4096;
   int devCount;

   if ( argc < 4)
   {
      printf("Usage: ./kmer [dataset.fasta] [k] <chunkSize: Default 4096>");
      return 1;
   }

   cudaDeviceReset();
   
   k = atoi(argv[2]);
   if (argc == 4)
      chunkSize = atoi(argv[3]);

   cudaGetDeviceCount(&devCount);
   device = SelectDevice(devCount);
   DeviceInfo(device);

   printf("\ndataset: %s, k: %d, chunkSize: %d\n", argv[1], k, chunkSize);

   lint st = time(NULL);
   puts("\n\n\t\tReading seqs!!!");
   struct read *rd;
   cudaMallocHost((void**)&rd, sizeof(struct read));
   ReadFASTASequences(argv[1], &gnN, &gnS, rd, 1);
   printf("\nnS: %ld, nN: %ld\n", gnS, gnN);
   lint et = time(NULL);

   printf("\n\t\tReading time: %ld\n", (et-st));

   int nChunk = floor(gnS/chunkSize);
   struct read *chunk;
   int i = 0;
   for (i = 0; i < nChunk; i++)
   {
      chunk = SelectChunk(rd, chunkSize, i, chunkSize, gnS, &nS, gnN, &nN);
      kmer_main(chunk, nN, nS, k, device);
      cudaFree(chunk->data);
      cudaFree(chunk->length);
      cudaFree(chunk->start);
      cudaFree(chunk);
   }
   int chunkRemain = abs(gnS - (nChunk*chunkSize));
   chunk = SelectChunk(rd, chunkSize, nChunk, chunkRemain, gnS, &nS, gnN, &nN);
   printf("\nnS: %ld, nN: %ld, chunkRemain: %d, it: %d\n", nS, nN, chunkRemain, nChunk);
   kmer_main(chunk, nN, nS, k, device);

return 0;
}
