/*
CFRK-MT - Contabilizador da Frequencia de Repetica de kmer (Multi GPU version)
Developer: Fabricio Gomes Vilasboas
Istitution: National Laboratory for Scientific Computing
*/
#define BILLION  1000000000.0
#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>
#include <string.h>
#include "kmer.cuh"
#include "tipos.h"
#include "fastaIO.h"

//*** Global Variables ***
struct read *chunk;
lint *nS, *nN;
int offset;
int device;
int k;
char file_out[512];
//************************

void PrintFreq(struct seq *seq, struct read *pchunk, int nChunk, int chunkSize, char *mode)
{
   FILE *out;
   int cont = 0;
   int cont_seq = 0;
   char str[256];
   lint fourk = pow(4,k);

   out = fopen(file_out, mode);

   for (int j = 0; j < nChunk; j++)
   {
      cont = 0;
      cont_seq = 0;
      for (int i = 0; i < (chunkSize*fourk); i++)
      {
         if (i % fourk == 0 && i != 0)
         {
            cont = 0;
            fprintf(out, "\t\t\t\t\n");
            cont_seq++;
         }
         if (pchunk[j].Freq[i] > 0)
         {
             fprintf(out, "%d:%d ", cont, pchunk[j].Freq[i]);
         }
         cont++;
      }
//      fprintf(out, "\t\t\t\t\n");
   }
   fclose(out);
}

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

struct read* SelectChunkRemain(struct read *rd, ushort chunkSize, ushort it, lint max, lint gnS, lint *nS, lint gnN, lint *nN, int nt)
{
   struct read *chunk;
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

   cudaMallocHost((void**)&chunk, sizeof(struct read));
   cudaMallocHost((void**)&chunk->data, sizeof(char)*length);
   cudaMallocHost((void**)&chunk->length, sizeof(int)*chunkSize);
   cudaMallocHost((void**)&chunk->start, sizeof(lint)*chunkSize);

   // Copy rd->data to chunk->data
   lint start = rd->start[chunkSize*it];
   lint end = start + (lint)length;

   #pragma omp parallel for num_threads(nt)
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

      cudaMallocHost((void**)&chunk[it].data, sizeof(char)*length);
      cudaMallocHost((void**)&chunk[it].length, sizeof(int)*chunkSize);
      cudaMallocHost((void**)&chunk[it].start, sizeof(lint)*chunkSize);

      // Copy rd->data to chunk->data
      lint start = rd->start[chunkSize*it];
      lint end = start + (lint)length;
      #pragma omp parallel for num_threads(nt)
      for (j = start; j < end; j++)
      {
         chunk[it].data[j-start] = rd->data[j];
      }

      chunk[it].length[0] = rd->length[chunkSize*it];
      chunk[it].start[0] = 0;

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

   lint tid = (lint)threadId;
   int start = tid * offset;
   int end = start+offset;

   //printf("\t\ttid: %d\n", tid);

   //DeviceInfo(tid);

   int i = 0;
   for (i = start; i < end; i++)
   {
      kmer_main(&chunk[i], nN[i], nS[i], k, device);
      cudaStreamSynchronize(0);
      cudaFreeHost(chunk[i].data);
      cudaFreeHost(chunk[i].length);
      cudaFreeHost(chunk[i].start);
   }

return NULL;
}

int main(int argc, char* argv[])
{

   lint gnN, gnS, chunkSize = 8192;
    struct timespec start, end;
    int devCount;
   int nt = 12;

   if ( argc < 4)
   {
      printf("Usage: ./cfrk [dataset.fasta] [file_out.cfrk] [k] <number of threads: Default 12> <chunkSize: Default 8192>");
      return 1;
   } 
   cudaDeviceReset();
   
   k = atoi(argv[3]);
   if (argc == 5)
      nt = atoi(argv[4]);
   if (argc == 6)
      chunkSize = atoi(argv[5]);

   cudaGetDeviceCount(&devCount);
   //DeviceInfo(device);

   strcpy(file_out, argv[2]);

   //printf("\ndataset: %s, out: %s, k: %d, chunkSize: %d\n", argv[1], file_out, k, chunkSize);

   lint st = time(NULL);
   struct read *rd;
    clock_gettime(CLOCK_REALTIME, &start);

    cudaMallocHost((void**)&rd, sizeof(struct read));
   struct seq *seq = ReadFASTASequences(argv[1], &gnN, &gnS, rd, 1);

    clock_gettime(CLOCK_REALTIME, &end);
    printf("%d,%d,%.10f,", k, nt, (end.tv_sec - start.tv_sec) +
                                  (end.tv_nsec - start.tv_nsec) / BILLION);

   int nChunk = floor(gnS/chunkSize);
   int i = 0;
   cudaMallocHost((void**)&chunk, sizeof(struct read)*nChunk);
   cudaMallocHost((void**)&nS, sizeof(lint)*nChunk);
   cudaMallocHost((void**)&nN, sizeof(lint)*nChunk);
    clock_gettime(CLOCK_REALTIME, &start);
   SelectChunk(chunk, nChunk, rd, chunkSize, chunkSize, gnS, nS, gnN, nN, nt);
    clock_gettime(CLOCK_REALTIME, &end);
    printf("%.10f,", (end.tv_sec - start.tv_sec) +
                     (end.tv_nsec - start.tv_nsec) / BILLION);
   device = SelectDevice(devCount);
   offset = floor(nChunk/devCount);
   pthread_t threads[devCount];

    clock_gettime(CLOCK_REALTIME, &start);
   for (i = 0; i < devCount; i++)
   {
      pthread_create(&threads[i], NULL, LaunchKmer, (void*)i);
   }

   for (i = 0; i < devCount; i++)
   {
      pthread_join(threads[i], NULL);
   }

   int threadRemain = nChunk - (offset*devCount);
   if (threadRemain > 0)
   {
      kmer_main(&chunk[nChunk-1], nN[nChunk-1], nS[nChunk-1], k, device);
   }
    clock_gettime(CLOCK_REALTIME, &end);
    printf("%.10f,", (end.tv_sec - start.tv_sec) +
                     (end.tv_nsec - start.tv_nsec) / BILLION);

   int chunkRemain = abs(gnS - (nChunk*chunkSize));
   lint rnS, rnN;
    clock_gettime(CLOCK_REALTIME, &start);
   struct read *chunk_remain = SelectChunkRemain(rd, chunkSize, nChunk, chunkRemain, gnS, &rnS, gnN, &rnN, nt);
   kmer_main(chunk_remain, rnN, rnS, k, device);
    clock_gettime(CLOCK_REALTIME, &end);
    printf("%.10f,", (end.tv_sec - start.tv_sec) +
                     (end.tv_nsec - start.tv_nsec) / BILLION);

    clock_gettime(CLOCK_REALTIME, &start);
   PrintFreq(seq, chunk, nChunk, chunkSize, "w");
   // et = time(NULL);
   PrintFreq(seq, chunk_remain, 1, rnS, "a");
    clock_gettime(CLOCK_REALTIME, &end);
    printf("%.10f", (end.tv_sec - start.tv_sec) +
                     (end.tv_nsec - start.tv_nsec) / BILLION);

return 0;
}
