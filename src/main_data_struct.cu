/*
CFRK-MT - Contabilizador da Frequencia de Repetica de kmer (Multi GPU version)
Developer: Fabricio Gomes Vilasboas
Istitution: National Laboratory for Scientific Computing
*/

#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>
#include <string.h>
#include "kmer_data_struct.cuh"
#include "tipos_data_struct.h"
#include "fastaIO_data_struct.h"

//*** Global Variables ***
struct chunk *chunk;
lint *n_sequence, *n_concat_sequence_length;
int offset;
int device;
int k;
char file_out[512];
//************************
/*
void PrintFreq(struct seq *seq, struct chunk *pchunk, int nChunk, int chunkSize)
{
   FILE *out;
   int cont = 0;
   int cont_seq = 0;
   char str[32];
   lint fourk = pow(4,k);

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
            // printf("%s", str);
            fwrite(str, sizeof(char), sizeof(str), out);
         }
         cont++;
      }
   }
   fclose(out);
}
 */

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

struct chunk* SelectChunkRemain(struct chunk *rd, ushort chunkSize, ushort id_chunk, lint max, lint gnS, lint *n_sequence, lint gnN, lint *n_concat_sequence_length, int nt, int k)
{
   struct chunk *chunk;
   lint i;
   lint j;
   lint length = 0;

   // Size to be allocated
   for (i = 0; i < max; i++)
   {
      lint id = chunkSize*id_chunk + i;
      if (id > gnS-1)
      {
         break;
      }
      length += rd->length[id]+1;
   }
   int size_of_index_vector = length - k + 1;


   cudaMallocHost((void**)&chunk, sizeof(struct chunk));
   cudaMallocHost((void**)&chunk->data, sizeof(char)*length);
   cudaMallocHost((void**)&chunk->length, sizeof(int)*chunkSize);
   cudaMallocHost((void**)&chunk->start, sizeof(lint)*chunkSize);
   cudaMallocHost((void**)&chunk[id_chunk].counter, sizeof(struct counter)*size_of_index_vector);


    // Copy rd->data to chunk->data
   lint start = rd->start[chunkSize*id_chunk];
   lint end = start + (lint)length;

   #pragma omp parallel for num_threads(nt)
   for (j = start; j < end; j++)
   {
      chunk->data[j-start] = rd->data[j];
   }

   chunk->length[0] = rd->length[chunkSize*id_chunk];
   chunk->start[0] = 0;

   // Copy start and length
   for (i = 1; i < max; i++)
   {
      lint id = chunkSize*id_chunk + i;
      chunk->length[i] = rd->length[id];
      chunk->start[i] = chunk->start[i-1]+(chunk->length[i-1]+1);
   }

   *n_concat_sequence_length = length;
   *n_sequence = max;
return chunk;
}


void SelectChunk(struct chunk *chunk, const int nChunk, struct chunk *rd, ushort chunkSize, lint max, lint gnS, lint *n_sequence, lint gnN, lint *n_concat_sequence_length, int nt, int k)
{
   lint i, j, id_chunk;

   for (id_chunk = 0; id_chunk < nChunk; id_chunk++)
   {
      lint length = 0;

      // Size to be allocated
      for (i = 0; i < max; i++)
      {
         lint id = chunkSize*id_chunk + i;
         if (id > gnS-1)
         {
            break;
         }
         length += rd->length[id]+1;
      }
//       int size_of_index_vector = length - k + 1; // c = L - k + 1 296/1024

//      cudaMallocHost((void**)&chunk[id_chunk].data, sizeof(char)*length);
//      cudaMallocHost((void**)&chunk[id_chunk].length, sizeof(int)*chunkSize);
//      cudaMallocHost((void**)&chunk[id_chunk].start, sizeof(lint)*chunkSize);
//      cudaMallocHost((void**)&chunk[id_chunk].n_combination, sizeof(int));
//      cudaMallocHost((void**)&chunk[id_chunk].counter, sizeof(struct counter)*chunkSize);

       chunk[id_chunk].data = (char *) malloc(sizeof(char)*length);
       chunk[id_chunk].length = (int *) malloc(sizeof(int)*chunkSize);
       chunk[id_chunk].start = (lint *) malloc(sizeof(lint)*chunkSize);
       chunk[id_chunk].n_combination = (int *) malloc(sizeof(int));
       chunk[id_chunk].counter = (struct counter *) malloc(sizeof(struct counter)*chunkSize);

      // Copy rd->data to chunk->data
      lint start = rd->start[chunkSize*id_chunk];
      lint end = start + (lint)length;
      #pragma omp parallel for num_threads(nt)
      for (j = start; j < end; j++)
      {
         chunk[id_chunk].data[j-start] = rd->data[j];
      }

      chunk[id_chunk].length[0] = rd->length[chunkSize*id_chunk];
      chunk[id_chunk].start[0] = 0;
      *(chunk[id_chunk].n_combination) = *(rd->n_combination);
      for (lint k = 0; k < max; k++)
      {
//         cudaMallocHost((void**)&chunk[id_chunk].counter[k].index, sizeof(int)**chunk[id_chunk].n_combination);
//         cudaMallocHost((void**)&chunk[id_chunk].counter[k].frequency, sizeof(int)**chunk[id_chunk].n_combination);

          chunk[id_chunk].counter[k].index = (int *) malloc(sizeof(int)**chunk[id_chunk].n_combination);
          chunk[id_chunk].counter[k].frequency = (int *) malloc(sizeof(int)**chunk[id_chunk].n_combination);

//         for(int d = 0; d < *chunk[id_chunk].n_combination; d++)
//         {
//            chunk[id_chunk].counter[k].index[d] = -1;
//            chunk[id_chunk].counter[k].frequency[d] = 0;
//         }
      }
      // Copy start and length
      for (i = 1; i < max; i++)
      {
         lint id = chunkSize*id_chunk + i;
         chunk[id_chunk].length[i] = rd->length[id];
         chunk[id_chunk].start[i] = chunk[id_chunk].start[i-1]+(chunk[id_chunk].length[i-1]+1);
      }

      n_concat_sequence_length[id_chunk] = length;
      n_sequence[id_chunk] = max;
   }
}

void *LaunchKmer(void* threadId)
{

   lint tid = (lint)threadId;
   int start = tid * offset;
   int end = start+offset;

   printf("\t\ttid: %d\n", (int) tid);

   DeviceInfo(tid);

   int i = 0;
   for (i = start; i < end; i++)
   {
      kmer_main(&chunk[i], n_concat_sequence_length[i], n_sequence[i], k, device);
      cudaStreamSynchronize(0);
//      cudaFreeHost(chunk[i].data);
//      cudaFreeHost(chunk[i].length);
//      cudaFreeHost(chunk[i].start);
//      cudaFreeHost(chunk[i].n_combination);
//      cudaFreeHost(chunk[i].counter);
   }

return NULL;
}

int main(int argc, char* argv[])
{

   lint gnN, gnS, chunkSize = 8192;
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
   DeviceInfo(device);

   strcpy(file_out, argv[2]);

   printf("\ndataset: %s, out: %s, k: %d, chunkSize: %d\n", argv[1], file_out, k, (int) chunkSize);

   lint st = time(NULL);
   puts("\n\n\t\tReading seqs!!!");
   struct chunk *rd;
//   cudaMallocHost((void**)&rd, sizeof(struct chunk));
    rd = (struct chunk*) malloc(sizeof(struct chunk));
   struct seq *seq = ReadFASTASequences(argv[1], &gnN, &gnS, rd, 1, k);
   printf("\nnS: %ld, nN: %ld\n", gnS, gnN);
   lint et = time(NULL);

   printf("\n\t\tReading time: %ld\n", (et-st));

   int nChunk = floor(gnS/chunkSize);
   lint i = 0;
   cudaMallocHost((void**)&chunk, sizeof(struct chunk)*nChunk);
   cudaMallocHost((void**)&n_sequence, sizeof(lint)*nChunk);
   cudaMallocHost((void**)&n_concat_sequence_length, sizeof(lint)*nChunk);
   SelectChunk(chunk, nChunk, rd, chunkSize, chunkSize, gnS, n_sequence, gnN, n_concat_sequence_length, nt, k);

   device = SelectDevice(devCount);
   offset = floor(nChunk/devCount);
   pthread_t threads[devCount];

   for (i = 0; i < devCount; i++)
   {
      pthread_create(&threads[i], NULL, LaunchKmer, (void*)i);
   }

   for (i = 0; i < devCount; i++)
   {
      pthread_join(threads[i], NULL);
   }
   for (int k = 0; k < nChunk; k++)
   {
       for (lint t = 0; t < chunkSize; t ++)
       {
           for (int q = 0; q < *(chunk[k].n_combination); q++)
           {
//               if ( chunk[k].counter[t].index[q] != -1)
                   printf("Index: %d Frequency: %d\n", chunk[k].counter[t].index[q], chunk[k].counter[t].frequency[q]);
           }
       }
   }

//   int threadRemain = nChunk - (offset*devCount);
//   if (threadRemain > 0)
//   {
//      kmer_main(&chunk[nChunk-1], n_concat_sequence_length[nChunk-1], n_sequence[nChunk-1], k, device);
//   }
   /*
   int chunkRemain = abs(gnS - (nChunk*chunkSize));
   lint rnS, rnN;
   struct chunk *chunk_remain = SelectChunkRemain(rd, chunkSize, nChunk, chunkRemain, gnS, &rnS, gnN, &rnN, nt, k);
   kmer_main(chunk_remain, rnN, rnS, k, device);
   */
   // st = time(NULL);
//    PrintFreq(seq, chunk, nChunk, chunkSize);
   // et = time(NULL);
   // PrintFreq(seq, chunk_remain, 1, rnS);
   // printf("\n\t\tWriting time: %ld\n", (et-st));

return 0;
}
