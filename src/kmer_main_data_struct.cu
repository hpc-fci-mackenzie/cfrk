#include <stdio.h>
#include <cuda.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "tipos_data_struct.h"
#include "kmer_data_struct.cuh"

void GetDeviceProp(uint8_t device, lint *maxGridSize, lint *maxThreadDim, lint *deviceMemory)
{
   cudaDeviceProp prop;

   cudaGetDeviceProperties(&prop, device);

   *maxThreadDim = prop.maxThreadsDim[0];
   *maxGridSize = prop.maxGridSize[0];
   *deviceMemory = prop.totalGlobalMem;
}

void kmer_main(struct chunk *rd, lint n_concat_sequence_length, lint n_sequence, int k, ushort device)
{
   char *d_Seq;// Seq matrix
   lint *d_start;
   int *d_n_combination;
   struct counter *d_counter;
   int *d_length;// The beginning and the length of each sequence
   lint block[4], grid[4];// Grid config; 0:n_concat_sequence_length, 1:n_sequence
   lint maxGridSize, maxThreadDim, deviceMemory;// Device config
   ushort offset[4] = {1,1,1,1};
   size_t size[5], totalsize;
   int sizeOfAllCounters = 0; //
   int i;

   cudaSetDevice(device);
   GetDeviceProp(device, &maxGridSize, &maxThreadDim, &deviceMemory);

//---------------------------------------------------------------------------
   size[0] = n_concat_sequence_length * sizeof(char);// d_Seq and Seq size
   size[1] = n_sequence * sizeof(int);  // d_length
   size[2] = n_sequence * sizeof(lint); // d_start
   size[3] = n_sequence * sizeof(struct counter); // d_reads
   totalsize = size[0] + (size[1] * 2) + size[2] + size[3];

   if (totalsize > deviceMemory)
   {
      printf("\n\n\t\t\t[Error] There is no enough space on GPU memory\n");
      printf("\t\t\t[Error] Required memory: %ld; Available memory: %ld\n", totalsize, deviceMemory);
      exit(1);
   }
//---------------------------------------------------------------------------

   if ( cudaMalloc ((void**)&d_Seq, size[0])    != cudaSuccess ) printf("\n[Error 1] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMalloc ((void**)&d_length, size[1]) != cudaSuccess ) printf("\n[Error 2] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMalloc ((void**)&d_start, size[2])  != cudaSuccess ) printf("\n[Error 3] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMalloc ((void**)&d_n_combination, sizeof(int)) != cudaSuccess ) printf("\n[Error 4] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMalloc ((void**)&d_counter, size[3]) != cudaSuccess )  printf("\n[Error 5] %s\n", cudaGetErrorString(cudaGetLastError()));

   for(i = 0; i < n_sequence; i++)
   {
      if( cudaMalloc((void**)&d_counter[i].index, sizeof(int)**rd->n_combination) != cudaSuccess )      printf("\n[Error 6-%d] %s\n", i, cudaGetErrorString(cudaGetLastError()));
      if( cudaMalloc((void**)&d_counter[i].frequence, sizeof(int)**rd->n_combination) != cudaSuccess )  printf("\n[Error 7-%d] %s\n", i, cudaGetErrorString(cudaGetLastError()));
   }
  
//************************************************
   // Thread mapping for raw data
   block[0] = maxThreadDim;
   grid[0] = floor(n_concat_sequence_length / block[0]) + 1;
   if (grid[0] > maxGridSize)
   {
      grid[0] = maxGridSize;
      offset[0] = (n_concat_sequence_length / (grid[0] * block[0])) + 1;
   }

   // Thread mapping for 
   block[1] = maxThreadDim;
   grid[1] = (n_sequence / block[1]) + 1;
   if (grid[1] > maxGridSize)
   {
      grid[1] = maxGridSize;
      offset[1] = (n_sequence / (grid[1] * block[1])) + 1;
   }

   block[2] = maxThreadDim;
   grid[2] = n_sequence;
   if (n_sequence > maxGridSize)
   {
      grid[2] = maxGridSize;
      offset[2] = (n_sequence / grid[2]) + 1;
   }

   block[3] = maxThreadDim;
   grid[3] = ((n_sequence*POW(k))/1024)+1;
   if (grid[3] > maxGridSize)
   {
      grid[3] = maxGridSize;
      offset[3] = (sizeOfAllCounters / (grid[3] * block[3])) + 1;
   }

//************************************************

   if ( cudaMemcpyAsync(d_Seq, rd->data, size[0], cudaMemcpyHostToDevice) != cudaSuccess)      printf("[Error 8] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMemcpyAsync(d_length, rd->length, size[1], cudaMemcpyHostToDevice) != cudaSuccess) printf("[Error 9] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMemcpyAsync(d_start, rd->start, size[2], cudaMemcpyHostToDevice) != cudaSuccess)   printf("[Error 10] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMemcpyAsync(d_n_combination, rd->n_combination, sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess)   printf("[Error 11] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMemcpyAsync(d_counter, rd->counter, size[3], cudaMemcpyHostToDevice) != cudaSuccess)   printf("[Error 12] %s\n", cudaGetErrorString(cudaGetLastError()));

   for(i = 0; i < n_sequence; i++)
   {
      if ( cudaMemcpyAsync(d_counter[i].index, rd->counter[i].index, sizeof(int)**rd->n_combination, cudaMemcpyHostToDevice) != cudaSuccess)           printf("\n[Error 13-%d] %s\n", i, cudaGetErrorString(cudaGetLastError()));
      if ( cudaMemcpyAsync(d_counter[i].frequence, rd->counter[i].frequence, sizeof(int)**rd->n_combination, cudaMemcpyHostToDevice) != cudaSuccess)   printf("\n[Error 14-%d] %s\n", i, cudaGetErrorString(cudaGetLastError()));
   }

//************************************************

   // SetMatrix<<<grid[0], block[0]>>>(d_counter, offset[0], n_concat_sequence_length);
   ComputeFrequence<<<grid[0], block[0]>>>(d_Seq, d_counter, d_start, d_length, k, n_concat_sequence_length, offset[0], n_sequence, *d_n_combination);


//************************************************

   for (i = 0; i < n_sequence; i++)
   {
      if ( cudaMemcpy(rd->counter[i].index, d_counter[i].index, sizeof(int)**rd->n_combination, cudaMemcpyHostToDevice) != cudaSuccess)           printf("\n[Error 15-%d] %s\n", i, cudaGetErrorString(cudaGetLastError()));
      if ( cudaMemcpy(rd->counter[i].frequence, d_counter[i].frequence, sizeof(int)**rd->n_combination, cudaMemcpyHostToDevice) != cudaSuccess)   printf("\n[Error 16-%d] %s\n", i, cudaGetErrorString(cudaGetLastError()));
   }

//************************************************

   cudaFree(d_Seq);
   cudaFree(d_counter);
   cudaFree(d_start);
   cudaFree(d_length);
   for (i = 0; i < n_sequence; i++)
   {
         cudaFree(d_counter[i].index);
         cudaFree(d_counter[i].frequence);
   }
   cudaFree(d_counter);

//---------------------------------------------------------------------------

   //printf("\nFim kmer_main\n");
}
