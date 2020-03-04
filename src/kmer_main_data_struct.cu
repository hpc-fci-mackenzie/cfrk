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

void kmer_main(struct read *rd, lint n_concat_sequence_length, lint n_sequence, int k, ushort device)
{

//   int *d_Index;// Index vector
   char *d_Seq;// Seq matrix
   // int *Freq;
//   int *d_Freq;// Frequency vector
//   int fourk;// 4 power k
int *d_start;
   int *d_length;// The beginning and the length of each sequence
   lint block[4], grid[4];// Grid config; 0:n_concat_sequence_length, 1:n_sequence
   lint maxGridSize, maxThreadDim, deviceMemory;// Device config
   ushort offset[4] = {1,1,1,1};
   size_t size[5], totalsize;

   struct counter *d_counter;
   int n_combination;
   if ( cudaMemcpyAsync(&n_combination, rd->n_combination, sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) printf("[Error 8] %s\n", cudaGetErrorString(cudaGetLastError()));
    if ( cudaMalloc ((void**)&d_counter, sizeof(struct counter)*n_combination) != cudaSuccess ) printf("\n[Error 1] %s\n", cudaGetErrorString(cudaGetLastError()));
    if ( cudaMalloc ((void**)&n_combination, sizeof(int)) != cudaSuccess ) printf("\n[Error 1] %s\n", cudaGetErrorString(cudaGetLastError()));
    if ( cudaMemcpyAsync(d_counter, rd->counter, sizeof(struct counter)*n_combination, cudaMemcpyHostToDevice) != cudaSuccess) printf("[Error 8] %s\n", cudaGetErrorString(cudaGetLastError()));
    for(int d = 0; d < n_combination; d++){
        cudaMalloc((void**)&d_counter[d].index, sizeof(int));
        *d_counter[d].index = -1;
        cudaMalloc((void**)&d_counter[d].Freq, sizeof(int));
        *d_counter[d].Freq = 0;
    }

//    fourk = POW(k);
   // n_possible_combinations = POW(k);

   cudaSetDevice(device);
   GetDeviceProp(device, &maxGridSize, &maxThreadDim, &deviceMemory);

//---------------------------------------------------------------------------
   size[0] = n_concat_sequence_length * sizeof(char);// d_Seq and Seq size
//   size[1] = n_concat_sequence_length * sizeof(int); // d_Index and Index size
   size[2] = n_sequence * sizeof(int);  // d_length
//   size[3] = n_sequence * n_possible_combinations * sizeof(int);// Freq and d_Freq
   size[4] = n_sequence * sizeof(lint); // d_start
   totalsize = size[0] + size[1] + (size[2] * 2) + size[3];

   if (totalsize > deviceMemory)
   {
      printf("\n\n\t\t\t[Error] There is no enough space on GPU memory\n");
      printf("\t\t\t[Error] Required memory: %ld; Available memory: %ld\n", totalsize, deviceMemory);
      exit(1);
   }
//---------------------------------------------------------------------------

   if ( cudaMalloc ((void**)&d_Seq, size[0])    != cudaSuccess ) printf("\n[Error 1] %s\n", cudaGetErrorString(cudaGetLastError()));
//   if ( cudaMalloc ((void**)&d_Index, size[1])  != cudaSuccess ) printf("\n[Error 2] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMalloc ((void**)&d_start, size[4])  != cudaSuccess ) printf("\n[Error 3] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMalloc ((void**)&d_length, size[2]) != cudaSuccess ) printf("\n[Error 4] %s\n", cudaGetErrorString(cudaGetLastError()));
//   if ( cudaMalloc ((void**)&d_Freq, size[3])   != cudaSuccess ) printf("\n[Error 5] %s\n", cudaGetErrorString(cudaGetLastError()));

//************************************************
   block[0] = maxThreadDim;
   grid[0] = floor(n_concat_sequence_length / block[0]) + 1;
   if (grid[0] > maxGridSize)
   {
      grid[0] = maxGridSize;
      offset[0] = (n_concat_sequence_length / (grid[0] * block[0])) + 1;
   }

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

   int nF = n_sequence*POW(k);
   block[3] = maxThreadDim;
   grid[3] = ((n_sequence*POW(k))/1024)+1;
   if (grid[3] > maxGridSize)
   {
      grid[3] = maxGridSize;
      offset[3] = (nF / (grid[3] * block[3])) + 1;
   }

//************************************************

   if ( cudaMemcpyAsync(d_Seq, rd->data, size[0], cudaMemcpyHostToDevice) != cudaSuccess) printf("[Error 6] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMemcpyAsync(d_start, rd->start, size[4], cudaMemcpyHostToDevice) != cudaSuccess) printf("[Error 7] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMemcpyAsync(d_length, rd->length, size[2], cudaMemcpyHostToDevice) != cudaSuccess) printf("[Error 8] %s\n", cudaGetErrorString(cudaGetLastError()));

//   if ( cudaMemcpyAsync(d_Index, rd->index, size[2], cudaMemcpyHostToDevice) != cudaSuccess) printf("[Error 9] %s\n", cudaGetErrorString(cudaGetLastError()));
//   if ( cudaMemcpyAsync(d_Freq, rd->Freq, size[2], cudaMemcpyHostToDevice) != cudaSuccess) printf("[Error 10] %s\n", cudaGetErrorString(cudaGetLastError()));

//************************************************

//   SetMatrix<<<grid[0], block[0]>>>(d_Index, offset[0], -1, n_concat_sequence_length);
//   SetMatrix<<<grid[3], block[3]>>>(d_Freq, offset[3], 0, nF);
   ComputeIndex<<<grid[0], block[0]>>>(d_Seq, d_counter, k, n_concat_sequence_length, offset[0], &n_combination);
   //ComputeFreq<<<grid[1], block[1]>>>(d_Index, d_Freq, d_start, d_length, offset[1], n_possible_combinations, n_sequence, n_concat_sequence_length);
//   ComputeFreqNew<<<grid[2],block[2]>>>(counter d_start, d_length, offset[2], k, n_sequence);


   //cudaFree(rd);

   if ( cudaMallocHost((void**)&rd->counter, sizeof(struct counter)*n_combination) != cudaSuccess) printf("\n[Error 11] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMemcpy(rd->counter, d_counter, sizeof(struct counter)*n_combination, cudaMemcpyDeviceToHost) != cudaSuccess) printf("\n[Error 12] %s\n", cudaGetErrorString(cudaGetLastError()));


//************************************************
   cudaFree(d_Seq);
   // cudaFree(d_Freq);
   // cudaFree(d_Index);
   cudaFree(d_counter);
   cudaFree(d_start);
   cudaFree(d_length);
//---------------------------------------------------------------------------

   //printf("\nFim kmer_main\n");
}
