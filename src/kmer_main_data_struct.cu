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
   int n_combination = *(rd->n_combination);
   cudaSetDevice(device);
   GetDeviceProp(device, &maxGridSize, &maxThreadDim, &deviceMemory);

//---------------------------------------------------------------------------

   size[0] = n_concat_sequence_length * sizeof(char);// d_Seq and Seq size
   size[1] = n_sequence * sizeof(int);  // d_length
   size[2] = n_sequence * sizeof(lint); // d_start
   size[3] = n_sequence * sizeof(struct counter); // d_reads
   size[4] = n_combination * sizeof(int);
   totalsize = size[0] + (size[1] * 2) + size[2] + size[3] + (size[3] * 2 * size[4]);

   if (totalsize > deviceMemory)
   {
      printf("\n\n\t\t\t[Error] There is no enough space on GPU memory\n");
      printf("\t\t\t[Error] Required memory: %ld; Available memory: %ld\n", totalsize, deviceMemory);
      exit(1);
   }

//---------------------------------------------------------------------------

    fprintf(stderr, "Memory Allocation\n");
   if ( cudaMalloc ((void**)&d_Seq, size[0])    != cudaSuccess )                                                        fprintf(stderr, "\n[Error 1] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMalloc ((void**)&d_length, size[1]) != cudaSuccess )                                                        fprintf(stderr, "\n[Error 2] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMalloc ((void**)&d_start, size[2])  != cudaSuccess )                                                        fprintf(stderr, "\n[Error 3] %s\n", cudaGetErrorString(cudaGetLastError()));

   struct counter *h_counter = (struct counter *) malloc(size[3]);
   memcpy(h_counter, rd->counter, size[3]);

   for(i = 0; i < n_sequence; i++)
   {
       if( cudaMalloc((void**)&(h_counter[i].index), size[5]) != cudaSuccess )                                          fprintf(stderr, "\n[Error 4-%d] %s\t", i, cudaGetErrorString(cudaGetLastError()));
       cudaMemcpyAsync(h_counter[i].index, rd->counter[i].index, size[4], cudaMemcpyHostToDevice);

       if( cudaMalloc((void**)&(h_counter[i].frequency), size[5]) != cudaSuccess )                                      fprintf(stderr, "[Error 5-%d] %s", i, cudaGetErrorString(cudaGetLastError()));
       cudaMemcpyAsync(h_counter[i].frequency, rd->counter[i].frequency, size[4], cudaMemcpyHostToDevice);

   }
    if ( cudaMalloc ((void**)&d_counter, size[3]) != cudaSuccess )                                                      fprintf(stderr, "\n[Error 6] %s\n", cudaGetErrorString(cudaGetLastError()));
    if ( cudaMemcpyAsync(d_counter, h_counter, size[3], cudaMemcpyHostToDevice) != cudaSuccess)                         fprintf(stderr, "\n[Error 7] %s\n", cudaGetErrorString(cudaGetLastError()));

//************************************************
    fprintf(stderr, "Memory Mapping\n");
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
    fprintf(stderr, "Memory Copy\n");
    if ( cudaMemcpyAsync(d_Seq, rd->data, size[0], cudaMemcpyHostToDevice) != cudaSuccess)                              fprintf(stderr, "[Error 8] %s\n", cudaGetErrorString(cudaGetLastError()));
    if ( cudaMemcpyAsync(d_length, rd->length, size[1], cudaMemcpyHostToDevice) != cudaSuccess)                         fprintf(stderr, "[Error 9] %s\n", cudaGetErrorString(cudaGetLastError()));
    if ( cudaMemcpyAsync(d_start, rd->start, size[2], cudaMemcpyHostToDevice) != cudaSuccess)                           fprintf(stderr, "[Error 10] %s\n", cudaGetErrorString(cudaGetLastError()));

//************************************************
    fprintf(stderr, "Kernel Execution\n");
   // SetMatrix<<<grid[0], block[0]>>>(d_counter, offset[0], n_concat_sequence_length);
   ComputeFrequency<<<grid[0], block[0]>>>(d_Seq, d_counter, d_start, d_length, k, n_concat_sequence_length, offset[0], n_sequence, n_combination);

//************************************************
    fprintf(stderr, "Memory Copy back\n");
    if ( cudaMemcpy(h_counter, d_counter, size[3], cudaMemcpyDeviceToHost) != cudaSuccess)                              fprintf(stderr, "[Error 11] %s\n", cudaGetErrorString(cudaGetLastError()));
    for(i = 0; i < n_sequence; i++)
    {
        if(cudaMemcpy(rd->counter[i].index, h_counter[i].index, size[4], cudaMemcpyDeviceToHost) != cudaSuccess)        fprintf(stderr, "\n[Error 12-%d] %s\t", i, cudaGetErrorString(cudaGetLastError()));
        if(cudaMemcpy(rd->counter[i].frequency, h_counter[i].frequency, size[4],cudaMemcpyDeviceToHost) != cudaSuccess) fprintf(stderr, "[Error 13-%d] %s", i, cudaGetErrorString(cudaGetLastError()));
    }
    memcpy(rd->counter, h_counter, size[3]);
//************************************************
    fprintf(stderr, "\nFree Memory\n");
    printf("\nTotalSize: %d, DeviceMemory: %d\n", totalsize, deviceMemory);
   cudaFree(d_Seq);
   cudaFree(d_counter);
   cudaFree(d_start);
   cudaFree(d_length);
   free(h_counter);

//---------------------------------------------------------------------------

   //printf("\nFim kmer_main\n");
}
