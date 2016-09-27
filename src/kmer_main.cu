#include <stdio.h>
#include <cuda.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "tipos.h"
#include "kmer.cuh"

void PrintFreq(int fourk, struct read *rd, int *Freq, lint nS)
{
   int cont = 0;
   int cont_seq = 0;
   for (int i = 0; i < (nS*fourk); i++)
   {
      if (i % fourk == 0)
      {
         cont = 0;
         printf("> %d\n", rd->length[cont_seq]);
         cont_seq++;
      }
      if (Freq[i] != 0)
      {
         printf("%d: %d\n", cont, Freq[i]);
      }
      cont++;
   }
   printf("\n");
}

void GetDeviceProp(uint8_t device, lint *maxGridSize, lint *maxThreadDim, lint *deviceMemory)
{
   cudaDeviceProp prop;

   cudaGetDeviceProperties(&prop, device);

   *maxThreadDim = prop.maxThreadsDim[0];
   *maxGridSize = prop.maxGridSize[0];
   *deviceMemory = prop.totalGlobalMem;
}

void kmer_main(struct read *rd, lint nN, lint nS, int k, ushort device)
{

   int *d_Index;// Index vector
   char *d_Seq;// Seq matrix
   int *Freq, *d_Freq;// Frequence vector
   int fourk;// 4 power k
   lint *d_start;
   int *d_length;// The beggining and the length of each sequence
   lint block[4], grid[4];// Grid config; 0:nN, 1:nS
   lint maxGridSize, maxThreadDim, deviceMemory;// Device config
   ushort offset[4] = {1,1,1,1};
   size_t size[5], totalsize;

   d_Index =  NULL;
   d_Seq = NULL;

   fourk = POW(k);

   cudaSetDevice(device);
   GetDeviceProp(device, &maxGridSize, &maxThreadDim, &deviceMemory);

//---------------------------------------------------------------------------
   size[0] = nN * sizeof(char);// d_Seq and Seq size
   size[1] = nN * sizeof(int); // d_Index and Index size
   size[2] = nS * sizeof(int);  // d_length
   size[3] = nS * fourk * sizeof(int);// Freq and d_Freq
   size[4] = nS * sizeof(lint); // d_start
   totalsize = size[0] + size[1] + (size[2] * 2) + size[3];

   if (totalsize > deviceMemory)
   {
      printf("\n\n\t\t\t[Error] There is no enough space on GPU memory\n");
      printf("\t\t\t[Error] Required memory: %ld; Available memory: %ld\n", totalsize, deviceMemory);
      exit(1);
   }
//---------------------------------------------------------------------------

   if ( cudaMalloc    ((void**)&d_Seq, size[0]) != cudaSuccess) printf("\n[Error 1] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMalloc    ((void**)&d_Index, size[1]) != cudaSuccess) printf("\n[Error 2] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMalloc    ((void**)&d_start, size[4]) != cudaSuccess) printf("\n[Error 3] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMalloc    ((void**)&d_length, size[2]) != cudaSuccess) printf("\n[Error 4] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMalloc    ((void**)&d_Freq, size[3]) != cudaSuccess) printf("\n[Error 5] %s\n", cudaGetErrorString(cudaGetLastError()));

//************************************************
   block[0] = maxThreadDim;
   grid[0] = floor(nN / block[0]) + 1;
   if (grid[0] > maxGridSize)
   {
      grid[0] = maxGridSize;
      offset[0] = (nN / (grid[0] * block[0])) + 1;
   }

   block[1] = maxThreadDim;
   grid[1] = (nS / block[1]) + 1;
   if (grid[1] > maxGridSize)
   {
      grid[1] = maxGridSize;
      offset[1] = (nS / (grid[1] * block[1])) + 1;
   }

   block[2] = maxThreadDim;
   grid[2] = nS;
   if (nS > maxGridSize)
   {
      grid[2] = maxGridSize;
      offset[2] = (nS / grid[2]) + 1;
   }

   int nF = nS*POW(k);
   block[3] = maxThreadDim;
   grid[3] = ((nS*POW(k))/1024)+1;
   if (grid[3] > maxGridSize)
   {
      grid[3] = maxGridSize;
      offset[3] = (nF / (grid[3] * block[3])) + 1;
   }

//************************************************

   if ( cudaMemcpyAsync(d_Seq, rd->data, size[0], cudaMemcpyHostToDevice) != cudaSuccess) printf("[Error 6] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMemcpyAsync(d_start, rd->start, size[4], cudaMemcpyHostToDevice) != cudaSuccess) printf("[Error 7] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMemcpyAsync(d_length, rd->length, size[2], cudaMemcpyHostToDevice) != cudaSuccess) printf("[Error 8] %s\n", cudaGetErrorString(cudaGetLastError()));

//************************************************

   SetMatrix<<<grid[0], block[0]>>>(d_Index, offset[0], -1, nN);
   SetMatrix<<<grid[3], block[3]>>>(d_Freq, offset[3], 0, nF);
   ComputeIndex<<<grid[0], block[0]>>>(d_Seq, d_Index, k, nN, offset[0]);
   //ComputeFreq<<<grid[1], block[1]>>>(d_Index, d_Freq, d_start, d_length, offset[1], fourk, nS, nN);
   ComputeFreqNew<<<grid[2],block[2]>>>(d_Index, d_Freq, d_start, d_length, offset[2], fourk, nS);

   cudaFree(rd);

   if ( cudaMallocHost((void**)&Freq, size[3]) != cudaSuccess) printf("\n[Error 9] %s\n", cudaGetErrorString(cudaGetLastError()));
   if ( cudaMemcpy(Freq, d_Freq, size[3], cudaMemcpyDeviceToHost) != cudaSuccess) printf("\n[Error 10] %s\n", cudaGetErrorString(cudaGetLastError()));

   //PrintFreq(fourk, rd, Freq, nS);

//************************************************
   cudaFree(d_Seq);
   cudaFree(d_Freq);
   cudaFree(d_Index);
   cudaFree(d_start);
   cudaFree(d_length);
   cudaFreeHost(Freq);
//---------------------------------------------------------------------------

   //printf("\nFim kmer_main\n");
}
