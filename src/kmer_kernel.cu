#include <stdio.h>
#include <cuda.h>
#include "tipos.h"

//Set Matrix values
__global__ void SetMatrix(int *Mat, ushort offset, int val, int nF)
{
   lint idx = threadIdx.x + (blockDim.x * blockIdx.x);

   lint start = idx * offset;
   lint end   = start + offset;

   for(lint id = start; id < end; id++)
   {
      if (id < nF)
         Mat[idx] = val;
   }
}

//Compute k-mer index
__global__ void ComputeIndex(char *Seq, int *Index, const int k, lint nN, ushort offset)
{
   lint idx = threadIdx.x + (blockDim.x * blockIdx.x);

   lint start = idx * offset;
   lint end   = start + offset;

   for(lint id = start; id < end; id++)
   {
      lint index = 0;
      if (id < nN)
      {
         for( lint i = 0; i < k; i++ )
         {
            if (Seq[i + id] != -1) //Verifica se ha alguem que nao deve ser processado
            {
               index += Seq[i + id] * POW( (k - 1) - i );
            }
            else
            {
               index = -1;
               break;
            }
         }//End for i
         Index[id] = index;// Value of the combination
      }
   }//End for id
}

//Compute k-mer frequency
__global__ void ComputeFreq(int *Index, int *Freq, lint *start, int *length, ushort offset, int fourk, lint nS, lint nN)
{

   int idx = threadIdx.x + (blockDim.x * blockIdx.x);

   if (idx < nS)
   {
      int end = start[idx] + (length[idx] + 1);

      for (int i = start[idx]; i < end; i++)
      {
         if (Index[i] != -1 && i < nN)
         {
            int pos = (fourk * idx) + Index[i];
            Freq[pos] += 1;
         }
      }
   }
}

//New way to compute k-mer frequency
__global__ void ComputeFreqNew(int *Index, int *Freq, lint *start, int *length, ushort offset, int fourk, lint nS)
{

   int blx = blockIdx.x;

   int st = blx * offset;
   int nd = st + offset;

   for (int i = st; i < nd; i++)
   {
      int idx = start[i] + threadIdx.x;
      int id_freq = (fourk * i) + Index[idx];
      if (threadIdx.x < length[i])
      {
         //Freq[idx] = 1;
         atomicAdd(&Freq[id_freq], 1);
      }
   }
}
