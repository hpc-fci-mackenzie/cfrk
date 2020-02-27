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
__global__ void ComputeIndex(char *Seq, struct counter *counter, const int k, lint nN, ushort offset, int *n_combination)
{
   lint idx = threadIdx.x + (blockDim.x * blockIdx.x);

   lint start = idx * offset;
   lint end   = start + offset;
   int n_combination_counter = 0;

   for(lint id = start; id < end; id++)
   {
      int index = 0;
      if (id < nN)
      {
         for( lint i = 0; i < k; i++ )
         {
            char nuc = Seq[i + id];
            if (nuc != -1) //Verifica se ha alguem que nao deve ser processado
            {
               index += nuc * powf(4, ((k-1)-i));
            }
            else
            {
//               index = -1;
               break;
            }
         }//End for i
         __threadfence();
         for (int t = 0; t < *n_combination; t++){
             if (*counter[t].index == -1){
                 atomicAdd(counter[t].index, index);// Value of the combination
                 atomicAdd(counter[t].Freq, 1);// Value of the combination
             } else {
                 atomicAdd(counter[t].index, 1);// Value of the combination
             }
         }
         __syncthreads();
         if (n_combination_counter >= *n_combination) printf("ERROR: Counter greater them combination");

      }
   }//End for id
}

//Compute k-mer frequency
// Change Index and Freq vectors by Struct
/*
 *  struct X {
 *      int index;
 *      int count = 0;
 *  }
 * */
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
//__global__ void ComputeFreqNew(struct counter *counter, lint *start, int *length, ushort offset, int k, lint nS)
//{
//
//   int blx = blockIdx.x;
//
//   int st = blx * offset;
//   int nd = st + offset;
//
//   for (int i = st; i < nd; i++)
//   {
//      int idx = start[i] + threadIdx.x;
//      int id_freq = (fourk * i) + Index[idx];
//      if (threadIdx.x < length[i]-1)
//      {
//         atomicAdd(&counter->Freq, 1);
//      }
//   }
//}
