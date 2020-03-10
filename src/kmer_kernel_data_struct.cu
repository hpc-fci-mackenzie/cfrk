#include <stdio.h>
#include <cuda.h>
#include "tipos_data_struct.h"

//Set Matrix values
__global__ void SetMatrix(struct read *Mat, ushort offset, int val, int nF)
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
__global__ void ComputeFrequence(char *Seq, struct read *d_read, lint *start, int d_length const int k, lint nN, ushort offset)
{
   lint idx =  blockIdx.x;

   lint start = idx * offset;
   lint end   = start + offset;
   int n_combination_counter = 0;

   for(lint id = start; id < end; id++)
   {
      int index = -1;
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
               index = -1;
               break;
            }
         }//End for i
         if(index != -1)
         {
            __threadfence();
            for (int t = 0; t < *n_combination; t++){
                if (*counter[t].index == -1){
                    atomicAdd(counter[t].index, index);// Value of the combination
                    atomicAdd(counter[t].frequence, 1);// Value of the combination
                } else {
                    atomicAdd(counter[t].frequence, 1);// Value of the combination
                }
            }
            __syncthreads();
         }
         
         if (n_combination_counter >= *n_combination) printf("ERROR: Counter greater them combination");

      }
   }//End for id
}

