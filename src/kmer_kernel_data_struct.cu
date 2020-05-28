#include <stdio.h>
#include <cuda.h>
//#include <string.h>
#include "tipos_data_struct.h"

//__device__ struct counter *sd_counter;

//Set Matrix values
__global__ void SetMatrix(struct counter *Mat, ushort offset, int nF, int id_sequence) {
    lint idx = threadIdx.x + (blockDim.x * blockIdx.x);

    lint start = idx * offset;
    lint end = start + offset;

    for (lint id = start; id < end; id++) {
        if (id < nF) {
            Mat[id_sequence].index[idx] = -1;
            Mat[id_sequence].frequency[idx] = 0;
        }
    }
}

//Compute k-mer index
__global__ void ComputeFrequency(char *Seq, struct counter *d_counter, lint *d_start, int *d_length, const int k, lint nN, ushort offset, lint n_sequence, int n_combination)
{

    int idx = threadIdx.x + (blockDim.x * blockIdx.x);

    if (idx < n_sequence)
    {
        int end = d_start[idx] + (d_length[idx] + 1);

        for (int i = d_start[idx]; i < end; i++)
        {
            int index = 0;
            for ( lint j = 0; j < k ; j++ )
            {
                if ( Seq[i + j] != -1 )
                {
                    index += Seq[i + j] * exp10f(j);
                    break;
                }
                else
                {
                    index = -1;
                    break;
                }
            }
            __threadfence();
            for (int t = 0; t < n_combination && index > -1; t++)
            {
                if (d_counter[idx].index[t] == -1)
                {
                    atomicAdd(&(d_counter[idx].index[t]), index + 1);// Index
                    atomicAdd(&(d_counter[idx].frequency[t]), 1);// Value of the combination
//                    d_counter[id_sequence].index[t] = index + 1;
//                    d_counter[id_sequence].frequency[t]++;
                    break;
                }
                else if (d_counter[idx].index[t] == index)
                {
                    atomicAdd(&(d_counter[idx].frequency[t]), 1);// Value of the combination
//                    d_counter[id_sequence].frequency[t]++;
                    break;
                }
            }
            __syncthreads();
        }
    }
    /*
   lint idx =  threadIdx.x + (blockDim.x * blockIdx.x);

   lint start = idx * offset;
   lint end   = start + offset;

   for(lint id = start; id < end; id++) // start ao end da Thread
   {
      
      lint index = -1;
      if (id < nN)
      {
         lint id_sequence;
         lint p;
         for (p = 0; p < n_sequence; p++) // Unnecessary
         {
            if(d_start[p] < id && id < (d_start[p] + d_length[p]))
            {
               id_sequence = p;
            }
         }
         // TODO: Change to Key-Value logic
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
//            __threadfence();
            for (int t = 0; t < n_combination; t++){
                if (d_counter[id_sequence].index[t] == -1)
                {
                    atomicAdd(&(d_counter[id_sequence].index[t]), index);// Value of the combination
                    atomicAdd(&(d_counter[id_sequence].frequency[t]), 1);// Value of the combination
//                    d_counter[id_sequence].index[t] = index + 1;
//                    d_counter[id_sequence].frequency[t]++;
                    break;
                } else if (d_counter[id_sequence].index[t] == index)
                {
                    atomicAdd(&(d_counter[id_sequence].frequency[t]), 1);// Value of the combination
//                    d_counter[id_sequence].frequency[t]++;
                    break;
                }
            }
//            __syncthreads();
         }
      }
   }//End for id
     */

}

