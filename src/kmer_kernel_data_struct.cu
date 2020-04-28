#include <stdio.h>
#include <cuda.h>
#include "tipos_data_struct.h"



//Compute k-mer index
__global__ void ComputeFrequency(char *Seq, struct counter *d_counter, lint *d_start, int *d_length, const int k, lint nN, ushort offset, int n_sequence, int n_combination)
{
   int idx =  blockIdx.x;

   int start = idx * offset;
   int end   = start + offset;

   for(lint id = start; id < end; id++)
   {
      
      int index = -1;
      if (id < nN)
      {
         lint id_sequence;
         lint p;
         for (p = 0; p < n_sequence; p++)
         {
            if(d_start[p] < id && id < (d_start[p] + d_length[p]))
            {
               id_sequence = p;
            }
         }
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
            for (int t = 0; t < n_combination; t++){
                if (d_counter[id_sequence].index[t] == -1){
                    atomicAdd(&d_counter[id_sequence].index[t], index);// Value of the combination
                    atomicAdd(&d_counter[id_sequence].frequency[t], 1);// Value of the combination
                    break;
                } else if (d_counter[id_sequence].index[t] == index) {
                    atomicAdd(&d_counter[id_sequence].frequency[t], 1);// Value of the combination
                    break;
                }
            }
            __syncthreads();
         }
      }
   }//End for id

}

