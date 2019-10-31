#include <stdio.h>
// #include <cuda.h>
#include "tipos.h"

/*
   The execution configuration is specified by inserting an expression of the form

      <<< Dg, Db, Ns, S >>> between the function name and the parenthesized argument list, where:

   - Dg is of type dim3 (see Section B.3.2) and specifies the dimension and size of the grid, such that Dg.x * Dg.y * Dg.z equals the number of blocks being launched; Dg.z must be equal to 1 for devices of compute capability 1.x;
   - Db is of type dim3 (see Section B.3.2) and specifies the dimension and size of each block, such that Db.x * Db.y * Db.z equals the number of threads per block;
   - Ns is of type size_t and specifies the number of bytes in shared memory that is dynamically allocated per block for this call in addition to the statically allocated memory; this dynamically allocated memory is used by any of the variables declared as an external array as mentioned in Section B.2.3; Ns is an optional argument which defaults to 0;
   - S is of type cudaStream_t and specifies the associated stream; S is an optional argument which defaults to 0.
*/

//Set Matrix values
void SetMatrix(int *Mat, ushort offset, int val, int nF)
{
   /*
   lint idx = threadIdx.x + (blockDim.x * blockIdx.x);

   lint start = idx * offset;
   lint end   = start + offset;

   for(lint id = start; id < end; id++)
   {
      if (id < nF)
         Mat[idx] = val;
   }
   */
}

//Compute k-mer index
void ComputeIndex(char *Seq, int *Index, const int k, lint nN, ushort offset)
{
   /*
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
         Index[id] = index;// Value of the combination
      }
   }//End for id
   */
}

//Compute k-mer frequency
void ComputeFreq(int *Index, int *Freq, lint *start, int *length, ushort offset, int fourk, lint nS, lint nN)
{
   /*
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
   */
}

//New way to compute k-mer frequency
void ComputeFreqNew(int *Index, int *Freq, lint *start, int *length, ushort offset, int fourk, lint nS)
{
   /*
   int blx = blockIdx.x;

   int st = blx * offset;
   int nd = st + offset;

   for (int i = st; i < nd; i++)
   {
      int idx = start[i] + threadIdx.x;
      int id_freq = (fourk * i) + Index[idx];
      if (threadIdx.x < length[i]-1)
      {
         atomicAdd(&Freq[id_freq], 1);
      }
   }
   */
}
