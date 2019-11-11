#include <stdio.h>
// #include <cuda.h>
#include <math.h>
#include "tipos.h"

//Set Matrix values
void SetMatrix(int *Mat, ushort offset, int val, int nF)
{
}

/*
   Compute k-mer index
   Na primeira fase é realizada a conversão dos valores numéricos dos k-mers,
   que estão na base 4, para a base 10.O kernel ComputeIndex é responsável por
   esta conversão. O resultado desta fase é armazenado em um vetor denominado
   vetor de conversão.

   Nota: No código abaixo, estão sendo transformados 2 caracteres de Seq em 1
         de Index. Entretando, acredito que não esteja pegando as posições
         corretas. Por exemplo:

         Seq   -> 1 1 3 2 1
         Index -> 5 (11), 7 (13), 14 (32)
*/
void ComputeIndex(char *Seq, int *Index, const int k, lint nN, ushort offset)
{
   lint start = 0;
   lint end = offset;

   for(lint id = start; id < end; id++)
   {
      lint index = 0;

      if(id < nN) {
         // printf("@deb | nuc: ");
         for( lint i = 0; i < k; i++ )
         {
            char nuc = Seq[i + id];

            if (nuc != -1) //Verifica se ha alguem que nao deve ser processado
            {
               //printf("%d", nuc);
               index += nuc * powf(4, ((k-1)-i));
            }
            else
            {
               index = -1;
               break;
            }
         }

         // printf("\n");

         Index[id] = index;// Value of the combination
         // printf("@deb | ComputeIndex | Index[%ld]: %d\n", id, Index[id]);
      }
   }
}

/*
   Compute k-mer frequency
   Na segunda fase é calculada a repetição dos valores convertidos na fase
   anterior. O kernel ComputeFreq é responsável pela contabilização da
   repetição dos k-mers. O resultado do cálculo é armazenado em um vetor
   denominado vetor de frequência.
*/
void ComputeFreq(int *Index, int *Freq, lint *start, int *length, ushort offset, int fourk, lint nS, lint nN)
{
   // int idx = threadIdx.x + (blockDim.x * blockIdx.x);
   int idx = 0;

   if (idx < nS)
   {
      int end = start[idx] + (length[idx] + 1);

      for (int i = start[idx]; i < end; i++)
      {
         if (Index[i] != -1 && i < nN)
         {
            int pos = (fourk * idx) + Index[i];
            Freq[pos] += 1;
            // printf("%d - Freq[%d]: %d\n", i, pos, Freq[pos]);
         }
      }
   }
}

//New way to compute k-mer frequency
void ComputeFreqNew(int *Index, int *Freq, lint *start, int *length, ushort offset, int fourk, lint nS)
{
}
