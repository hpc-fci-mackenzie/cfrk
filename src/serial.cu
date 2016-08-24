#include <stdio.h>
#include <stdlib.h>
#include "fastaIO.h"

#define POW(k) (1U << 2*(k))

void ComputeIndex(short *S, int *M, const ushort k, lint nN)
{
   int i, j;
   int index;

   for (i = 0; i < nN; i++)
   {
      index = 0;
      for (j = 0; j < k; j++)
      {
         if (S[j + i] != -1 && i < nN-(k-1)) //Verifica se ha alguem que nao deve ser processado
         {
            index += S[j + i] * POW( (k - 1) - j );
         }
         else
         {
            index = -1;
            break;
         }
      }//End for j
      M[i] = index;
   }//End for i
}

void ComputeFreq(int *M, ushort *Freq, int *start, int *length, ushort fourk, lint nS)
{
   int i, j;

   for (i = 0; i < nS; i++)
   {
      int end = start[i] + (length[i] + 1);
      for (j = start[i]; j < end; j++)
      {
         if (M[j] != -1)
         {
            int pos = (fourk * i) + M[j];
            Freq[pos] = Freq[pos] + 1;
         }
      }//End for j
   }//End for i
}


void Print(int *M, lint nN)
{
   int i;

   puts("Seq");
   for (i = 0; i < nN; i++)
   {
     printf("%d ", M[i]);
   }
}

void PrintFreq(ushort *Freq, lint nS, ushort fourk)
{
   //int cont = 0;
   for (int i = 0; i <= nS*fourk; i++)
   {
/*
      if (i % (fourk + 1) != 0 || i == 0)
      {
         printf("%5d: %d, ", cont, Freq[i]);
         cont++;
      }
      else
      {
         printf("\n");
         cont = 0;
      }
*/
      printf("%d: %d\n", i, Freq[i]);
   }
}

struct read* SelectChunk(struct read *rd, ushort chunkSize, ushort it, lint max, lint gnS, lint *nS, lint gnN, lint *nN)
{
   struct read *chunk;
   int i;
   lint length = 0;

   // Size to be allocated
   for (i = 0; i < max; i++)
   {
      int id = chunkSize*it + i;
      if (id > gnS-1)
      {
         break;
      }
      //printf("id: %d\n", id);
      length += rd->length[id]+1;
   }

   cudaMallocHost((void**)&chunk, sizeof(struct read));
   cudaMallocHost((void**)&chunk->data, sizeof(short)*length);
   cudaMallocHost((void**)&chunk->length, sizeof(int)*chunkSize);
   cudaMallocHost((void**)&chunk->start, sizeof(int)*chunkSize);

   // Copy rd->data to chunk->data
   lint start = rd->start[chunkSize*it];
   lint end = start + length;
   for (i = start; i < end; i++)
   {
      chunk->data[i-start] = rd->data[i];
   }

   chunk->length[0] = rd->length[chunkSize*it];
   chunk->start[0] = 0;
   // Copy start and length
   for (i = 1; i < max; i++)
   {
      int id = chunkSize*it + i;
      chunk->length[i] = rd->length[id];
      chunk->start[i] = chunk->start[i-1]+(chunk->length[i-1]+1);
   }

   *nN = length;
   *nS = max;
return chunk;
}

int main(int argc, char *argv[])
{

   long int gnN, gnS, nS, nN;
   struct read *rd;
   int *Index;
   ushort *Freq;;
   size_t size1, size2;

   if( argc < 2 )
   {
      printf("Usage: %s << <dataset> <k> ", argv[0]);
      return 1;
   }

   int k = atoi( argv[2] ); if( k <= 0 ) return 2;
   int fourk = POW(k);

   puts("\n\t\tReading seqs!!!");
   rd = (struct read*)malloc(sizeof(struct read));
   ReadFASTASequences(argv[1], &gnN, &gnS, rd, 0);
   printf("\nnS: %ld, nN: %ld\n", gnS, gnN);


   int chunkSize = 4096;
   int nChunk = floor(gnS/chunkSize);
   struct read *chunk;
   int i = 0;
   for (i = 0; i < nChunk; i++)
   {
      chunk = SelectChunk(rd, chunkSize, i, chunkSize, gnS, &nS, gnN, &nN);
      size1 = nN * sizeof(int);// Index size
      size2 = pow(4,k) * nS * sizeof(ushort);// Freq size
      Index = (int*)malloc(size1);
      Freq = (ushort*)malloc(size2);

      memset(Index, -1, size1);

      ComputeIndex(chunk->data, Index, k, nN);
      ComputeFreq(Index, Freq, chunk->start, chunk->length, fourk, nS);
      cudaFree(chunk->data);
      cudaFree(chunk->length);
      cudaFree(chunk->start);
      cudaFree(chunk);
   }
   int chunkRemain = abs(gnS - (nChunk*chunkSize));
   chunk = SelectChunk(rd, chunkSize, nChunk, chunkRemain, gnS, &nS, gnN, &nN);
   ComputeIndex(chunk->data, Index, k, nN);
   ComputeFreq(Index, Freq, chunk->start, chunk->length, fourk, nS);

   //Print(Index, nN);

   ushort fourk1 = pow(4,k);
   //ushort fourk2 = pow(4,k+1);
   //ComputeFreq(M, Freq2, rd->start, rd->length, fourk2, nS);


   //printf("Freq1:\n");
   //PrintFreq(Freq, nS, fourk);

   printf("\nfim\n");

return 0;
}
