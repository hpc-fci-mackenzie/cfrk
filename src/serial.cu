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

int main(int argc, char *argv[])
{

   long int nN, nS;
   struct read *rd;
   int *M;
   ushort *Freq1;//, *Freq2;
   //double **M_Freq;
   size_t size1, size2;//, size3;

   if( argc < 2 )
   {
      printf("Usage: %s << <dataset> <k> ", argv[0]);
      return 1;
   }

   int k = atoi( argv[2] ); if( k <= 0 ) return 2;

   puts("\n\t\tReading seqs!!!");
   rd = (struct read*)malloc(sizeof(struct read));
   ReadFASTASequences(argv[1], &nN, &nS, rd, 0);
   printf("\nnS: %ld, nN: %ld\n", nS, nN);

   size1 = nN * sizeof(int);// M size
   size2 = pow(4,k) * nS * sizeof(ushort);// Freq size
   //size3 = pow(4,k+1) * nS * sizeof(ushort);// Freq size with k+1
   
   M = (int*)malloc(size1);
   Freq1 = (ushort*)malloc(size2);
   //Freq2 = (ushort*)malloc(size3);

   memset(M, -1, size1);

   ComputeIndex(rd->data, M, k, nN);

   //Print(M, nN);

   ushort fourk1 = pow(4,k);
   ComputeFreq(M, Freq1, rd->start, rd->length, fourk1, nS);
   //ushort fourk2 = pow(4,k+1);
   //ComputeFreq(M, Freq2, rd->start, rd->length, fourk2, nS);

/*
   printf("Freq1:\n");
   PrintFreq(Freq1, nS, fourk1);
   printf("Freq2:\n");
   PrintFreq(Freq2, nS, fourk2);
*/
   printf("\nfim\n");

return 0;
}
