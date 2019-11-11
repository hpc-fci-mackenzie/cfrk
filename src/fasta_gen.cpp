/*
	Gerador de arquivos FASTA

	Fontes:
	https://zhanglab.ccmb.med.umich.edu/FASTA/
	https://pt.wikipedia.org/wiki/Formato_FASTA
	http://www.acgt.me/blog/2013/6/25/the-fasta-file-format-a-showcase-for-the-best-and-worst-of-b.html
	https://www.uniprot.org/downloads
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SEQUENCE_SIZE 1024
#define NUMBER_OF_SEQUENCES 10

int main() {
   int ssIndex, nofIndex;
   char dna[4] = {'A', 'T', 'C', 'G'};
   FILE *fp;

   srand(time(0));
   fp = fopen("teste.fasta", "w");

   if(fp != NULL)
   {
      for(nofIndex = 0; nofIndex < NUMBER_OF_SEQUENCES; nofIndex++)
      {
         fprintf(fp, ">gi|DNA para teste do algoritmo CFRK\n");

         for(ssIndex = 0; ssIndex < SEQUENCE_SIZE; ssIndex++)
         {
            putc(dna[rand() % 4], fp);

            if(!((ssIndex + 1) % 75))
               putc('\n', fp);
         }

         putc('\n', fp);
      }

      fclose(fp);
   }

   return 0;
}