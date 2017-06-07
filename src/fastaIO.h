#ifndef _fastaIO_h
#define _fastaIO_h
#pragma once

#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <time.h>
#include <cuda.h>
#include "kseq.h"
#include "tipos.h"

KSEQ_INIT(gzFile, gzread)

//-------------------------------------------------------------------------------------------
void ProcessTmpData(struct tmp_data *tdfirst, struct read *rd, lint nN, lint nS, ushort flag)
{
   struct tmp_data *aux = tdfirst;
   lint i, pos = 0, seq = 0;
   if (flag == 0) //CPU
   {
      rd->data = (char*)malloc(sizeof(char)*(nN + nS));
      rd->length = (int*)malloc(sizeof(int)*nS);
      rd->start = (lint*)malloc(sizeof(lint)*nS);
   }
   if (flag == 1) //GPU
   {
      cudaMallocHost((void**)&rd->data, sizeof(char)*(nN + nS));
      cudaMallocHost((void**)&rd->length, sizeof(int)*nS);
      cudaMallocHost((void**)&rd->start, sizeof(lint)*nS);
   }
   rd->start[0] = 0;
   while (aux != NULL && seq < nS)
   {
      for(i = 0; i < aux->length; i++)
      {
         rd->data[pos] = aux->data[i];
         pos++;
      }
      rd->data[pos] = -1;
      pos++;
      rd->length[seq] = aux->length;
      seq++;
      rd->start[seq] = pos;
      aux = aux->next;
   }
}

//-------------------------------------------------------------------------
void ReadFASTASequences(char *file, lint *nN, lint *nS, struct read *rd, ushort flag)
{
   gzFile fp;
   kseq_t *seq;
   struct tmp_data *tdfirst = NULL, *td = NULL, *aux = NULL;
   int len;
   lint lnN = 0, lnS = 0;
   int i;
   //char letter;
   struct timespec start, stop;
   double seconds;

   fp = gzopen(file, "r");
   seq = kseq_init(fp);
   tdfirst = (struct tmp_data*)malloc(sizeof(struct tmp_data));
   td = tdfirst;
   td->next = NULL;
   while ((len = kseq_read(seq)) >= 0)
   {
      clock_gettime(CLOCK_REALTIME, &start);
      lnN += len; //Count the total number of nucleotides read
      lnS++;//count the total number of seqs read
      td->data = (char*)malloc(sizeof(char) * len);
      td->length = len;
      for (i = 0; i < len; i++)
      {   
         //letter = toupper(seq->seq.s[i]);
         switch(seq->seq.s[i])
         {   
            case 'a':
            case 'A':
               td->data[i] = 0; break;
            case 'c':
            case 'C':
               td->data[i] = 1; break;
            case 'g':
            case 'G':
               td->data[i] = 2; break;
            case 't':
            case 'T':
               td->data[i] = 3; break;
            default:
               td->data[i] = -1; break;
         }
      }

      aux = (struct tmp_data*)malloc(sizeof(struct tmp_data));
      td->next = aux;
      aux->next=NULL;
      td = aux;
      if (DBG == 1)
      {
         printf("[util]Seq %ld read!-> ", lnS);
         clock_gettime(CLOCK_REALTIME, &stop);
         seconds = (double)((stop.tv_sec+stop.tv_nsec*1e-9) - (double)(start.tv_sec+start.tv_nsec*1e-9));
         printf("wall time %fs\n", seconds);
      }
   }

   ProcessTmpData(tdfirst, rd, lnN, lnS, flag);

   if (DBG == 1)
     for (i = 0;i < lnS; i++)
     {
       printf("[util] %d\n", rd->length[i]);
     }

   *nN = lnN + lnS;
   *nS = lnS;
   
   gzclose(fp);

}

#endif
