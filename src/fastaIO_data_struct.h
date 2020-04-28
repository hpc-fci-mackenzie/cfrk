#ifndef _fastaIO_h
#define _fastaIO_h
#pragma once

#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <time.h>
#include <cuda.h>
#include "tipos.h"

int GetNs(char *FileName)
/*Get the length of the sequence
 * */
{
   char temp[64], cmd[512];
   FILE *in;
   strcpy(cmd, "grep -c \">\" ");
   strcat(cmd, FileName);
   in = popen(cmd, "r");
   fgets(temp, 64, in);
   fclose(in);
   return atoi(temp);
}

struct seq *ReadFasta(char *fileName, lint *n_sequence)
/*
 * Create sequence struct from file information
 * */
{
   FILE *fastaFile;
   char *line = NULL, *aux;
   size_t len = 0;
   ssize_t size, oldRead;
   struct seq *seq;
   int count = -1, flag = 0; 

   *n_sequence = GetNs(fileName); // Gets the sequence as Integer
   seq = (struct seq*)malloc(sizeof(struct seq) * *n_sequence); // Allocate sequence memory on Main Memory
       
   if ((fastaFile = fopen(fileName, "r")) == NULL) exit(EXIT_FAILURE);
          
   while ((size = getline(&line, &len, fastaFile)) != -1)
   {
      if (line[0] == '>')
      {
         count++;
         seq[count].header = (char*)malloc(sizeof(char)*len);
         strcpy(seq[count].header, line);
         flag = 0;
      }
      else
      {
         if (flag == 0)
         {
            seq[count].read = (char*)malloc(sizeof(char)*len);
            strcat(seq[count].read, line);
	         seq[count].len = strlen(seq[count].read) - 1;
            flag = 1;
         }
         else
         {
            oldRead = strlen(seq[count].read);
            aux = (char*)malloc(sizeof(char)*oldRead);
            strcpy(aux, seq[count].read);
            seq[count].read = NULL;
            seq[count].read = (char*)malloc(sizeof(char)*(len+oldRead));
            strcat(seq[count].read, aux);
            strcat(seq[count].read, line);
	         seq[count].len = strlen(seq[count].read) - 1; // len variable have no use since the code uses strlen function to get Length information
            aux = NULL;
         }
      }
   }
   return seq;
}

//-------------------------------------------------------------------------------------------
void ProcessData(struct seq *seq, struct chunk *chunk, lint n_concat_sequence_length, lint n_sequence, ushort flag, int k)
{
   lint i, j, pos = 0, seqCount = 0, w;

   cudaMallocHost((void**)&chunk->data, sizeof(char)*(n_concat_sequence_length + n_sequence));
   cudaMallocHost((void**)&chunk->length, sizeof(int)*n_sequence);
   cudaMallocHost((void**)&chunk->start, sizeof(lint)*n_sequence);
   cudaMallocHost((void**)&chunk->n_combination, sizeof(int));
   cudaMallocHost((void**)&chunk->counter, sizeof(struct counter)*n_sequence);

   //rd->data = (char*)malloc(sizeof(char)*(n_concat_sequence_length + n_sequence));
   //rd->length = (int*)malloc(sizeof(int)*n_sequence);
   //rd->start = (lint*)malloc(sizeof(lint)*n_sequence);


   chunk->start[0] = 0;

   for (j = 0; j < n_sequence; j++)
   {
      for(i = 0; i < seq[j].len; i++)
      {
         chunk->data[pos] = seq[j].data[i];
         pos++;
      }
      chunk->data[pos] = -1;
      pos++;
      chunk->length[seqCount] = seq[j].len;
      // Initializing Counters X Read
      if ((chunk->length[seqCount] - k +1) > *chunk->n_combination)
      {
         *chunk->n_combination = chunk->length[seqCount] - k +1;
      }
      cudaMallocHost((void**)&chunk->counter[j].index, sizeof(int)**chunk->n_combination);
      cudaMallocHost((void**)&chunk->counter[j].frequency, sizeof(int)**chunk->n_combination);
      for (w = 0; w < *chunk->n_combination; w++)
      {
         chunk->counter[j].index[w] = -1;
         chunk->counter[j].frequency[w] = 0;
      }
      seqCount++;
      chunk->start[seqCount] = pos;
   }
}

//-------------------------------------------------------------------------
struct seq *ReadFASTASequences(char *file, lint *n_concat_sequence_length, lint *n_sequence, struct chunk *chunk, ushort flag, int k)
{
   struct seq *seq;
   int len;
   lint n_concat_sequence_length_partial = 0; // total length of all sequence concatenated
   int i, j;

   seq = ReadFasta(file, n_sequence);

   for (i = 0; i < *n_sequence; i++)
   {
      len = seq[i].len;
      n_concat_sequence_length_partial += len;

      seq[i].data = (char*)malloc(sizeof(char)*len);

      for (j = 0; j < len; j++)
      {
         switch(seq[i].read[j])
         {   
            case 'a':
            case 'A':
               seq[i].data[j] = 0; break;
            case 'c':
            case 'C':
               seq[i].data[j] = 1; break;
            case 'g':
            case 'G':
               seq[i].data[j] = 2; break;
            case 't':
            case 'T':
               seq[i].data[j] = 3; break;
            default:
               seq[i].data[j] = -1; break;
         }
      }
   }

   ProcessData(seq, chunk, n_concat_sequence_length_partial, *n_sequence, flag, k);

   *n_concat_sequence_length = n_concat_sequence_length_partial + *n_sequence;

   return seq;
}

#endif
