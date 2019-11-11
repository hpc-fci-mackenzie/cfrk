#ifndef _fastaIO_h
#define _fastaIO_h
#pragma once

#include <stdio.h>
#include <string.h>
#include <time.h>
#include "tipos.h"

#ifndef CPU
#include <zlib.h>
#include <cuda.h>
#endif

int GetNs(char *FileName)
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

struct seq *ReadFasta(char *fileName, lint *nS)
{
   FILE *fastaFile;
   char *line = NULL, *aux;
   size_t len = 0;
   ssize_t size, oldRead;
   struct seq *seq;
   int count = -1, flag = 0;

   *nS = GetNs(fileName);
   seq = (struct seq*)calloc(*nS, sizeof(struct seq));

   if ((fastaFile = fopen(fileName, "r")) == NULL) exit(EXIT_FAILURE);

   while ((size = getline(&line, &len, fastaFile)) != -1)
   {
      /*
         NOTA: Segundo a especificação padrão de C, getline retorna o tamanho
               do buffer lido sem contar o caractere de terminação nulo ('\0').
               Utilizar funções como strcpy e strcat pode dar problema no
               código abaixo porque está sendo alocando memória sem considerar
               os caracteres de terminação nulo. Além disso, foram corrigidos
               pontos de memory leak.
      */

      // printf("@deb | GETLINE | size: %ld - len: %ld\n", size, len);

      if (line[0] == '>')
      {
         count++;
         // seq[count].header = (char*)malloc(sizeof(char)*size);
         seq[count].header = (char*)calloc(size + 1, sizeof(char));
         strcpy(seq[count].header, line);
         flag = 0;
         // printf("@deb | HEADER |seq[%d].header:\n%s", count, seq[count].header);
      }
      else
      {
         if (flag == 0)
         {
            // seq[count].read = (char*)malloc(sizeof(char)*size);
            seq[count].read = (char*)calloc(size + 1, sizeof(char));
            // strcat(seq[count].read, line);
            strcpy(seq[count].read, line);
            seq[count].len = strlen(seq[count].read) - 1;
            flag = 1;
            // printf("@deb | FLAG0 | seq[%d].read:\n%s", count, seq[count].read);
         }
         else
         {
            oldRead = strlen(seq[count].read);
            // aux = (char*)malloc(sizeof(char)*oldRead);
            aux = (char*)calloc(oldRead + 1, sizeof(char));
            strcpy(aux, seq[count].read);
            // seq[count].read = NULL;
            free(seq[count].read);
            // seq[count].read = (char*)malloc(sizeof(char)*(size+oldRead));
            seq[count].read = (char*)calloc((size + oldRead + 1), sizeof(char));
            // strcat(seq[count].read, aux);
            strcpy(seq[count].read, aux);
            strcat(seq[count].read, line);
            seq[count].len = strlen(seq[count].read) - 1;
            // aux = NULL;
            free(aux);
            // printf("@deb | FLAG1 | seq[%d].read:\n%s", count, seq[count].read);
         }
      }
   }

   return seq;
}

//-------------------------------------------------------------------------------------------
void ProcessData(struct seq *seq, struct read *rd, lint nN, lint nS, ushort flag)
{
   lint i, j, pos = 0, seqCount = 0;

#ifndef CPU
   cudaMallocHost((void**)&rd->data, sizeof(char)*(nN + nS));
   cudaMallocHost((void**)&rd->length, sizeof(int)*nS);
   cudaMallocHost((void**)&rd->start, sizeof(lint)*nS);
#else
   rd->data = (char*)calloc((nN + nS), sizeof(char));
   rd->length = (int*)calloc(nS, sizeof(int));
   rd->start = (lint*)calloc(nS, sizeof(lint));
#endif

   // rd->start[0] = 0;

   for (j = 0; j < nS; j++)
   {
      for(i = 0; i < seq[j].len; i++)
      {
         rd->data[pos] = seq[j].data[i];
         pos++;
      }
      rd->data[pos] = -1;
      pos++;
      rd->length[seqCount] = seq[j].len;
      seqCount++;
      rd->start[seqCount] = pos;
      // printf("@deb | ProcessData | rd->start[%ld]: %ld\n", seqCount, rd->start[seqCount]);
   }
}

//-------------------------------------------------------------------------
struct seq *ReadFASTASequences(char *file, lint *nN, lint *nS, struct read *rd, ushort flag)
{
   struct seq *seq;
   int len; // tamanho de cada sequência
   lint lnN = 0; // soma do tamanho de todas as sequências
   int i, j;

   seq = ReadFasta(file, nS);

   for (i = 0; i < *nS; i++)
   {
      len = seq[i].len;
      lnN += len;

      // seq[i].data = (char*)malloc(sizeof(char)*len);
      seq[i].data = (char*)calloc(len, sizeof(char));

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

   ProcessData(seq, rd, lnN, *nS, flag);

   *nN = lnN + *nS;

   return seq;
}

#endif
