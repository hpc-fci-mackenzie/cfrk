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
   ssize_t read, oldRead;
   struct seq *seq;
   int count = -1, flag = 0; 

   *nS = GetNs(fileName);
   seq = (struct seq*)malloc(sizeof(struct seq) * *nS);
       
   if ((fastaFile = fopen(fileName, "r")) == NULL) exit(EXIT_FAILURE);
          
   while ((read = getline(&line, &len, fastaFile)) != -1)
   {      
       if (line[0] == '>')
       {
          count++;
          seq[count].header = (char*)malloc(sizeof(char)*read);
          strcpy(seq[count].header, line);
          flag = 0;
       }     
       else  
       {  
          if (flag == 0)
          {
             seq[count].data = (char*)malloc(sizeof(char)*read);
             strcat(seq[count].data, line);
             flag = 1;
          }
          else
          {
             oldRead = strlen(seq[count].data);
             aux = (char*)malloc(sizeof(char)*oldRead);
             strcpy(aux, seq[count].data);
             seq[count].data = NULL;
             seq[count].data = (char*)malloc(sizeof(char)*(read+oldRead));
             strcat(seq[count].data, aux);
             strcat(seq[count].data, line);
             aux = NULL;
          }
       }
   }
   return seq;
}

//-------------------------------------------------------------------------------------------
void ProcessData(struct seq *seq, struct read *rd, lint nN, lint nS, ushort flag)
{
   lint i, j, pos = 0, seqCount = 0;

   cudaMallocHost((void**)&rd->data, sizeof(char)*(nN + nS));
   cudaMallocHost((void**)&rd->length, sizeof(int)*nS);
   cudaMallocHost((void**)&rd->start, sizeof(lint)*nS);

   rd->start[0] = 0;

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
   }
}

//-------------------------------------------------------------------------
void ReadFASTASequences(char *file, lint *nN, lint *nS, struct read *rd, ushort flag)
{
   struct seq *seq;
   int len;
   lint lnN = 0;
   int i, j;

   seq = ReadFasta(file, nS);
   for (i = 0; i < *nS; i++)
   {
      len = strlen(seq[i].data);
      lnN += len;
      for (j = 0; j < len; j++)
      {   
         //char letter = toupper(seq->seq.s[i]);
         switch(seq->data[j])
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
   
}

#endif
