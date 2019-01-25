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
             seq[count].seq = (char*)malloc(sizeof(char)*read);
             strcat(seq[count].seq, line);
             flag = 1;
          }
          else
          {
             oldRead = strlen(seq[count].seq);
             aux = (char*)malloc(sizeof(char)*oldRead);
             strcpy(aux, seq[count].seq);
             seq[count].seq = NULL;
             seq[count].seq = (char*)malloc(sizeof(char)*(read+oldRead));
             strcat(seq[count].seq, aux);
             strcat(seq[count].seq, line);
             aux = NULL;
          }
       }
   }
   return seq;
}

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
   struct seq *seq;
   struct tmp_data *tdfirst = NULL, *td = NULL, *aux = NULL;
   int len;
   lint lnN = 0;
   int i, j;

   seq = ReadFasta(file, nS);
   tdfirst = (struct tmp_data*)malloc(sizeof(struct tmp_data));
   td = tdfirst;
   td->next = NULL;
   for (i = 0; i < *nS; i++)
   {
      len = strlen(seq[i].seq);
      lnN += len;
      td->data = (char*)malloc(sizeof(char) * len);
      td->length = len;
      for (j = 0; j < len; j++)
      {   
         //char letter = toupper(seq->seq.s[i]);
         switch(seq->seq[j])
         {   
            case 'a':
            case 'A':
               td->data[j] = 0; break;
            case 'c':
            case 'C':
               td->data[j] = 1; break;
            case 'g':
            case 'G':
               td->data[j] = 2; break;
            case 't':
            case 'T':
               td->data[j] = 3; break;
            default:
               td->data[j] = -1; break;
         }
      }

      aux = (struct tmp_data*)malloc(sizeof(struct tmp_data));
      td->next = aux;
      aux->next=NULL;
      td = aux;
   }

   ProcessTmpData(tdfirst, rd, lnN, *nS, flag);


   *nN = lnN + *nS;
   
}

#endif
