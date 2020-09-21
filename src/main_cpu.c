/*
   CFRK - Contabilizador da Frequência de Repetição de K-mer (Versão CPU)
   Instiuição: Universidade Presbiteriana Mackenzie
   
   Este código foi baseado no código original CFRK-MT desenvolvido por
   Fabrício Gomes Vilasboas (LNCC | Laboratório Nacional de Computação Científica)
*/
#define NUM_THREADS 1
#define BILLION  1000000000.0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <omp.h>
#include "tipos_cpu.h"
#include "fastaIO_cpu.h"
#include "kmer_cpu.h"

#define FILENAME_LENGTH 512

//*** Global Variables ***
struct read *chunk;
lint *nS, *nN;
int offset;
int device;
int k;
char file_out[FILENAME_LENGTH];

void PrintFreqCSV(struct seq *seq, struct read *pchunk, int nChunk, lint chunkSize, struct read *rchunk, lint rChunkSize)
{
    FILE *out;
    out = fopen(file_out, "w");
    int start = 0;
    int end = offset;
    fprintf(out, "CHUNK,SEQUENCE,INDEX,FREQUENCY\n");
    for (int j = start; j < end; j++)
    {
        for (lint i = 0; i < (chunkSize); i++)
        {
            for (int w = 0; w < (pchunk[j].length[i] - k + 1); w++)
            {
                if (pchunk[j].counter[i].frequency[w] != 0)
                {
                    fprintf(out, "%d,%ld,%ld,%d\n", end, i, pchunk[j].counter[i].kmer[w], pchunk[j].counter[i].frequency[w]);
                }
            }
        }
    }
    for (lint i = 0; i < (rChunkSize); i++)
    {
        for (int w = 0; w < (rchunk->length[i] - k + 1); w++)
        {
            if (rchunk->counter[i].frequency[w] != 0)
            {
                fprintf(out, "%d,%ld,%ld,%d\n", end, i, rchunk->counter[i].kmer[w], rchunk->counter[i].frequency[w]);
            }
        }
    }
    fclose(out);
}

void OrderedPrintFreq(struct seq *seq, struct read *pchunk, int nChunk, lint chunkSize, struct read *rchunk, lint rChunkSize)
{
    FILE *out;
    lint fourk = pow(4, k);
    out = fopen(file_out, "w");
    int start = 0;
    int end = nChunk;
    // Same size Chunks
    for (int j = start; j < end; j++)
    {
        for (lint i = 0; i < (chunkSize); i++)
        {
            int n_combination = (pchunk[j].length[i] - k + 1);
            lint *line_index = (lint *) malloc(fourk * sizeof(lint));
            int *line_frequency = (int *) malloc(fourk * sizeof(int));
            #pragma omp parallel for
            for (int t = 0; t < fourk; t++)
            {
                line_index[t] = -1;
            }
            #pragma omp parallel for
            for (int w = 0; w < n_combination; w++)
            {
                if (pchunk[j].counter[i].frequency[w] != 0)
                {
                    line_index[pchunk[j].counter[i].kmer[w]] = pchunk[j].counter[i].kmer[w];
                    line_frequency[pchunk[j].counter[i].kmer[w]] = pchunk[j].counter[i].frequency[w];
                }
            }
            for (int w = 0; w < fourk; w++)
            {
                if (line_index[w] > -1)
                {
                    fprintf(out, "%ld:%d ", line_index[w], line_frequency[w]);
                }
                else
                {
                    fprintf(out, "%ld:%d ", (lint) w, 0);
                }
            }
            fprintf(out, "\t\t\t\t\n");
        }
    }
    // Remain sequences
    if (rChunkSize > 0)
    {
        for (lint i = 0; i < (rChunkSize); i++)
        {
            int n_combination = (rchunk->length[i] - k + 1);
            lint *line_index = (lint *) malloc(fourk * sizeof(lint));
            int *line_frequency = (int *) malloc(fourk * sizeof(int));
            #pragma omp parallel for
            for (int t = 0; t < fourk; t++)
            {
                line_index[t] = -1;
            }
            #pragma omp parallel for
            for (int w = 0; w < n_combination; w++)
            {
                if (rchunk->counter[i].frequency[w] != 0)
                {
                    line_index[rchunk->counter[i].kmer[w]] = rchunk->counter[i].kmer[w];
                    line_frequency[rchunk->counter[i].kmer[w]] = rchunk->counter[i].frequency[w];
                }
            }
            for (int w = 0; w < fourk; w++)
            {
                if (line_index[w] > -1)
                {
                    fprintf(out, "%ld:%d ", line_index[w], line_frequency[w]);
                }
                else
                {
                    fprintf(out, "%ld:%d ", (lint) w, 0);
                }
            }
            fprintf(out, "\t\t\t\t\n");
        }
    }
    fclose(out);
}

void PrintFreq(struct seq *seq, struct read *pchunk, int nChunk, lint chunkSize, struct read *rchunk, lint rChunkSize)
{
    FILE *out;
    out = fopen(file_out, "w");
    int start = 0;
    int end = nChunk;
    // Same size Chunks
    for (int j = start; j < end; j++)
    {
        for (lint i = 0; i < (chunkSize); i++)
        {
            int n_combination = (pchunk[j].length[i] - k + 1);
            for (int w = 0; w < n_combination; w++)
            {
                if (pchunk[j].counter[i].frequency[w] != 0)
                {
                    fprintf(out, "%ld:%d ", pchunk[j].counter[i].kmer[w], pchunk[j].counter[i].frequency[w]);
                }
            }
            fprintf(out, "\t\t\t\t\n");
        }
    }
    // Remain sequences
    if (rChunkSize > 0)
    {
        for (lint i = 0; i < (rChunkSize); i++)
        {
            int n_combination = (rchunk->length[i] - k + 1);
            for (int w = 0; w < n_combination; w++)
            {
                if (rchunk->counter[i].frequency[w] != 0)
                {
                    fprintf(out, "%ld:%d ", rchunk->counter[i].kmer[w], rchunk->counter[i].frequency[w]);
                }
            }
            fprintf(out, "\t\t\t\t\n");
        }
    }
    fclose(out);
}

struct read *SelectChunkRemain(struct read *rd, ushort chunkSize, ushort it, lint max, lint gnS, lint *nS, lint gnN, lint *nN)
{
    struct read *chunk;
    lint i;
    lint j;
    lint length = 0;

    // Size to be allocated
    for (i = 0; i < max; i++)
    {
        lint id = chunkSize * it + i;
        if (id > gnS - 1)
        {
            break;
        }
        length += rd->length[id] + 1;
    }


    chunk = (struct read *) malloc(sizeof(struct read));
    chunk->data = (char *) malloc(sizeof(char) * length);
    chunk->length = (int *) malloc(sizeof(int) * chunkSize);
    chunk->start = (lint *) malloc(sizeof(lint) * chunkSize);

    // Copy rd->data to chunk->data
    lint start = rd->start[chunkSize * it];
    lint end = start + (lint) length;


    for (j = start; j < end; j++)
    {
        chunk->data[j - start] = rd->data[j];
    }

    chunk->length[0] = rd->length[chunkSize * it];
    chunk->start[0] = 0;

    // Copy start and length
    for (i = 1; i < max; i++)
    {
        lint id = chunkSize * it + i;
        chunk->length[i] = rd->length[id];
        chunk->start[i] = chunk->start[i - 1] + (chunk->length[i - 1] + 1);
    }

    *nN = length;
    *nS = max;

    return chunk;
}

void SelectChunk(struct read *chunk, const int nChunk, struct read *rd, ushort chunkSize, lint max, lint gnS, lint *nS, lint gnN, lint *nN)
{
    lint i, j, it;

    for (it = 0; it < nChunk; it++)
    {
        lint length = 0;

        // Size to be allocated
        for (i = 0; i < max; i++)
        {
            lint id = chunkSize * it + i;

            if (id > gnS - 1)
            {
                break;
            }

            length += rd->length[id] + 1;
        }

        chunk[it].data = (char *) malloc(length * sizeof(char));
        chunk[it].length = (int *) malloc(chunkSize * sizeof(int));
        chunk[it].start = (lint *) malloc(chunkSize * sizeof(lint));

        lint start = rd->start[chunkSize * it];
        lint end = start + (lint) length;

        for (j = start; j < end; j++)
        {
            chunk[it].data[j - start] = rd->data[j];
        }

        chunk[it].length[0] = rd->length[chunkSize * it];
        chunk[it].start[0] = 0;

        // Copy start and length
        for (i = 1; i < max; i++)
        {
            lint id = chunkSize * it + i;
            chunk[it].length[i] = rd->length[id];
            chunk[it].start[i] = chunk[it].start[i - 1] + (chunk[it].length[i - 1] + 1);
        }

        nN[it] = length;
        nS[it] = max;
    }
}

void *LaunchKmer(void *threadId)
{
    int *start = (int *) threadId;
    int end = offset;

    for (int i = *start; i < end; i++)
    {
        printf("\tkmer_main - chunk(%d)\n", i);
        kmer_main(&chunk[i], nN[i], nS[i], k);
    }

    return NULL;
}

void memoryUsage(int k, int chunkSize, int n_sequence)
{
    FILE *out;
    out = fopen("memory.csv", "w");
    fprintf(out, "K,Memory\n");
    for (int x = 1; x <= 300; x++)
    {
        size_t integer = sizeof(int);
        size_t character = sizeof(char);
        size_t long_integer = sizeof(lint);
        int nChunk = floor(n_sequence/chunkSize);
        int remaining_sequences = n_sequence - (nChunk * chunkSize);
        size_t data_vector = (chunkSize) * character;
        size_t length_vector = chunkSize * character;
        size_t start_vector = chunkSize * long_integer;
        size_t counter_vector = 0;
        int sequence_length = 300;
        for (int i = 0; i < chunkSize; i++)
        {
            counter_vector += (sequence_length - x +1) * ((character * k) + character );
        }
        size_t each_chunk_size = data_vector + length_vector + start_vector + counter_vector;
//    printf("Size of a Chunk: %d\n", each_chunk_size);
//    printf("All Chunks(%d): %d\n", nChunk, nChunk * each_chunk_size);
        data_vector = (remaining_sequences) * character;
        length_vector = remaining_sequences * character;
        start_vector = remaining_sequences * long_integer;
        counter_vector = 0;
        for (int i = 0; i < remaining_sequences; i++)
        {
            counter_vector += (sequence_length - x +1) * ((character * k) + character);
        }
        size_t remaining_chunk_size = data_vector + length_vector + start_vector + counter_vector;
//    printf("Remaining Sequences: %d\n", remaining_sequences);
//    printf("Remaining Sequences Size: %d\n", remaining_chunk_size);
        fprintf(out, "%d,%ld\n", x, remaining_chunk_size + (nChunk * each_chunk_size));
    }
    fclose(out);
}

int main(int argc, char *argv[])
{
    lint gnN, gnS, chunkSize = 8192;
    struct timespec start, end;

    int devCount = 1;
    int nt = 4;
    char dataset[FILENAME_LENGTH] = {0};

    /*
       (!) CONFIGURAÇÃO PADRÃO
          dataset: teste.fasta
          outfile: resultado.cfrk
          número de threads: 12
          tamanho do chunksize: 8192
          k: 2
          nt: 12
    */

    memset(file_out, 0, FILENAME_LENGTH);


    if (argc > 1)
    {
        strcpy(dataset, argv[1]);
        strcpy(file_out, argv[2]);
        k = atoi(argv[3]);

        if (argc == 5)
        {
            nt = atoi(argv[4]);
        }

        if (argc == 6)
        {
            chunkSize = atoi(argv[5]);
        }
    }
    else
    {
        strcpy(dataset, "teste.fasta");
        strcpy(file_out, "resultado.cfrk");
        k = 2;
    }

//    printf("Usage:\n> dataset: %s\n> file_out: %s\n> k: %d\n> nt: %d\n> chunkSize: %d\n", dataset, file_out, k, nt, (int) chunkSize);
    omp_set_num_threads(nt);
    /*
       (!) LEITURA DE SEQUÊNCIAS (FASTA)
    */
    struct read *rd;
    clock_gettime(CLOCK_REALTIME, &start);
//    puts("\t... Reading sequences ...");
    rd = (struct read *) malloc( sizeof(struct read));
    struct seq *seq = ReadFASTASequences(argv[1], &gnN, &gnS, rd, 1);

    clock_gettime(CLOCK_REALTIME, &end);
//    printf("> Reading time: %1f\n", (et - st));
    printf("%d,%d,%.10f,", k, nt, (end.tv_sec - start.tv_sec) +
                                 (end.tv_nsec - start.tv_nsec) / BILLION);
//    printf("> nS: %ld, nN: %ld\n", gnS, gnN);
//    for (int x = 1; x <= 300; x++)
//    {
//        memoryUsage(x, chunkSize, gnS);
//    }
    /*
       (!) DIVIDINDO OS DADOS EM CHUNKS
    */
//    int nChunk = floor(gnS / chunkSize);
    int nChunk = 1;
    chunkSize = (lint) ceil((gnS / nChunk));
//    printf("> nChunk: %d\n", nChunk);
//    printf("> chunkSize: %ld\n", chunkSize);
    int i = 0;

    chunk = (struct read *) malloc(nChunk * sizeof(struct read));
    nS = (lint *) malloc(nChunk * sizeof(lint));
    nN = (lint *) malloc(nChunk * sizeof(lint));

    clock_gettime(CLOCK_REALTIME, &start);
    SelectChunk(chunk, nChunk, rd, chunkSize, chunkSize, gnS, nS, gnN, nN);
    clock_gettime(CLOCK_REALTIME, &end);

//    printf("> Chunk Construction Time: %1f\n", (et - st));
    printf("%.10f,", (end.tv_sec - start.tv_sec) +
                    (end.tv_nsec - start.tv_nsec) / BILLION);

    offset = floor(nChunk / devCount);
//    printf("> offset: %d\n", offset);


    clock_gettime(CLOCK_REALTIME, &start);
//    #pragma omp parallel for
    for (i = 0; i < nChunk; i++)
    {
//        printf("\tkmer_main - chunk(%d)\n", i);
        kmer_main(&chunk[i], nN[i], nS[i], k);
    }
    clock_gettime(CLOCK_REALTIME, &end);
//    printf("> Processing Time: %1f\n", (et - st));
    printf("%.10f,", (end.tv_sec - start.tv_sec) +
                    (end.tv_nsec - start.tv_nsec) / BILLION);

    int chunkRemain = abs(gnS - (nChunk * chunkSize));
    lint rnS, rnN;
//    printf("> chunkRemain: %d\n", chunkRemain);
    struct read *chunk_remain;

    if (chunkRemain)
    {
        chunk_remain = SelectChunkRemain(rd, chunkSize, nChunk, chunkRemain, gnS, &rnS, gnN, &rnN);
        clock_gettime(CLOCK_REALTIME, &start);
        kmer_main(chunk_remain, rnN, rnS, k);
        clock_gettime(CLOCK_REALTIME, &end);
//        printf("> Remain Processing Time: %1f\n", (et - st));
        printf("%.10f,", (end.tv_sec - start.tv_sec) +
                        (end.tv_nsec - start.tv_nsec) / BILLION);
    }
    else
    {
        rnS = 0;
        printf("%.10f,", 0.0);
    }

    clock_gettime(CLOCK_REALTIME, &start);
    PrintFreq(seq, chunk, nChunk, chunkSize, chunk_remain, rnS);
    clock_gettime(CLOCK_REALTIME, &end);
//    printf("> Writing time: %1f\n", (et - st));
    printf("%.10f", (end.tv_sec - start.tv_sec) +
                     (end.tv_nsec - start.tv_nsec) / BILLION);
    free(rd);
    free(nS);
    free(nN);
    free(chunk);
    if (chunkRemain) free(chunk_remain);
    free(seq);
    return 0;
}
