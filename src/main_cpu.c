/*
   CFRK - Contabilizador da Frequência de Repetição de K-mer (Versão CPU)
   Instiuição: Universidade Presbiteriana Mackenzie
   
   Este código foi baseado no código original CFRK-MT desenvolvido por
   Fabrício Gomes Vilasboas (LNCC | Laboratório Nacional de Computação Científica)
*/
#define NUM_THREADS 1

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

void PrintFreqCSV(struct seq *seq, struct read *pchunk, int nChunk, lint chunkSize)
{
    FILE *out;
    int cont = 0;
    int cont_seq = 0;
    char *str = (char *) calloc(32, sizeof(char));
    lint fourk = pow(4, k);

    out = fopen(file_out, "w");
    char *index = (char *) malloc(sizeof(char) * k);
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
                    for (int c = 0; c < k; c++)
                    {
                        index[c] = pchunk[j].counter[i].index[w][c] + 48;
                    }
                    fprintf(out, "%d,%d,%s,%d\n", j, i, index, pchunk[j].counter[i].frequency[w]);
                }
            }
        }
    }
    fclose(out);
}

void PrintFreq(struct seq *seq, struct read *pchunk, int nChunk, lint chunkSize, struct read *rchunk, lint rChunkSize)
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
            int *line_index = (int *) malloc(fourk * sizeof(int));
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
                    int index = 0;
                    for (int c = 0; c < k; c++)
                    {
                        index += pchunk[j].counter[i].index[w][c] * pow(4, (k - 1) - c);
                    }
                    line_index[index] = index;
                    line_frequency[index] = pchunk[j].counter[i].frequency[w];
                }
            }
            char *str = (char *) malloc(32 * sizeof(char));
            #pragma omp parallel for
            for (int z = 0; z < 32; z++)
            {
                str[z] = 0;
            }
            for (int w = 0; w < fourk; w++)
            {
                if (line_index[w] > -1)
                {
                    fprintf(out, "%d:%d ", line_index[w], line_frequency[w]);
                }
                else
                {
                    fprintf(out, "%d:%d ", w, 0);
                }
            }
            fprintf(out, "\t\t\t\t\n");
            free(str);
        }
    }
    // Remain sequences
    if (rChunkSize > 0)
    {
        for (lint i = 0; i < (rChunkSize); i++)
        {
            int n_combination = (rchunk->length[i] - k + 1);
            int *line_index = (int *) malloc(fourk * sizeof(int));
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
                    int index = 0;
                    for (int c = 0; c < k; c++)
                    {
                        index += rchunk->counter[i].index[w][c] * pow(4, (k - 1) - c);
                    }
                    line_index[index] = index;
                    line_frequency[index] = rchunk->counter[i].frequency[w];
                }
            }
            char *str = (char *) malloc(32 * sizeof(char));
            #pragma omp parallel for
            for (int z = 0; z < 32; z++)
            {
                str[z] = 0;
            }
            for (int w = 0; w < fourk; w++)
            {
                if (line_index[w] > -1)
                {
                    fprintf(out, "%d:%d ", line_index[w], line_frequency[w]);
                }
                else
                {
                    fprintf(out, "%d:%d ", w, 0);
                }
            }
            fprintf(out, "\t\t\t\t\n");
            free(str);
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
        fprintf(out, "%d,%d\n", x, remaining_chunk_size + (nChunk * each_chunk_size));
    }
    fclose(out);
}

int main(int argc, char *argv[])
{
    lint gnN, gnS, chunkSize = 8192;
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
    double st = time(NULL);
//    puts("\t... Reading sequences ...");
    rd = (struct read *) malloc( sizeof(struct read));
    struct seq *seq = ReadFASTASequences(argv[1], &gnN, &gnS, rd, 1);

    double et = time(NULL);
//    printf("> Reading time: %1f\n", (et - st));
    printf("%d,%d,%.4f,", k, nt, (et - st));
//    printf("> nS: %ld, nN: %ld\n", gnS, gnN);
//    for (int x = 1; x <= 300; x++)
//    {
//        memoryUsage(x, chunkSize, gnS);
//    }
    /*
       (!) DIVIDINDO OS DADOS EM CHUNKS
    */
    int nChunk = floor(gnS / chunkSize);
//    printf("> nChunk: %d\n", nChunk);
    int i = 0;

    chunk = (struct read *) malloc(nChunk * sizeof(struct read));
    nS = (lint *) malloc(nChunk * sizeof(lint));
    nN = (lint *) malloc(nChunk * sizeof(lint));

    st = time(NULL);
    SelectChunk(chunk, nChunk, rd, chunkSize, chunkSize, gnS, nS, gnN, nN);
    et = time(NULL);

//    printf("> Chunk Construction Time: %1f\n", (et - st));
    printf("%.4f,", (et - st));

    offset = floor(nChunk / devCount);
//    printf("> offset: %d\n", offset);


    st = time(NULL);
//    #pragma omp parallel for
    for (i = 0; i < nChunk; i++)
    {
//        printf("\tkmer_main - chunk(%d)\n", i);
        kmer_main(&chunk[i], nN[i], nS[i], k);
    }
    et = time(NULL);
//    printf("> Processing Time: %1f\n", (et - st));
    printf("%.4f,", (et - st));

    int chunkRemain = abs(gnS - (nChunk * chunkSize));
    lint rnS, rnN;
//    printf("> chunkRemain: %d\n", chunkRemain);
    struct read *chunk_remain;

    if (chunkRemain)
    {
        chunk_remain = SelectChunkRemain(rd, chunkSize, nChunk, chunkRemain, gnS, &rnS, gnN, &rnN);
        st = time(NULL);
        kmer_main(chunk_remain, rnN, rnS, k);
        et = time(NULL);
//        printf("> Remain Processing Time: %1f\n", (et - st));
        printf("%.4f,", (et - st));
    }
    else
    {
        printf("%.4f", 0.0);
    }

    st = time(NULL);
//    PrintFreq(seq, chunk, nChunk, chunkSize, chunk_remain, rnS);
    et = time(NULL);
//    printf("> Writing time: %1f\n", (et - st));
    printf("%.4f\n", (et - st));
    free(rd);
    free(nS);
    free(nN);
    free(chunk);
    free(chunk_remain);
    free(seq);
    return 0;
}
