#include <stdio.h>
#include <cuda.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "tipos_data_struct.h"
#include "kmer_data_struct.cuh"


void GetDeviceProp(uint8_t device, lint *maxGridSize, lint *maxThreadDim, lint *deviceMemory) {
    cudaDeviceProp prop;

    cudaGetDeviceProperties(&prop, device);

    *maxThreadDim = prop.maxThreadsDim[0];
    *maxGridSize = prop.maxGridSize[0];
    *deviceMemory = prop.totalGlobalMem;
}

void kmer_main(struct chunk *rd, lint n_concat_sequence_length, lint n_sequence, int k, ushort device) {
    lint block[3], grid[3];// Grid config; 0:n_concat_sequence_length, 1:n_sequence
    lint maxGridSize, maxThreadDim, deviceMemory;// Device config
    ushort offset[3] = {1, 1, 1};
    size_t size[7], totalsize;
    int sizeOfAllCounters = 0; //
    lint i;
    int n_combination = *(rd->n_combination);
    cudaSetDevice(device);
    GetDeviceProp(device, &maxGridSize, &maxThreadDim, &deviceMemory);
    fprintf(stderr, "Combination: %d\n", n_combination);
    fprintf(stderr, "N Sequence: %d\n", n_sequence);


//---------------------------------------------------------------------------

    size[0] = n_concat_sequence_length * sizeof(char);// d_Seq and Seq size
    size[1] = n_sequence * sizeof(int);  // d_length
    size[2] = n_sequence * sizeof(lint); // d_start
    size[3] = n_sequence * sizeof(struct counter); // d_reads
    size[4] = n_combination * sizeof(int);
    size[5] = n_combination * sizeof(int);
    totalsize = size[0] + (size[1] * 2) + size[2] + size[3] + (n_sequence * (sizeof(struct counter) + size[4] + size[5]));
//    fprintf(stderr, "Size[0]: %d\n", size[0]);
//    fprintf(stderr, "Size[1]: %d\n", size[1] * 2);
//    fprintf(stderr, "Size[2]: %d\n", size[2]);
//    fprintf(stderr, "Size[3]: %d\n", size[3]);
//    fprintf(stderr, "Size[4]: %d\n", size[4]);
//    fprintf(stderr, "Size[5]: %d\n", size[5]);
    fprintf(stderr, "TotalSize: %d\n", totalsize);

    if (totalsize > deviceMemory) {
        printf("\n\n\t\t\t[Error] There is no enough space on GPU memory\n");
        printf("\t\t\t[Error] Required memory: %ld; Available memory: %ld\n", totalsize, deviceMemory);
        exit(1);
    }

//---------------------------------------------------------------------------

    fprintf(stderr, "Memory Allocation\n");
    if (cudaMallocManaged(&rd->data, size[0]) != cudaSuccess)
        fprintf(stderr, "\n[Error 1] %s\n", cudaGetErrorString(cudaGetLastError()));
    cudaDeviceSynchronize();

    if (cudaMallocManaged(&rd->length, size[1]) != cudaSuccess)
        fprintf(stderr, "\n[Error 2] %s\n", cudaGetErrorString(cudaGetLastError()));
    cudaDeviceSynchronize();

    if (cudaMallocManaged(&rd->start, size[2]) != cudaSuccess)
        fprintf(stderr, "\n[Error 3] %s\n", cudaGetErrorString(cudaGetLastError()));
    cudaDeviceSynchronize();

    if (cudaMallocManaged(&rd->counter, size[3]) != cudaSuccess)
        fprintf(stderr, "\n[Error 4] %s\n", cudaGetErrorString(cudaGetLastError()));
    cudaDeviceSynchronize();

    for (i = 0; i < n_sequence; i++) {

        if (cudaMallocManaged(&(rd->counter[i].index), size[4]) != cudaSuccess)
            fprintf(stderr, "\n[Error 6-%d] %s\t", i, cudaGetErrorString(cudaGetLastError()));
        cudaDeviceSynchronize();
//        for (int w = 0; w < n_combination; w++) {
//            if (cudaMallocManaged(&(rd->counter[i].index[w]), size[6]) != cudaSuccess)
//                fprintf(stderr, "\n[Error 5-%d-%d] %s\t", i, w, cudaGetErrorString(cudaGetLastError()));
//            cudaDeviceSynchronize();
//        }
        if (cudaMallocManaged(&(rd->counter[i].frequency), size[5]) != cudaSuccess)
            fprintf(stderr, "[Error 7-%d] %s", i, cudaGetErrorString(cudaGetLastError()));

        cudaDeviceSynchronize();
    }

// ************************************************
    fprintf(stderr, "Memory Mapping\n");
    // Thread mapping for raw data
    block[0] = maxThreadDim;
    grid[0] = floor(n_concat_sequence_length / block[0]) + 1;
    if (grid[0] > maxGridSize)
    {
        grid[0] = maxGridSize;
        offset[0] = (n_concat_sequence_length / (grid[0] * block[0])) + 1;
    }

    // Thread mapping for
    block[1] = maxThreadDim;
    grid[1] = (n_combination / block[1]) + 1;
    if (grid[1] > maxGridSize)
    {
        grid[1] = maxGridSize;
        offset[1] = (n_combination / (grid[1] * block[1])) + 1;
    }

    // Thread mapping for
    block[2] = 1;
    grid[2] =  n_sequence;
    if (grid[2] > maxGridSize)
    {
        grid[2] = maxGridSize;
        offset[2] = (n_sequence / (grid[2] * block[2])) + 1;
    }

// ************************************************
    fprintf(stderr, "Kernel Execution Matrix\n");
    for (i = 0; i < n_sequence; i++) {
        SetMatrix<<<grid[1], block[1]>>>(rd->counter, offset[1], n_combination, i);
    }
    fprintf(stderr, "Kernel Execution Compute\n");
    cudaDeviceSynchronize();
//    cudaMemcpyToSymbol("sd_counter", rd->counter, (n_sequence * (sizeof(struct counter) + size[4] + size[5])), size_t(0), cudaMemcpyHostToDevice);
//    cudaDeviceSynchronize();
    ComputeFrequency<<<grid[2], block[2]>>>(rd->data, rd->counter, rd->start, rd->length, k, n_concat_sequence_length,
                                            offset[2], n_sequence, n_combination);
    cudaDeviceSynchronize();
//    cudaMemcpyFromSymbol(rd->counter, "sd_counter", (n_sequence * (sizeof(struct counter) + size[4] + size[5])), size_t(0), cudaMemcpyDeviceToHost);

    for (lint t = 0; t < n_sequence; t ++)
    {
        for (int q = 0; q < n_combination; q++)
        {
           if (rd->counter[t].frequency[q] > 0)
                printf("Sequence: %d Index: %d Frequency: %d\n", t, rd->counter[t].index[q], rd->counter[t].frequency[q]);
        }
    }
//---------------------------------------------------------------------------

    //printf("\nFim kmer_main\n");
}
