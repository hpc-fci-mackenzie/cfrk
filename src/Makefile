NVCCFLAGS=-lm -m64 -arch=sm_35 -O3

default: kmer split serial

kmer:
	nvcc $(NVCCFLAGS) kmer_kernel.cu -lz -c
	nvcc $(NVCCFLAGS) kmer_main.cu -lz -c
	nvcc $(NVCCFLAGS) main.cu -c
	nvcc main.o kmer_kernel.o kmer_main.o -lz -o ../bin/cfrk
	rm main.o kmer_main.o kmer_kernel.o
split:
	gcc -Wall -std=c99 -pedantic -o ../bin/fastasplitn fastasplitn.c
serial:
	nvcc serial.cu -o ../bin/serial -lz
clean:
	rm -f ../bin/* ./*.o
