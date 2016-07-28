NVCCFLAGS=-lm -m64 -arch=sm_20 -O3

default: kmer serial

kmer:
	nvcc $(NVCCFLAGS) kmer_kernel.cu -lz -c
	nvcc $(NVCCFLAGS) kmer_main.cu -lz -c
	nvcc $(NVCCFLAGS) main.cu -c
	nvcc main.o kmer_kernel.o kmer_main.o -lz -o ../bin/cfrk
	rm main.o kmer_main.o kmer_kernel.o
serial:
	nvcc serial.cu -o ../bin/serial -lz
clean:
	rm -f ../bin/* ./*.o
