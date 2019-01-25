#NVCCFLAGS=-ccbin=icpc -lm -m64 -arch=sm_35 --use_fast_math --default-stream per-thread --ptxas-options=-O3 
#INTELFLAGS=-Xcompiler="-qopenmp -Ofast -pthread"
NVCCFLAGS=-lm -m64 -arch=sm_50 --use_fast_math --default-stream per-thread --ptxas-options=-O3
INTELFLAGS=-Xcompiler -fopenmp -lpthread
default: cfrk split serial

cfrk:
	nvcc $(INTELFLAGS) $(NVCCFLAGS) kmer_kernel.cu -lz -c
	nvcc $(INTELFLAGS) $(NVCCFLAGS) kmer_main.cu -lz -c
	nvcc main.cu $(INTELFLAGS) $(NVCCFLAGS) -c
	nvcc main.o kmer_kernel.o kmer_main.o $(NVCCFLAGS) $(INTELFLAGS) -lz -o ../bin/cfrk
	rm main.o kmer_main.o kmer_kernel.o
split:
	gcc -Wall -std=c99 -pedantic -o ../bin/fastasplitn fastasplitn.c
serial:
	nvcc $(NVCCFLAGS) serial.cu -o ../bin/serial -lz
clean:
	rm -f ../bin/* ./*.o
