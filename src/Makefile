#NVCCFLAGS=-ccbin=icpc -lm -m64 -arch=sm_35 --use_fast_math --default-stream per-thread --ptxas-options=-O3 
#INTELFLAGS=-Xcompiler="-qopenmp -Ofast -pthread"
NVCCFLAGS=-lm -m64 -arch=sm_50 --default-stream per-thread --ptxas-options=-O3 #--use_fast_math
INTELFLAGS=-Xcompiler -fopenmp -lpthread
default: cfrk

cfrk:
	nvcc $(INTELFLAGS) $(NVCCFLAGS) kmer_kernel.cu -c
	nvcc $(INTELFLAGS) $(NVCCFLAGS) kmer_main.cu -c
	nvcc main.cu $(INTELFLAGS) $(NVCCFLAGS) -c
	nvcc main.o kmer_kernel.o kmer_main.o $(NVCCFLAGS) $(INTELFLAGS) -o ../bin/cfrk
	rm main.o kmer_main.o kmer_kernel.o
clean:
	rm -f ../bin/* ./*.o
