#NVCCFLAGS=-ccbin=icpc -lm -m64 -arch=sm_35 --use_fast_math --default-stream per-thread --ptxas-options=-O3 
#INTELFLAGS=-Xcompiler="-qopenmp -Ofast -pthread"
NVCCFLAGS=-lm -m64 -arch=sm_50 --use_fast_math --default-stream per-thread --ptxas-options=-O3
INTELFLAGS=-Xcompiler -fopenmp -lpthread
default: cfrk

cfrk:
	nvcc $(INTELFLAGS) $(NVCCFLAGS) kmer_kernel_data_struct.cu -c
	nvcc $(INTELFLAGS) $(NVCCFLAGS) kmer_main_data_struct.cu -c
	nvcc main_data_struct.cu $(INTELFLAGS) $(NVCCFLAGS) -c
	nvcc main_data_struct.o kmer_kernel_data_struct.o kmer_main_data_struct.o $(NVCCFLAGS) $(INTELFLAGS) -o ../bin/cfrk
	rm main_data_struct.o kmer_main_data_struct.o kmer_kernel_data_struct.o
clean:
	rm -f ../bin/* ./*.o
