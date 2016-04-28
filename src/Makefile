NVCCFLAGS=-lm -m64 -arch=sm_35 -O3

default: kmer

kmer:
	nvcc $(NVCCFLAGS) kmer_kernel.cu -lz -c
	nvcc $(NVCCFLAGS) kmer_main.cu -lz -c
	nvcc $(NVCCFLAGS) main.cu -c
	nvcc main.o kmer_kernel.o kmer_main.o -lz -o ../bin/kmer
	rm main.o kmer_main.o kmer_kernel.o
clean:
	rm -f ../bin/* ./*.o
