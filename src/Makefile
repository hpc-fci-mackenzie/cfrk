NVCCFLAGS=-lm -m64 -arch=sm_35 -O3

<<<<<<< HEAD
default: kmer split serial
=======
default: kmer serial
>>>>>>> 4cff1901a4ac984b3789fd9ba6f9e912754af930

kmer:
	nvcc $(NVCCFLAGS) kmer_kernel.cu -lz -c
	nvcc $(NVCCFLAGS) kmer_main.cu -lz -c
	nvcc $(NVCCFLAGS) main.cu -c
	nvcc main.o kmer_kernel.o kmer_main.o -lz -o ../bin/cfrk
	rm main.o kmer_main.o kmer_kernel.o
<<<<<<< HEAD
split:
	gcc -Wall -std=c99 -pedantic -o ../bin/fastasplitn fastasplitn.c
serial:
	nvcc serial.cu -o ../bin/serial -lz	
=======
serial:
	nvcc serial.cu -o ../bin/serial -lz
>>>>>>> 4cff1901a4ac984b3789fd9ba6f9e912754af930
clean:
	rm -f ../bin/* ./*.o
