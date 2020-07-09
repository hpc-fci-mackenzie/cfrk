#NVCCFLAGS=-ccbin=icpc -lm -m64 -arch=sm_35 --use_fast_math --default-stream per-thread --ptxas-options=-O3 
#INTELFLAGS=-Xcompiler="-qopenmp -Ofast -pthread"
NVCCFLAGS=-lm -m64 -arch=sm_50 --default-stream per-thread --ptxas-options=-O3
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

# 
# Compilação C Tradicional
#

# Caminho dos diretórios
PROJ_DIR=$(patsubst %/src/,%, $(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
SRC_DIR=$(PROJ_DIR)/src
BIN_DIR=$(PROJ_DIR)/bin
OUTPUT=cfrk-cpu.out
#OUTPUT=crfk-cpu.exe

# Dependências para os TARGETS
CSOURCES=$(wildcard $(SRC_DIR)/*.c)
OBJFILES=$(patsubst $(SRC_DIR)/%.c,$(SRC_DIR)/%.o,$(CSOURCES))

# Flags de compilação
CFLAGS=-Wall -DCPU -fopenmp
LFLAGS=-Wall -lpthread -lm -fopenmp

cpu: $(BIN_DIR)/$(OUTPUT)
	@rm $(SRC_DIR)/*.o

$(BIN_DIR)/$(OUTPUT): $(OBJFILES)
	@echo BULDING: $@
	@gcc $^ $(LFLAGS) -o $@

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	@echo COMPILANDO: $<
	@gcc -c $< $(CFLAGS) -o $@

clean-cpu:
	@rm $(BIN_DIR)/*.out

fasta-gen:
	@g++ $(SRC_DIR)/fasta_gen.cpp -o $(BIN_DIR)/fasta-gen.out