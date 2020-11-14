#NVCCFLAGS=-ccbin=icpc -lm -m64 -arch=sm_35 --use_fast_math --default-stream per-thread --ptxas-options=-O3 
#INTELFLAGS=-Xcompiler="-qopenmp -Ofast -pthread"

# Caminho dos diretórios
PROJ_DIR=$(patsubst %/src/,%, $(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
SRC_DIR=$(PROJ_DIR)/src
BIN_DIR=$(PROJ_DIR)/bin
CPU_OUTPUT=cfrk-cpu
CPU_PKG_OUTPUT=cfrk-cpu-pkg
GPU_OUTPUT=cfrk-gpu

# Dependências para os TARGETS de GPU
GPU_CSOURCES=$(wildcard $(SRC_DIR)/base/*.cu)
GPU_OBJFILES=$(patsubst $(SRC_DIR)/base/%.cu,$(SRC_DIR)/base/%.o,$(GPU_CSOURCES))

# Dependências para os TARGETS de CPU
CPU_CSOURCES=$(wildcard $(SRC_DIR)/fork/*.c)
CPU_OBJFILES=$(patsubst $(SRC_DIR)/fork/%.c,$(SRC_DIR)/fork/%.o,$(CPU_CSOURCES))

# Dependências para os TARGETS de CPU
CPU_PKG_CSOURCES=$(wildcard $(SRC_DIR)/fork_pkg/*.c)
CPU_PKG_OBJFILES=$(patsubst $(SRC_DIR)/fork_pkg/%.c,$(SRC_DIR)/fork_pkg/%.o,$(CPU_PKG_CSOURCES))

# Compilação em GPU
NVCCFLAGS=-lm -m64 -arch=sm_50 --default-stream per-thread --ptxas-options=-O3
INTELFLAGS=-Xcompiler -fopenmp -lpthread

# Compilação em CPU
CFLAGS=-Wall -DCPU -fopenmp -O3
LFLAGS=-Wall -lpthread -lm -fopenmp


default: cfrk

cfrk: gpu cpu_seq cpu_pkg

gpu: $(BIN_DIR)/$(GPU_OUTPUT)
	rm $(SRC_DIR)/base/*.o

$(BIN_DIR)/$(GPU_OUTPUT): $(GPU_OBJFILES)
	@echo BULDING: $@
	@nvcc $^ $(NVCCFLAGS) $(INTELFLAGS) -o $@

$(SRC_DIR)/base/%.o: $(SRC_DIR)/base/%.cu
	@echo COMPILANDO: $<
	@nvcc -c $< $(INTELFLAGS) $(NVCCFLAGS) -o $@

cpu_seq: $(BIN_DIR)/$(CPU_OUTPUT)
	rm $(SRC_DIR)/fork/*.o

$(BIN_DIR)/$(CPU_OUTPUT): $(CPU_OBJFILES)
	@echo BULDING: $@
	@gcc $^ $(LFLAGS) -o $@

$(SRC_DIR)/fork/%.o: $(SRC_DIR)/fork/%.c
	@echo COMPILANDO: $<
	@gcc -c $< $(CFLAGS) -o $@

cpu_pkg: $(BIN_DIR)/$(CPU_PKG_OUTPUT)
	rm $(SRC_DIR)/fork_pkg/*.o

$(BIN_DIR)/$(CPU_PKG_OUTPUT): $(CPU_PKG_OBJFILES)
	@echo BULDING: $@
	@gcc $^ $(LFLAGS) -o $@

$(SRC_DIR)/fork_pkg/%.o: $(SRC_DIR)/fork_pkg/%.c
	@echo COMPILANDO: $<
	@gcc -c $< $(CFLAGS) -o $@

clean:
	rm -f ../bin/* ./*.o
