#!/usr/bin/make -f

# Linux Makefile for CUDALucas

COMMON_INCLUDES =
COMMON_DEFINES =
COMMON_LIBS = 

CXX = g++
#CXX = g++-4.6 -floop-parallelize-all -ftree-parallelize-loops=4 -ftree-loop-distribution -Wall -Wextra -pedantic -march=native
CFLAGS = -O3 $(COMMON_INCLUDES) $(COMMON_DEFINES)

NVCC = nvcc

# uncomment the relevant line for your hardware,
# or leave all uncommented for a generic binary
# sm_21 is actually slower than sm_20 on sm_21 hardware...
NVCC_ARCHES = -gencode arch=compute_13,code=sm_13
NVCC_ARCHES += -gencode arch=compute_20,code=sm_20
NVCC_ARCHES += -gencode arch=compute_20,code=sm_21

NVCC_CFLAGS = -O3 $(COMMON_INCLUDES) $(COMMON_DEFINES) $(NVCC_ARCHES) --compiler-options="$(CFLAGS) -fno-strict-aliasing" -use_fast_math --ptxas-options="-v"

CUDALucas: CUDALucas.o setup.o rw.o balance.o zero.o main.o
	$(CXX) -fPIC -O3 -o CUDALucas CUDALucas.o setup.o rw.o balance.o zero.o main.o $(COMMON_LIBS) -Wl,-O1 -Wl,--as-needed -lcudart -lcufft -lm
CUDALucas.o: CUDALucas.cu cuda_safecalls.h setup.h
	$(NVCC) $(NVCC_CFLAGS) -c CUDALucas.cu
setup.o: setup.c setup.h
	$(CXX) $(CFLAGS) -c setup.c
rw.o: rw.c rw.h setup.h
	$(CXX) $(CFLAGS) -c rw.c
balance.o: balance.c balance.h globals.h setup.h
	$(CXX) $(CFLAGS) -c balance.c
zero.o: zero.c zero.h setup.h
	$(CXX) $(CFLAGS) -c zero.c
main.o: main.c setup.h balance.h rw.h
	$(CXX) $(CFLAGS) -c main.c

test: CUDALucas
	-rm c216091 t216091
	@echo "Iteration 1000, 0xd222ff4c0604633e <== expected first residue at 1000 iterations"
	./CUDALucas -o- -c 1000 -t 216091

clean:
	-rm *.o *~
	-rm CUDALucas
	-rm c216091 t216091
