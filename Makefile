NVIDIA_SDK = $(HOME)/NVIDIA_GPU_Computing_SDK

COMMON_INCLUDES =
COMMON_DEFINES =
COMMON_LIBS = 

CXX = g++
CFLAGS = -O3 $(COMMON_INCLUDES) $(COMMON_DEFINES)

NVCC = nvcc

# uncomment the relevant line for your hardware,
# or leave all uncommented for a generic binary
NVCC_ARCHES = -gencode arch=compute_13,code=sm_13
NVCC_ARCHES += -gencode arch=compute_20,code=sm_20
#NVCC_ARCHES += -gencode arch=compute_20,code=sm_21

NVCC_CFLAGS = -O3 $(COMMON_INCLUDES) $(COMMON_DEFINES) $(NVCC_ARCHES) --compiler-options="$(CFLAGS) -fno-strict-aliasing" -use_fast_math --ptxas-options="-v"

CUDALucas: CUDALucas.o setup.o rw.o balance.o zero.o
	$(CXX) -fPIC -O3 -o CUDALucas CUDALucas.o setup.o rw.o balance.o zero.o $(COMMON_LIBS) -lcudart -lcufft -lm
CUDALucas.o: CUDALucas.cu
	$(NVCC) $(NVCC_CFLAGS) -c CUDALucas.cu
setup.o: setup.c
	$(CXX) $(CFLAGS) -c setup.c
rw.o: rw.c
	$(CXX) $(CFLAGS) -c rw.c
balance.o: balance.c
	$(CXX) $(CFLAGS) -c balance.c
zero.o: zero.c
	$(CXX) $(CFLAGS) -c zero.c

test: CUDALucas
	-rm c216091 t216091
	@echo "Iteration 10000 M( 216091 )C, 0x00000000758f6786 < expected"
	./CUDALucas -o- -c 10000 -t 216091

clean:
	-rm *.o  *~
	-rm CUDALucas
	-rm c216091 t216091
