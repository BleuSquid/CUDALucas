NVIDIA_SDK = $(HOME)/NVIDIA_GPU_Computing_SDK

COMMON_INCLUDES = -I/usr/local/include -I$(NVIDIA_SDK)/C/common/inc
COMMON_DEFINES =

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
	$(CXX) -fPIC -O3 -o CUDALucas CUDALucas.o setup.o rw.o balance.o zero.o -L/usr/local/cuda/lib64 -lcudart -L/usr/local/cuda/lib64 -lcufft -lm
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

clean:
	-rm *.o CUDALucas *~
