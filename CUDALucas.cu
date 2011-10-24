/* CUDALucas.cu
Shoichiro Yamada Oct. 2010 */

/* The following comment is by Guillermo Ballester Valor.  Feel free to send me email
about problems as well.  --Will Edgington, wedgingt@acm.org */

/* MacLucasFFTW.c

lucaslfftw v.1.0.0
Guillermo Ballester Valor, Oct. 1999
gbv@ctv.es

This is an adaptation of Richard Crandall lucdwt.c and Sweeney 
MacLucasUNIX.c code.  There are few things mine own.

The memory requirements is about q bytes, where q is the exponent of
mersenne number being tested. (m(q))
*/

#include <cuda.h>
#include "cuda_safecalls.h"
#include "setup.h"

/************************ definitions ************************************/

/* global variables needed */
double     *g_ttp,*g_ttmp;
float          *g_inv;
double *g_x;
double *g_maxerr;
int *g_carry;

/* rint is not ANSI compatible, so we need a definition for 
* WIN32 and other platforms with rint.
* Also we use that to write the trick to rint()
*/

#if defined(__x86_32__)
#define RINT(x) (floor(x+0.5))	
#else
#define RINT(x) (((x) + 6755399441055744.0 ) - 6755399441055744.0)
#endif

/*
http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html

base code is Mr. ooura's FFT.

Fast Fourier/Cosine/Sine Transform
dimension   :one
data length :power of 2
decimation  :frequency
radix       :4, 2
data        :inplace
table       :use
functions
rdft: Real Discrete Fourier Transform
Appendix :
The cos/sin table is recalculated when the larger table required.
w[] and ip[] are compatible with all routines.
*/

__global__ void rftfsub_kernel(const int n, double *a) {
	const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	double wkr, wki, xr, xi, yr, yi,cc,d,aj,aj1,ak,ak1, *c ;
	
	if(threadID != 0) {
		const int j = threadID * 2;
		c = &a[n+n/4+512*512];
		const int nc = n >> 2 ;
		
		aj  = a[j];
		ak  = a[(n-j)];
		
		wkr = c[nc-threadID];
		wki = c[threadID];
		
		aj1 = a[1+j];
		ak1 = a[1+(n-j)];
		
		xr  = aj - ak;
		wkr = 0.5 - wkr;
		xi  = aj1 + ak1;
		yr  = wkr * xr;
		yi  = wki * xr;
		yr -= wki * xi;
		yi += wkr * xi;
		
		cc  = yr  - aj;
		d   = yi  - aj1;
		aj1 = 2.0 * cc;
		aj  = (cc+d)*(cc-d);
		aj1 = aj1 * d;
		
		cc  = ak  + yr;
		d   = ak1 - yi;
		ak1 = 2.0 * cc;
		ak  = (cc+d)*(cc-d);
		ak1 = ak1 * d;
		xr  = aj  - ak;
		xi  = aj1 + ak1;
		
		yr = wkr * xr;
		yi = -wki * xr;
		yr+= wki * xi;
		yi+= wkr * xi;
		
		aj  = aj - yr;
		ak  = ak + yr;
		aj1 = yi - aj1;
		ak1 = yi - ak1;
		
		a[j]=aj;
		a[(n-j)]=ak;
		a[1+j]=aj1;
		a[1+(n-j)]=ak1;
	} else {
		const int m = n >> 1;
		aj  = a[0];
		ak  = a[1];
		cc  = a[0+m];
		d   = -a[1+m];
		xi  = aj - ak;
		aj += ak;
		ak  = xi;
		aj *= aj;
		
		if ((n & 1) == 0) ak *= ak;
		
		a[1+m] = -2.0*cc;
		ak  = 0.5 * (aj - ak);
		a[1+m] *= d;
		aj -= ak;
		ak  = -ak;
		a[1+m] = -a[1+m];
		a[0+m] = (cc+d)*(cc-d);
		a[0] = aj;
		a[1] = ak;
	}
}

void rftfsub(const int n, double *g_x) {
	dim3 grid(n/512,1,1);
	dim3 threads(128,1,1);
	rftfsub_kernel<<<grid,threads>>>(n,g_x);
}

#define BLOCK_DIM 16

// This kernel is optimized to ensure all global reads and writes are coalesced,
// and to avoid bank conflicts in shared memory.  This kernel is up to 11x faster
// than the naive kernel below.  Note that the shared memory array is sized to 
// (BLOCK_DIM+1)*BLOCK_DIM.  This pads each row of the 2D block in shared memory 
// so that bank conflicts do not occur when threads address the array column-wise.

// This is templated because the transpose is square. We can eliminate the branching
// from the if statements because we know height and width are identical
template <int width>
__global__ void square_transpose(double* __restrict__ odata, const double* __restrict__ idata) {
	__shared__ double block[BLOCK_DIM][BLOCK_DIM+1];

	// read the matrix tile into shared memory	
	unsigned int index_in = (blockIdx.y * (BLOCK_DIM * width)) + (threadIdx.y * width) + (blockIdx.x * BLOCK_DIM) + threadIdx.x;
	block[threadIdx.y][threadIdx.x] = idata[index_in];
	
	__syncthreads();
	
	// write the transposed matrix tile to global memory
	unsigned int index_out = (blockIdx.x * (BLOCK_DIM * width)) + (threadIdx.y * width) + (blockIdx.y * BLOCK_DIM) + threadIdx.x;
	odata[index_out] = block[threadIdx.x][threadIdx.y];
}

template <int width>
__global__ void square_transposef(float* __restrict__ odata, const double* __restrict__ idata) {
	__shared__ float block[BLOCK_DIM][BLOCK_DIM+1];
	
	// read the matrix tile into shared memory
	unsigned int xIndex = blockIdx.x * BLOCK_DIM + threadIdx.x;
	unsigned int yIndex = blockIdx.y * BLOCK_DIM + threadIdx.y;
	
	unsigned int index_in = yIndex * width + xIndex;
	block[threadIdx.y][threadIdx.x] = idata[index_in];
	
	__syncthreads();
	
	// write the transposed matrix tile to global memory
	xIndex = blockIdx.y * BLOCK_DIM + threadIdx.x;
	yIndex = blockIdx.x * BLOCK_DIM + threadIdx.y;
	unsigned int index_out = yIndex * width + xIndex;
	odata[index_out] = block[threadIdx.x][threadIdx.y];
}

/****************************************************************************
*           Lucas Test - specific routines                                 *
***************************************************************************/

void init_lucas_cu(double *s_inv, double *s_ttp, double *s_ttmp, UL n) {
	int i;
	dim3 grid(512 / BLOCK_DIM, 512 / BLOCK_DIM, 1);
	dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);
	
	cudaMemcpy(g_x,s_inv,sizeof(double)*n,cudaMemcpyHostToDevice);
	for(i=0;i<n;i+=(512*512))
		square_transposef<512><<< grid, threads >>>((float *)&g_inv[i],(double *)&g_x[i]);

	cudaMemcpy(g_x,s_ttp,sizeof(double)*n,cudaMemcpyHostToDevice);
	for(i=0;i<n;i+=(512*512))
		square_transpose<512><<< grid, threads >>>((double *)&g_ttp[i],(double *)&g_x[i]);

	cudaMemcpy(g_x,s_ttmp,sizeof(double)*n,cudaMemcpyHostToDevice);
	for(i=0;i<n;i+=(512*512))
		square_transpose<512><<< grid, threads >>>((double *)&g_ttmp[i],(double *)&g_x[i]);
}

#define IDX(i) ((((i) >> 18) << 18) + (((i) & (512*512-1)) >> 9)  + ((i & 511) << 9))
template <bool g_err_flag, int stride>
__global__ void normalize_kernel(double *g_xx, volatile double *g_maxerr, int *g_carry,
		const float *g_inv, const double *g_ttp, const double *g_ttmp) {
	const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	const int js= stride * threadID;
	
	int carry = (threadID==0) ? -2.0 : 0.0; /* this is the -2 of the LL x*x - 2 */
	
	if (g_err_flag) {
		double err=0.0, tempErr, temp0;
		int idx;
#pragma unroll 32
		for (int j=0; j < stride; j++) {
			idx = IDX(j + js);
			temp0 = g_xx[idx];
			tempErr = RINT( temp0 * g_ttmp[idx] );
			err = fmax(err,fabs( (temp0 * g_ttmp[idx]) - tempErr));
			temp0 = tempErr + carry;
			temp0 *= g_inv[idx];
			carry = RINT(temp0);
			g_xx[idx] = (temp0-carry) * g_ttp[idx];
		}

		g_maxerr[0] = fmax(err, g_maxerr[0]);
	} else {
		double4 buf[4];
		int4 idx;
		double2 temp0;
#pragma unroll 16
		for (int j=0; j < stride; j+=4) {
			// Store first part of IDX calculation for re-use. This saves many registers
			idx.w = (((j + js) >> 18) << 18) + (((j + js) & (512*512-1)) >> 9);
			// Then just add on the bits unique to each idx
			idx.x = idx.w + (((j+js  ) & 511) << 9);
			idx.y = idx.w + (((j+js+1) & 511) << 9);
			idx.z = idx.w + (((j+js+2) & 511) << 9);
			idx.w = idx.w + (((j+js+3) & 511) << 9);

			buf[0].x = g_xx[idx.x];
			buf[1].x = g_ttmp[idx.x];
			buf[2].x = g_inv[idx.x];
			buf[3].x = g_ttp[idx.x];
			
			buf[0].y = g_xx[idx.y];
			buf[1].y = g_ttmp[idx.y];
			buf[2].y = g_inv[idx.y];
			buf[3].y = g_ttp[idx.y];
			
			buf[0].z = g_xx[idx.z];
			buf[1].z = g_ttmp[idx.z];
			buf[2].z = g_inv[idx.z];
			buf[3].z = g_ttp[idx.z];
			
			buf[0].w = g_xx[idx.w];
			buf[1].w = g_ttmp[idx.w];
			buf[2].w = g_inv[idx.w];
			buf[3].w = g_ttp[idx.w];

			temp0.x  = RINT(buf[0].x*buf[1].x);
			temp0.y  = RINT(buf[0].y*buf[1].y);
			temp0.x += carry;
			temp0.x *= buf[2].x;
			carry    = RINT(temp0.x);
			temp0.y += carry;
			temp0.y *= buf[2].y;
			temp0.x  = (temp0.x-carry) * buf[3].x;
			g_xx[idx.x] = temp0.x;
			
			carry    = RINT(temp0.y);
			temp0.y  = (temp0.y-carry) * buf[3].y;
			g_xx[idx.y] = temp0.y;
			
			temp0.x  = RINT(buf[0].z*buf[1].z);
			temp0.y  = RINT(buf[0].w*buf[1].w);
			temp0.x += carry;
			temp0.x *= buf[2].z;
			carry    = RINT(temp0.x);
			temp0.y += carry;
			temp0.y *= buf[2].w;
			temp0.x  = (temp0.x-carry) * buf[3].z;
			g_xx[idx.z] = temp0.x;
			
			carry    = RINT(temp0.y);
			temp0.y  = (temp0.y-carry) * buf[3].w;
			g_xx[idx.w] = temp0.y;
		}
	}

	g_carry[threadID]=carry;
}

__global__ void normalize2_kernel(double *g_xx, const int *g_carry, const UL g_N,
		const float *g_inv, const double *g_ttp, const double *g_ttmp) {
	const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	const int stride = 512;
	const int js= stride * threadID;
	const int je= js + stride;
	UL j;
	double temp0,tempErr;
	int carry;
	int k,ke;
	
	k = (je == g_N) ? 0 : je;
	ke = k + stride;
	
	carry = g_carry[threadID];
	
	for (j=k; carry != 0.0 && j<ke; j+=2) {
		temp0 = g_xx[IDX(j)];
		tempErr = RINT( temp0*g_ttmp[IDX(j)]*(0.5f*g_N));
		temp0 = tempErr + carry;
		temp0 *= g_inv[IDX(j)];
		carry = RINT(temp0);
		g_xx[IDX(j)] = (temp0-carry) * g_ttp[IDX(j)];
		
		temp0 = g_xx[IDX(j+1)];
		tempErr = RINT( temp0*g_ttmp[IDX(j+1)]*-(0.5f*g_N));
		temp0 = tempErr + carry;
		temp0 *= g_inv[IDX(j+1)];
		carry = RINT(temp0);
		g_xx[IDX(j+1)] = (temp0-carry) * g_ttp[IDX(j+1)];
	}
}

void lucas_square_cu(UL N, UL error_log) {
	
	int i;
	
	dim3 grid(512 / BLOCK_DIM, 512 / BLOCK_DIM, 1);
	dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);
	
	for(i=N;i>0;i-=(512*512))
		square_transpose<512><<< grid, threads >>>(&g_x[i],&g_x[i-512*512]);
	
	if (error_log) {
		cutilSafeCall(cudaMemset(g_maxerr, 0, sizeof(double)));
		normalize_kernel<1,512><<< N/512/128,128 >>>(&g_x[512*512],g_maxerr,g_carry,g_inv,g_ttp,g_ttmp);
	} else {
		normalize_kernel<0,512><<< N/512/128,128 >>>(&g_x[512*512],g_maxerr,g_carry,g_inv,g_ttp,g_ttmp);
	}
	
	normalize2_kernel<<< N/512/256,256 >>>(&g_x[512*512],g_carry,N,g_inv,g_ttp,g_ttmp);
	
	for(i=(512*512);i<(N+512*512);i+=(512*512))
		square_transpose<512><<< grid, threads >>>(&g_x[i-512*512],&g_x[i]);
}
