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
#include <cufft.h>
#include "cuda_safecalls.h"
#include "setup.h"

/************************ definitions ************************************/


/* This used to try to align to an even multiple of 16 BIG_DOUBLEs in */
/*  Guillermo's original version, but adjusting pointers from malloc() */
/*  causes core dumps when they're later passed to free().  Is there */
/*  another way to do the same thing?  --wedgingt@acm.org */
#ifdef linux
#define ALLOC_DOUBLES(n) ((double *)memalign(128,(n)*sizeof(double)))
#else
#define ALLOC_DOUBLES(n) ((double *)malloc((n)*sizeof(double)))
#endif
/* global variables needed */
double     *two_to_phi, *two_to_minusphi;
double     *g_ttp,*g_ttmp;
float          *g_inv;
double     high,low,highinv,lowinv;
double     Gsmall,Gbig,Hsmall,Hbig;
UL             b, c;
cufftHandle    plan;
double *g_x;
double *g_maxerr;
double *g_carry;

/* rint is not ANSI compatible, so we need a definition for 
* WIN32 and other platforms with rint.
* Also we use that to write the trick to rint()
*/

#if defined(__x86_32__)
#define RINT(x) (floor(x+0.5))	
#else
#define RINT(x) (((x) + A ) - B)
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

__global__ void rftfsub_kernel(int n, double *a) {
	const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	double wkr, wki, xr, xi, yr, yi,cc,d,aj,aj1,ak,ak1, *c ;
	
	if(threadID != 0) {
		int j = threadID * 2;
		c = &a[n+n/4+512*512];
		int nc = n >> 2 ;
		
		aj  = a[j];
		ak  = a[(n-j)];
		
		wkr = c[nc-j/2];
		wki = c[j/2];
		
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
		int m = n >> 1;
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

void rftfsub(int n, double *g_x) {
	dim3 grid(n/1024,1,1);
	dim3 threads(256,1,1);
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
__global__ void square_transpose(double* __restrict__ odata, double* __restrict__ idata) {
	__shared__ double block[BLOCK_DIM][BLOCK_DIM+1];

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

template <int width>
__global__ void square_transposef(float* __restrict__ odata, double* __restrict__ idata) {
	__shared__ double block[BLOCK_DIM][BLOCK_DIM+1];
	
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

void init_lucas(UL q, UL n) {
	UL j,qn,a,i,done;
	UL size0,bj;
	double log2 = log(2.0);
	double ttp,ttmp;
	double *s_inv,*s_ttp,*s_ttmp;
	
	two_to_phi = ALLOC_DOUBLES(n/2);
	two_to_minusphi = ALLOC_DOUBLES(n/2);
	s_inv = ALLOC_DOUBLES(n);
	s_ttp = ALLOC_DOUBLES(n);
	s_ttmp = ALLOC_DOUBLES(n);
	
	if (g_x != NULL) {
		cutilSafeCall(cudaFree((char *)g_x));
		cutilSafeCall(cudaFree((char *)g_maxerr));
		cutilSafeCall(cudaFree((char *)g_carry));
		cutilSafeCall(cudaFree((char *)g_inv));
		cutilSafeCall(cudaFree((char *)g_ttp));
		cutilSafeCall(cudaFree((char *)g_ttmp));
	}
	
/*	printf("Attempting to malloc:\n");
	printf("g_x:      %dmb\n", sizeof(double)*(n/2*3+512*512)/1024/1024);
	printf("g_maxerr: %db \n", sizeof(double));
	printf("g_carry:  %dkb\n", sizeof(double)*n/512/1024);
	printf("g_inv:    %dmb\n", sizeof(double)*n/2/1024/1024);
	printf("g_ttp:    %dmb\n", sizeof(double)*n/1024/1024);
	printf("g_ttmp:   %dmb\n", sizeof(double)*n/1024/1024);
*/
	cutilSafeCall(cudaMalloc((void**)&g_x, sizeof(double)*(n/2*3+512*512)));
	cutilSafeCall(cudaMalloc((void**)&g_maxerr, sizeof(double)));
	cutilSafeCall(cudaMalloc((void**)&g_carry, sizeof(double)*n/512));
	cutilSafeCall(cudaMalloc((void**)&g_inv,sizeof(double)*n/2));
	cutilSafeCall(cudaMalloc((void**)&g_ttp,sizeof(double)*n));
	cutilSafeCall(cudaMalloc((void**)&g_ttmp,sizeof(double)*n));
	
	cufftSafeCall(cufftPlan1d(&plan, n/2, CUFFT_Z2Z, 1));
	
	low = floor((exp(floor((double)q/n)*log2))+0.5);
	high = low+low;
	lowinv = 1.0/low;
	highinv = 1.0/high;
	b = q % n;
	c = n-b;
	
	two_to_phi[0] = 1.0;
	two_to_minusphi[0] = 1.0/(double)(n);
	qn = (b*2)%n;
	
	for(i=1,j=2; j<n; j+=2,i++) {
		a = n - qn;
		two_to_phi[i] = exp(a*log2/n);
		two_to_minusphi[i] = 1.0/(two_to_phi[i]*n);
		qn+=b*2;
		qn%=n;
	}
	
	Hbig = exp(c*log2/n);
	Gbig = 1/Hbig;
	done = 0;
	j = 0;
	while (!done) {
		if (!((j*b) % n >=c || j==0)) {
			a = n -((j+1)*b)%n;
			i = n -(j*b)%n;
			Hsmall = exp(a*log2/n)/exp(i*log2/n);
			Gsmall = 1/Hsmall;
			done = 1;
		}
		j++;
	}
	bj = n;
	size0 = 1;
	bj = n - 1 * b;
	
	for (j=0,i=0; j<n; j=j+2,i++) {
		ttmp = two_to_minusphi[i];
		ttp = two_to_phi[i];
		
		bj += b;
		bj = bj & (n-1);
		size0 = (bj>=c);

		if(j == 0) size0 = 1;

		s_ttmp[j]=ttmp*2.0;

		if (size0) {
			s_inv[j]=highinv;
			ttmp *=Gbig;
			s_ttp[j]=ttp * high;
			ttp *= Hbig;
		} else {
			s_inv[j]=lowinv;
			ttmp *=Gsmall;
			s_ttp[j]=ttp * low;
			ttp *= Hsmall;
		}

		bj += b;
		bj = bj & (n-1);
		size0 = (bj>=c);
		
		if (j==(n-2)) size0 = 0;

		s_ttmp[j+1]=ttmp*-2.0;

		if (size0) {
			s_inv[j+1]=highinv;
			s_ttp[j+1]=ttp * high;
		} else {
			s_inv[j+1]=lowinv;
			s_ttp[j+1]=ttp * low;
		}
	}

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

	if (s_inv != NULL) free((char *)s_inv);
	if (s_ttp != NULL) free((char *)s_ttp);
	if (s_ttmp != NULL) free((char *)s_ttmp);
}

#define IDX(i) ((((i) >> 18) << 18) + (((i) & (512*512-1)) >> 9)  + ((i & 511) << 9))
template <bool g_err_flag, int stride>
__global__ void normalize_kernel(double *g_xx, double A, double B, volatile double *g_maxerr,
		double *g_carry, float *g_inv, double *g_ttp, double *g_ttmp) {
	const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	const int js= stride * threadID;
	
	double carry = (threadID==0) ? -2.0 : 0.0; /* this is the -2 of the LL x*x - 2 */
	
	if (g_err_flag) {
		double err=0.0, tempErr, temp0;
#pragma unroll 4
		for (int j=0; j < stride; j++) {
			temp0 = g_xx[IDX(j + js)];
			tempErr = RINT( temp0 * g_ttmp[IDX(j + js)] );
			err = max(err,fabs( (temp0 * g_ttmp[IDX(j + js)]) - tempErr));
			temp0 = tempErr + carry;
			temp0 *= g_inv[IDX(j + js)];
			carry = RINT(temp0);
			g_xx[IDX(j + js)] = (temp0-carry) * g_ttp[IDX(j + js)];
		}
		
		if (err > g_maxerr[0]) {
			g_maxerr[0] = err;
		}
	} else {
		double4 buf[4];
		int4 idx;
		double2 temp0;
#pragma unroll 1
		for (int j=0; j < stride; j+=4) {
			idx.x = IDX((j + js));
			idx.y = IDX((j + js) + 1);
			idx.z = IDX((j + js) + 2);
			idx.w = IDX((j + js) + 3);
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
			
			carry    = RINT(temp0.y);
			temp0.y  = (temp0.y-carry) * buf[3].y;
			
			g_xx[idx.x] = temp0.x;
			g_xx[idx.y] = temp0.y;
			
			temp0.x  = RINT(buf[0].z*buf[1].z);
			temp0.y  = RINT(buf[0].w*buf[1].w);
			temp0.x += carry;
			temp0.x *= buf[2].z;
			carry    = RINT(temp0.x);
			temp0.y += carry;
			temp0.y *= buf[2].w;
			temp0.x  = (temp0.x-carry) * buf[3].z;
			
			carry    = RINT(temp0.y);
			temp0.y  = (temp0.y-carry) * buf[3].w;
			
			g_xx[idx.z] = temp0.x;
			g_xx[idx.w] = temp0.y;
		}
	}
	
	g_carry[threadID]=carry;
}

__global__ void normalize2_kernel(double *g_xx, double A, double B,
		double *g_carry, UL g_N, float *g_inv, double *g_ttp, double *g_ttmp) {
	const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
	const int stride = 512;
	const int js= stride * threadID;
	const int je= js + stride;
	UL j;
	double temp0,tempErr;
	double carry;
	int k,ke;
	
	k = (je == g_N) ? 0 : je;
	ke = k + stride;
	
	carry = g_carry[threadID];
	
	for (j=k; carry != 0.0 && j<ke; j+=2) {
		temp0 = g_xx[IDX(j)];
		tempErr = RINT( temp0*g_ttmp[IDX(j)]*0.5*g_N );
		temp0 = tempErr + carry;
		temp0 *= g_inv[IDX(j)];
		carry = RINT(temp0);
		g_xx[IDX(j)] = (temp0-carry) * g_ttp[IDX(j)];
		
		temp0 = g_xx[IDX(j+1)];
		tempErr = RINT( temp0*g_ttmp[IDX(j+1)]*(-0.5)*g_N );
		temp0 = tempErr + carry;
		temp0 *= g_inv[IDX(j+1)];
		carry = RINT(temp0);
		g_xx[IDX(j+1)] = (temp0-carry) * g_ttp[IDX(j+1)];
	}
}

void lucas_square_cu(double *x, UL N,UL iter, UL last,UL error_log,int *ip) {
	
	unsigned int i;
	double bigAB=6755399441055744.0;
	
	dim3 grid(512 / BLOCK_DIM, 512 / BLOCK_DIM, 1);
	dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);
	
	for(i=N;i>0;i-=(512*512))
		square_transpose<512><<< grid, threads >>>(&g_x[i],&g_x[i-512*512]);
	
	if (error_log) {
		normalize_kernel<1,512><<< N/512/128,128 >>>(&g_x[512*512], bigAB,bigAB,g_maxerr,g_carry,g_inv,g_ttp,g_ttmp);
	} else {
		normalize_kernel<0,512><<< N/512/128,128 >>>(&g_x[512*512], bigAB,bigAB,g_maxerr,g_carry,g_inv,g_ttp,g_ttmp);
	}
	
	normalize2_kernel<<< N/512/256,256 >>>(&g_x[512*512], bigAB,bigAB,g_carry,N,g_inv,g_ttp,g_ttmp);
	
	for(i=(512*512);i<(N+512*512);i+=(512*512))
		square_transpose<512><<< grid, threads >>>(&g_x[i-512*512],&g_x[i]);
}
