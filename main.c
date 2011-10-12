const char *program_name = "CUDALucas"; /* for perror() and similar */
const char program_revision[] = "$Revision: 1.3alpha_ah42$";
char version[sizeof(program_name) + sizeof(program_revision)]; /* overly long, but certain */

/* CUDALucas.c
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

/* Include Files */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include "cuda_safecalls.h"

/* some include needed for mers package */
#include "setup.h"
#include "balance.h"
#include "rw.h"

/* some definitions needed by mers package */
#define kErrLimit (0.35)
#define kErrChkFreq (100)
#define kErrChk (1)

/************************ definitions ********************************
 * This used to try to align to an even multiple of 16 BIG_DOUBLEs in
 *  Guillermo's original version, but adjusting pointers from malloc()
 *  causes core dumps when they're later passed to free().  Is there
 *  another way to do the same thing?  --wedgingt@acm.org */
#ifdef linux
#define ALLOC_DOUBLES(n) ((double *)memalign(128,(n)*sizeof(double)))
#else
#define ALLOC_DOUBLES(n) ((double *)malloc((n)*sizeof(double)))
#endif

/* global variables needed */
double     *two_to_phi, *two_to_minusphi;
extern double     *g_ttp,*g_ttmp;
extern float          *g_inv;
double     high,low,highinv,lowinv;
double     Gsmall,Gbig,Hsmall,Hbig;
UL             b, c;
cufftHandle    plan;
extern double *g_x;
extern double *g_maxerr;
extern double *g_carry;

/********  The TRICK to round to nearest *******************************/
/* This plays with the internal hardward round to nearest when we add a
small number to a very big number, called bigA number. For intel fpu
this is 3*2^62. Then you have to sustract the same number (named different 
to confuse the clever compilers). This makes rint() very fast.

*/

double     bigA,bigB;

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

/* -------- child routines for rdft() -------- */
inline void bitrv2(int n, int *ip, double *a) {
	int j,j1,k,k1,l,m,m2;
	double xr,xi,yr,yi;
	
	ip[0] = 0;
	l = n;
	m = 1;
	
	while ((m << 3) < l) {
		l >>= 1;
		for (j = 0; j < m; j++) {
			ip[m + j] = ip[j] + l;
		}
		m <<= 1;
	}
	m2 = 2 * m;
	if ((m << 3) == l) {
		for (k = 0; k < m; k++) {
			for (j = 0; j < k; j++) {
				j1 = 2 * j + ip[k];
				k1 = 2 * k + ip[j];
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += m2;
				k1 += 2 * m2;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += m2;
				k1 -= m2;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += m2;
				k1 += 2 * m2;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
			}
			
			j1 = 2 * k + m2 + ip[k];
			k1 = j1 + m2;
			xr = a[j1];
			xi = a[j1 + 1];
			yr = a[k1];
			yi = a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
		}
	} else {
		for (k = 1; k < m; k++) {
			for (j = 0; j < k; j++) {
				j1 = 2 * j + ip[k];
				k1 = 2 * k + ip[j];
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += m2;
				k1 += m2;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
			}
		}
	}
}

inline void makewt(int nw, int *ip, double *w) {
	int j, nwh;
	double delta, x, y;
	
	ip[0] = nw;
	ip[1] = 1;
	if (nw > 2) {
		nwh = nw >> 1;
		delta = atan(1.0) / nwh;
		w[0] = 1;
		w[1] = 0;
		w[nwh] = cos(delta * nwh);
		w[nwh + 1] = w[nwh];
		if (nwh > 2) {
			for (j = 2; j < nwh; j += 2) {
				x = cos(delta * j);
				y = sin(delta * j);
				w[j] = x;
				w[j + 1] = y;
				w[nw - j] = y;
				w[nw - j + 1] = x;
			}
			bitrv2(nw, ip + 2, w);
		}
	}
}

inline void makect(int nc, int *ip, double *c) {
	int j,nch;
	double delta;
	
	ip[1] = nc;
	if (nc > 1) {
		nch = nc >> 1;
		delta = atan(1.0) / nch;
		c[0] = cos(delta * nch);
		c[nch] = 0.5 * c[0];
		for (j = 1; j < nch; j++) {
			c[j] = 0.5 * cos(delta * j);
			c[nc - j] = 0.5 * sin(delta * j);
		}
	}
}

extern void rftfsub(int n, double *g_x);

void rdft(int n, double *a, int *ip) {
	int nw, nc;
	
	if(ip[0] == 0){
		nw = n >> 2;
		makewt(nw, ip, &a[n+512*512]);
		nc = ip[1];
		nc = n >> 2;
		makect(nc, ip, &a[n+512*512] + nw);
		cutilSafeCall(cudaMemcpy(g_x, a, sizeof(double)*(n/2*3+512*512), cudaMemcpyHostToDevice));
	}

	cufftSafeCall(cufftExecZ2Z(plan,(cufftDoubleComplex *)g_x,(cufftDoubleComplex *)g_x, CUFFT_INVERSE));
	rftfsub(n,g_x);
	cufftSafeCall(cufftExecZ2Z(plan,(cufftDoubleComplex *)g_x,(cufftDoubleComplex *)g_x, CUFFT_INVERSE));
	return;
}

extern void init_lucas_cu(double *s_inv, double *s_ttp, double *s_ttmp, UL n);

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
	printf("g_inv:    %dmb\n", sizeof(float)*n/1024/1024);
	printf("g_ttp:    %dmb\n", sizeof(double)*n/1024/1024);
	printf("g_ttmp:   %dmb\n", sizeof(double)*n/1024/1024);
*/
	cutilSafeCall(cudaMalloc((void**)&g_x, sizeof(double)*(n/2*3+512*512)));
	cutilSafeCall(cudaMalloc((void**)&g_maxerr, sizeof(double)));
	cutilSafeCall(cudaMalloc((void**)&g_carry, sizeof(double)*n/512));
	cutilSafeCall(cudaMalloc((void**)&g_inv,sizeof(float)*n));
	cutilSafeCall(cudaMalloc((void**)&g_ttp,sizeof(double)*n));
	cutilSafeCall(cudaMalloc((void**)&g_ttmp,sizeof(double)*n));

	cutilSafeCall(cudaMemset(g_carry, 0, sizeof(double)));
	cutilSafeCall(cudaMemset(g_maxerr, 0, sizeof(double)));
	
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

	init_lucas_cu(s_inv, s_ttp, s_ttmp, n);
	
	if (s_inv != NULL) free((char *)s_inv);
	if (s_ttp != NULL) free((char *)s_ttp);
	if (s_ttmp != NULL) free((char *)s_ttmp);
}

inline double last_normalize(double *x,UL N,UL err_flag ) {
	UL i,j,k,bj;
	UL size0;
	double hi=high, hiinv=highinv, lo=low, loinv=lowinv;
	double temp0,tempErr;
	double maxerr=0.0,err=0.0,ttmpSmall=Gsmall,ttmpBig=Gbig,ttmp;
	double carry;
	double A=bigA,B=bigB;
	
	carry = - 2.0; /* this is the -2 of the LL x*x - 2 */
	bj=N;
	size0 = 1;
	
	for (j=0,i=0;j<N;j+=2,i++) {
		ttmp = two_to_minusphi[i];
		temp0 = x[j];
		temp0 *=2.0;
		tempErr = RINT(temp0*ttmp);
		if (err_flag) 
		{
			err = fabs(temp0*ttmp-tempErr);
			if (err>maxerr) maxerr=err;
		}
		temp0 = tempErr + carry;
		if (size0) {
			temp0 *= hiinv;
			carry = RINT(temp0);
			bj+=b;
			ttmp *=ttmpBig;
			if(bj>=N)
				bj -= N; 
			x[j] = (temp0-carry) * hi;
			size0 = (bj>=c);
		} else {
			temp0 *= loinv;
			carry = RINT(temp0);
			bj+=b;
			ttmp *=ttmpSmall;
			if(bj>=N)
				bj -= N; 
			x[j] = (temp0-carry) * lo;
			size0 = (bj>=c);
		}
		temp0 = x[j+1];
		temp0 *=-2.0;
		
		if (j==N-2)
			size0 = 0;
		tempErr = RINT(temp0*ttmp);
		if (err_flag) {
			err = fabs(temp0*ttmp-tempErr);
			if (err>maxerr) maxerr=err;
		}
		temp0 = tempErr + carry;
		if (size0) {
			temp0 *= hiinv;
			carry = RINT(temp0);
			bj+=b;
			ttmp *=ttmpBig;
			if(bj>=N)
				bj -= N; 
			x[j+1] = (temp0-carry) * hi;
			size0 = (bj>=c);
		} else {
			temp0 *= loinv;
			carry = RINT(temp0);
			bj+=b;
			ttmp *=ttmpSmall;
			if(bj>=N)
				bj -= N; 
			x[j+1] = (temp0-carry) * lo;
			size0 = (bj>=c);
		}
	}
	bj = N;
	k=0;
	while(carry != 0) {
		size0 = (bj>=c);
		bj += b;
		temp0 = (x[k] + carry);
		if(bj >= N)
			bj-=N;  
		if (size0) {
			temp0 *= hiinv;
			carry = RINT(temp0);
			x[k] = (temp0 -carry) * hi;
		} else {
			temp0 *= loinv;
			carry = RINT(temp0);
			x[k] = (temp0 -carry)* lo;
		}
		k++;
	}
	return(maxerr);
}

extern void lucas_square_cu(UL N,UL error_log);

double lucas_square(double *x, UL N,UL iter, UL last,UL error_log,int *ip) {
	double err = 0.0;

	rdft(N,x,ip);
	if( iter == last) {
		cutilSafeCall(cudaMemcpy(x,g_x, sizeof(double)*N, cudaMemcpyDeviceToHost));
		err=last_normalize(x,N,error_log);
	} else {
		lucas_square_cu(N, error_log);
		
		if(error_log) {
			double *c_maxerr = (double *)malloc(sizeof(double));

			cutilSafeCall(cudaMemcpy(c_maxerr,g_maxerr, sizeof(double), cudaMemcpyDeviceToHost));
			err=c_maxerr[0];
			
			free (c_maxerr);
		}
	}

	return(err);
}

/* This gives smallest power of two great or equal n */

UL power_of_two_length(UL n) {
	UL i = (((UL)1) << 31), k = 32;
	do {
		k--;
		i >>= 1;
	} while ((i & n) == 0);
	return k;
}

/* Choose the lenght of FFT , n is the power of two preliminary solution,
N is the minimum required length, the return is the estimation of optimal 
lenght.

The estimation is made very rougly. I suposse a prime k pass cost about
k*lengthFFT cpu time (in some units)
*/
UL choose_length(int n) {
	UL bestN=1<<n;
	if (bestN < 524288)
		return(524288);
	return bestN;
}

void init_device() {
	int device_count=0;
	struct cudaDeviceProp properties;
	
	cudaGetDeviceCount( &device_count);
	if (device_number >= device_count) {
		printf("device_number >=  device_count ... exiting\n");
		exit(2);
	}
	
#if CUDART_VERSION >= 4000
	cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync);
	cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
#else
	cudaSetDeviceFlags(cudaDeviceBlockingSync);
#endif
	
	cudaSetDevice(device_number);
	// From Iain
	cudaGetDeviceProperties(&properties, device_number);
	
	if (properties.major == 1 && properties.minor < 3) {
		printf("A GPU with compute capability >= 1.3 is required for double precision arithmetic\n");
		exit(2);
	}
}

/**************************************************************
*
*      Main Function
*
**************************************************************/

int main(int argc, char *argv[]) {
	UL q = 0L, n,  j = 1L, last = 2L, flag;
	size_t k;
	double err, *x =NULL;
	int restarting = 0;
	int *ip = NULL;
	char M = '\0';
	FILE *infp = NULL, *outfp = NULL, *dupfp = NULL;
	
	two_to_phi = NULL; 
	two_to_minusphi = NULL;
	
#if defined(__x86_64__)
	bigA=6755399441055744.0;
#else
	bigA=(((6.0)*0x2000000L)*0x2000000L)*0x800;
#endif
	bigB=bigA;
	
	g_x = NULL;
	
	UL averbits = (sizeof(double)*5)/2 - 1;
	
	strcpy(version, program_name);
	strcat(version, " v");
	strncat(version, program_revision + strlen("4Revision: "), strlen(program_revision) - strlen("4Revison:  4"));
	//strcat(version, " Ballester");
	setup();
	while (!shouldTerminate) {
		do { /* while (restarting) */
			switch ((restarting != 0) ? 3 : input(argc, argv, &q, &n, &j, &err, &x, last, &M, &infp, &outfp, &dupfp, version)) {
				case 0: /* no more input */
					printf("no more input\n");
					return(0);
					
				case 1: /* something wrong */
				default:
					printf("something wrong; error message, if any, already printed\n");
					print_time();
					return(1);
					
				case 2: /* continuing work from a partial result */
					init_device(); //msft
					printf("Resuming from iteration " PRINTF_FMT_UL "\n", j);
					fprintf(outfp, "Testing: M(" PRINTF_FMT_UL ") using FFT size = %dk\n", q, n/2/1024);
					restarting = 1; 
					/* not the usual sense of restarting (FFT errors too high) */
					/* size = n; */ /* supressed */
					break;
					
				case 3:
					init_device(); //msft
					n = (q-1)/averbits +1;
					j = power_of_two_length(n);
					n = choose_length(j);
			
					fprintf(outfp, "Testing: M(" PRINTF_FMT_UL ") using FFT size = %dk\n", q, n/2/1024);
					if (x != NULL)
						cutilSafeCall(cudaFreeHost((char *)x));
					cutilSafeCall(cudaMallocHost((void**) &x,(n+n)*sizeof(double)));
					for (k=1;k<n;k++)
						x[k]=0.0;
					x[0] = 4.0;
					j = 1;
					break;
			}
#ifdef _MSC_VER
			fflush (stdout);
#endif
			if (q <  216091) {
				printf(" too small Exponent\n");
				return 0;
			}

			if (two_to_phi != NULL) free((char *)two_to_phi);
			if (two_to_minusphi != NULL) free((char *)two_to_minusphi);
			
			if (!restarting) {
				if (M != 'U' && M != 'I' && last != q - 2)
					(void)fprintf(stderr, "%s: exponent " PRINTF_FMT_UL " should not need Lucas-Lehmer testing\n", program_name, q);
				if ((M != 'U' && M != 'I') || last == q - 2)
					continue;
			}
			restarting = 0;
			init_lucas(q, n);
			
			cufftSafeCall(cufftPlan1d(&plan, n/2, CUFFT_Z2Z, 1));

			if (ip != NULL)
				free((char *)ip);
			ip = (int *)ALLOC_DOUBLES(((2+(size_t)sqrt((float)n/2))*sizeof(int)));
			ip[0]=0;
			
			if(j==1)
				x[0]*=two_to_minusphi[0]*n;
			last = q - 2; /* the last iteration done in the primary loop */
			int output_frequency = chkpnt_iterations ? chkpnt_iterations : 10000;
			for ( ; !shouldTerminate && !restarting && j <= last; j++) {
#if (kErrChkFreq > 1)
				if ((j % kErrChkFreq) == 1 || j < 1000)
#endif
					flag = kErrChk;
#if (kErrChkFreq > 1)
				else 
					flag = 0;
#endif
				err = lucas_square(x, n, j, last, flag,ip);
				if (chkpnt_iterations != 0 && j % chkpnt_iterations == 2 && j < last) {
					cutilSafeCall(cudaMemcpy(x, g_x, sizeof(double)*n, cudaMemcpyDeviceToHost));
					if(check_point(q, n, (UL) (j+1), err, x) < 0) {
						print_time();
						return(errno);
					}
				}
				
				if (err > kErrLimit) { /* n is not big enough; increase it and start over */
					printf("err = %g, increasing n from %d\n",(double)err,(int)n);
					averbits--;
					restarting = 1;
				}
				
				if ((j % output_frequency) == 0) { 
					if (( j % (20 * output_frequency)) == 0) { /* 25 lines on standard console */
						fprintf(outfp, "Testing: M(" PRINTF_FMT_UL ") using FFT size = %dk\n", q, n/2/1024);
					}

					cutilSafeCall(cudaMemcpy(x,g_x, sizeof(double)*n, cudaMemcpyDeviceToHost));
					printbits(x, q, n, (UL)((q > 64L) ? 64L : q), b, c, high, low, version, outfp, dupfp, output_frequency, j);
				}
			}
			cufftSafeCall(cufftDestroy(plan));
			
			if (restarting) continue;
			else if (j < last) { /* not done, but need to quit, for example, by a SIGTERM signal*/
				chkpnt_iterations = 0L;
				
				cutilSafeCall(cudaMemcpy(x,g_x, sizeof(double)*n, cudaMemcpyDeviceToHost));
				
				return((check_point(q, n, j, err, x) <= 0) ? errno : 0);
			}
		} while (restarting);
		printbits(x, q, n, (UL)((q > 64L) ? 64L : q), b, c, high, low, version, outfp, dupfp, 0, j, true);
	}
	return(0);
}
