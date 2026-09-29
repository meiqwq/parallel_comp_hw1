#define _POSIX_C_SOURCE 199309L
#include <stdlib.h> // For: exit, drand48, malloc, free, NULL, EXIT_FAILURE
#include <stdio.h>  // For: perror
#include <string.h> // For: memset

#include <float.h>  // For: DBL_EPSILON
#include <math.h>   // For: fabs

#include <time.h> // For struct timespec, clock_gettime, CLOCK_MONOTONIC
#include <sys/time.h> // For struct timeval, gettimeofday
#include "benchmark-run.h"

/* BLAS reference for general M,N,K */
#define SGEMM sgemm_
extern void SGEMM(char*, char*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
static inline void reference_sgemm_mnk (int M, int N, int K, float alpha, float* A, float* B, float* C)
{
  char TRANSA = 'N';
  char TRANSB = 'N';
  int LDA = M;
  int LDB = K;
  int LDC = M;
  float beta = 1.0f;
  SGEMM(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}

/* Unified SGEMM API */
extern const char* sgemm_desc;
extern void sgemm (int M, int N, int K, float* A, float* B, float* C);

double wall_time ()
{
  struct timespec t;
  clock_gettime (CLOCK_MONOTONIC, &t);
  return 1.*t.tv_sec + 1.e-9*t.tv_nsec;
}


int randint(int l,int u)
{
  int temp;
  srand((unsigned)time(NULL));
  temp = floor(l + (1.0*rand()/RAND_MAX)*(u - l + 1 ));
  return temp;
}


void die (const char* message)
{
  perror (message);
  exit (EXIT_FAILURE);
}

void fill (float* p, int n)
{
  int tt;
  float tmp;
  for (int i = 0; i < n; ++i) {
    tt = rand();
    tmp = (float)tt / (float)(RAND_MAX);
    //printf("%.2lf\n", tmp);
    p[i] = 2 * tmp - 1; // Uniformly distributed over [-1, 1]
  }
}

void absolute_value (float *p, int n)
{
  for (int i = 0; i < n; ++i)
    p[i] = fabs (p[i]);
}

/* The benchmarking program */
int main (int argc, char **argv)
{
  printf ("Description:\t%s\n\n", sgemm_desc);

  /* Square test subset */
  const int square_sizes[] = {31,32,33,63,64,65,95,96,97,127,128,129,159,160,161,191,192,193,223,224,225,255,256,257,287,288,289,319,320,321,384,392,452,472,496,511,512,513,528,575,576,577,640,641,767,768,769,895,896,897,1023,1024,1025};

  /* Non-square quick set */
  const int S = 64, L = 4096;
  const int nonsq[][3] = {{S,S,L},{L,S,S},{S,L,S},{L,S,L},{S,L,L}};

  /* Shared macro now lives in benchmark-run.h */

  /* Run square cases */
  for (unsigned si = 0; si < sizeof(square_sizes)/sizeof(square_sizes[0]); ++si) {
    int n = square_sizes[si];
    RUN_CASE(n,n,n);
  }

  /* Run non-square quick cases */
  for (unsigned ci = 0; ci < sizeof(nonsq)/sizeof(nonsq[0]); ++ci) {
    RUN_CASE(nonsq[ci][0], nonsq[ci][1], nonsq[ci][2]);
  }

#undef RUN_CASE

  return 0;
}
