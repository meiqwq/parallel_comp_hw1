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

  /* Test sizes should highlight performance dips at multiples of certain powers-of-two */
  float initial = randint(1,10);
  int test_sizes[] =

  /* Multiples-of-32, +/- 1. for final benchmarking. */
   {31,32,33,63,64,65,95,96,97,127,128,129,159,160,161,191,192,193,223,224,225,255,256,257,287,288,289,319,320,321,351,352,353,383,384,385,415,416,417,447,448,449,479,480,481,511,512,513,543,544,545,575,576,577,607,608,609,639,640,641,671,672,673,703,704,705,735,736,737,767,768,769,799,800,801,831,832,833,863,864,865,895,896,897,927,928,929,959,960,961,991,992,993,1023,1024,1025};

  /* A representative subset of the first list for initial test. Currently uncommented. */
 // { 31, 32, 96, 97, 127, 128, 129, 191, 192, 229, 255, 256, 257,
 //   319, 320, 321, 417, 479, 480, 511, 512, 639, 640, 767, 768, 769 };

  int nsizes = sizeof(test_sizes)/sizeof(test_sizes[0]);

  /* assume last size is also the largest size */
  int nmax = test_sizes[nsizes-1];

  /* allocate memory for all problems */
  float* buf = NULL;
  buf = (float*) malloc (3 * nmax * nmax * sizeof(float));
  if (buf == NULL) die ("failed to allocate largest problem size");

  /* For each test size (square) */
  for (int isize = 0; isize < sizeof(test_sizes)/sizeof(test_sizes[0]); ++isize)
  {
    int n = test_sizes[isize];
    RUN_CASE(n, n, n);
  }

  /* Additional non-square full suite */
  const int s_vals[] = {63,64,65};
  const int l_vals[] = {4095,4096,4097};
  for (int si = 0; si < 3; ++si) {
    for (int li = 0; li < 3; ++li) {
      int S = s_vals[si];
      int L = l_vals[li];
      int cases[][3] = {{ S, S, L }, { L, S, S }, { S, L, S }, { L, S, L }, { S, L, L }};
      for (int ci = 0; ci < 5; ++ci) {
        int M = cases[ci][0], N = cases[ci][1], K = cases[ci][2];
        RUN_CASE(M, N, K);
      }
    }
  }

  free (buf);

  return 0;
}
