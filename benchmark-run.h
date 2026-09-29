#ifndef BENCHMARK_RUN_H
#define BENCHMARK_RUN_H

/* Shared RUN_CASE macro used by benchmark.c and benchmark-test.c
 * Expects the following symbols to be available in the including TU:
 *   - wall_time(), fill(), absolute_value(), reference_sgemm_mnk(), die()
 *   - sgemm(M,N,K,A,B,C)
 */
#define RUN_CASE(M,N,K) \
  do { \
    float* A = (float*) malloc ((M)*(K) * sizeof(float)); \
    float* B = (float*) malloc ((K)*(N) * sizeof(float)); \
    float* C = (float*) malloc ((M)*(N) * sizeof(float)); \
    if (!A || !B || !C) die("alloc failed"); \
    fill(A,(M)*(K)); fill(B,(K)*(N)); fill(C,(M)*(N)); \
    double Gflops_s, seconds=-1.0; double timeout=0.1; int niter=0; \
    for (niter=1; seconds<timeout; ) { \
      niter*=2; sgemm((M),(N),(K),A,B,C); \
      seconds=-wall_time(); \
      for (int it=0; it<niter; ++it) sgemm((M),(N),(K),A,B,C); \
      seconds+=wall_time(); \
      Gflops_s = 2.e-9 * (double)niter * (double)(M) * (double)(N) * (double)(K) / seconds; \
    } \
    printf("M=%d N=%d K=%d\tGflop/s: %.3g (%d iter, %.3f s)\n", (M),(N),(K), Gflops_s, niter, seconds); \
    memset(C,0,(M)*(N)*sizeof(float)); \
    sgemm((M),(N),(K),A,B,C); \
    reference_sgemm_mnk((M),(N),(K), -1.0f, A,B,C); \
    absolute_value(A,(M)*(K)); absolute_value(B,(K)*(N)); absolute_value(C,(M)*(N)); \
    reference_sgemm_mnk((M),(N),(K), -3.0f*FLT_EPSILON*(float)(K), A,B,C); \
    for (int ii=0; ii<(M)*(N); ++ii) if (C[ii] > 0.0f) die("*** FAILURE *** Error exceeds bounds.\n"); \
    free(A); free(B); free(C); \
  } while(0)

#endif


