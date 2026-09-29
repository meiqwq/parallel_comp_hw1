const char* sgemm_desc = "Naive, three-loop sgemm.";

/* This routine performs a sgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */    
void square_sgemm (int n, float* A, float* B, float* C)
{
  /* For each row i of A */
  for (int i = 0; i < n; ++i)
    /* For each column j of B */
    for (int j = 0; j < n; ++j) 
    {
      /* Compute C(i,j) */
      float cij = C[i+j*n];
      for( int k = 0; k < n; k++ )
	       cij += A[i+k*n] * B[k+j*n];
      C[i+j*n] = cij;
    }
}

/* Unified general SGEMM API (column-major):
 * C(M,N) += A(M,K) * B(K,N)
 */
void sgemm (int M, int N, int K, float* A, float* B, float* C)
{
  for (int j = 0; j < N; ++j)
    for (int i = 0; i < M; ++i)
    {
      float cij = C[i + j*M];
      for (int k = 0; k < K; ++k)
        cij += A[i + k*M] * B[k + j*K];
      C[i + j*M] = cij;
    }
}
