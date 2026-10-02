const char* sgemm_desc = "Simple blocked sgemm.";

#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE 41
#endif

#define min(a,b) (((a)<(b))?(a):(b))

/* This auxiliary subroutine performs a smaller sgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
static void do_block (int lda, int M, int N, int K, float* A, float* B, float* C)
{
  for(int j=0;j<N;++j){
    for(int k=0;k<K;++k){
      for(int i=0;i<M;++i){
        C[i+j*lda] += A[i+k*lda] * B[k+j*lda];
      }   
    }
  }
}

/* General-block kernel for non-square matrices with distinct leading dimensions */
static void do_block_nonsquare (int ldaA, int ldaB, int ldaC,
                                int M, int N, int K,
                                float* A, float* B, float* C)
{
  for(int j=0;j<N;++j){
    for(int k=0;k<K;++k){
      for(int i=0;i<M;++i){
        C[i+j*ldaC] += A[i+k*ldaA] * B[k+j*ldaB];
      }   
    }
  }
}

/* This routine performs a sgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format. 
 * On exit, A and B maintain their input values. */  
void square_sgemm (int lda, float* A, float* B, float* C)
{
  /* For each block-row of A */ 
  for (int i = 0; i < lda; i += BLOCK_SIZE)
    /* For each block-column of B */
    for (int j = 0; j < lda; j += BLOCK_SIZE)
      /* Accumulate block sgemms into block of C */
      for (int k = 0; k < lda; k += BLOCK_SIZE)
      {
	/* Correct block dimensions if block "goes off edge of" the matrix */
	int M = min (BLOCK_SIZE, lda-i);
	int N = min (BLOCK_SIZE, lda-j);
	int K = min (BLOCK_SIZE, lda-k);

	/* Perform individual block sgemm */
	do_block(lda, M, N, K, A + i + k*lda, B + k + j*lda, C + i + j*lda);
      }
}

void sgemm(int M, int N, int K, float* A, float* B, float* C)
{
  if (M == N && K == N) {
    /* Use the optimized blocked routine for square shapes */
    square_sgemm(M, A, B, C);
    return;
  }

  /* General non-square case: use a straightforward column-major triple loop */
  for (int i = 0; i < M; i += BLOCK_SIZE)
    for (int j = 0; j < N; j += BLOCK_SIZE)
      for (int k = 0; k < K; k += BLOCK_SIZE)
      {
        int blockM = min (BLOCK_SIZE, M - i);
        int blockN = min (BLOCK_SIZE, N - j);
        int blockK = min (BLOCK_SIZE, K - k);
        do_block_nonsquare(M, K, M, blockM, blockN, blockK,
                           A + i + k*M,
                           B + k + j*K,
                           C + i + j*M);
      }
}