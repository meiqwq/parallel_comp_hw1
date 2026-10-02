const char* sgemm_desc = "Blocked sgemm with 2x2 register tiles.";

#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE 41
#endif

#define min(a,b) (((a)<(b))?(a):(b))

/* This auxiliary subroutine performs a smaller sgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
static void do_block (int lda, int M, int N, int K, float* A, float* B, float* C)
{
  int M2 = M - M % 2;
  int N2 = N - N % 2;

  /* Compute complete 2x2 tiles. */
  for(int j=0;j<N2;j+=2)
    for(int i=0;i<M2;i+=2){
      float cij00 = C[i+j*lda];
      float cij01 = C[i+(j+1)*lda];
      float cij10 = C[(i+1)+j*lda];
      float cij11 = C[(i+1)+(j+1)*lda];

      int k = 0;
      for(;k+1<K;k+=2){
        float a00 = A[i+k*lda];
        float a01 = A[i+(k+1)*lda];
        float a10 = A[(i+1)+k*lda];
        float a11 = A[(i+1)+(k+1)*lda];

        float b00 = B[k+j*lda];
        float b01 = B[k+(j+1)*lda];
        float b10 = B[(k+1)+j*lda];
        float b11 = B[(k+1)+(j+1)*lda];

        cij00 += a00*b00 + a01*b10;
        cij01 += a00*b01 + a01*b11;
        cij10 += a10*b00 + a11*b10;
        cij11 += a10*b01 + a11*b11;
      }

      /* If K is odd, accumulate the last k before storing C. */
      if(k<K){
        float a0 = A[i+k*lda];
        float a1 = A[(i+1)+k*lda];
        float b0 = B[k+j*lda];
        float b1 = B[k+(j+1)*lda];

        cij00 += a0*b0;
        cij01 += a0*b1;
        cij10 += a1*b0;
        cij11 += a1*b1;
      }

      C[i+j*lda] = cij00;
      C[i+(j+1)*lda] = cij01;
      C[(i+1)+j*lda] = cij10;
      C[(i+1)+(j+1)*lda] = cij11;
    }

  /* The remaining row covers only the already paired columns. */
  if(M2<M){
    int i = M2;
    for(int j=0;j<N2;++j){
      float cij = C[i+j*lda];
      for(int k=0;k<K;++k)
        cij += A[i+k*lda] * B[k+j*lda];
      C[i+j*lda] = cij;
    }
  }

  /* The remaining column covers all rows, including the corner once. */
  if(N2<N){
    int j = N2;
    for(int i=0;i<M;++i){
      float cij = C[i+j*lda];
      for(int k=0;k<K;++k)
        cij += A[i+k*lda] * B[k+j*lda];
      C[i+j*lda] = cij;
    }
  }
}

/* General-block kernel for non-square matrices with distinct leading dimensions */
static void do_block_nonsquare (int ldaA, int ldaB, int ldaC,
                                int M, int N, int K,
                                float* A, float* B, float* C)
{
  int M2 = M - M % 2;
  int N2 = N - N % 2;

  /* Each matrix keeps its own leading dimension inside the tile. */
  for(int j=0;j<N2;j+=2)
    for(int i=0;i<M2;i+=2){
      float cij00 = C[i+j*ldaC];
      float cij01 = C[i+(j+1)*ldaC];
      float cij10 = C[(i+1)+j*ldaC];
      float cij11 = C[(i+1)+(j+1)*ldaC];

      int k = 0;
      for(;k+1<K;k+=2){
        float a00 = A[i+k*ldaA];
        float a01 = A[i+(k+1)*ldaA];
        float a10 = A[(i+1)+k*ldaA];
        float a11 = A[(i+1)+(k+1)*ldaA];

        float b00 = B[k+j*ldaB];
        float b01 = B[k+(j+1)*ldaB];
        float b10 = B[(k+1)+j*ldaB];
        float b11 = B[(k+1)+(j+1)*ldaB];

        cij00 += a00*b00 + a01*b10;
        cij01 += a00*b01 + a01*b11;
        cij10 += a10*b00 + a11*b10;
        cij11 += a10*b01 + a11*b11;
      }

      if(k<K){
        float a0 = A[i+k*ldaA];
        float a1 = A[(i+1)+k*ldaA];
        float b0 = B[k+j*ldaB];
        float b1 = B[k+(j+1)*ldaB];

        cij00 += a0*b0;
        cij01 += a0*b1;
        cij10 += a1*b0;
        cij11 += a1*b1;
      }

      C[i+j*ldaC] = cij00;
      C[i+(j+1)*ldaC] = cij01;
      C[(i+1)+j*ldaC] = cij10;
      C[(i+1)+(j+1)*ldaC] = cij11;
    }

  /* Remaining row, excluding the possible last column. */
  if(M2<M){
    int i = M2;
    for(int j=0;j<N2;++j){
      float cij = C[i+j*ldaC];
      for(int k=0;k<K;++k)
        cij += A[i+k*ldaA] * B[k+j*ldaB];
      C[i+j*ldaC] = cij;
    }
  }

  /* Remaining column, including the possible bottom-right corner. */
  if(N2<N){
    int j = N2;
    for(int i=0;i<M;++i){
      float cij = C[i+j*ldaC];
      for(int k=0;k<K;++k)
        cij += A[i+k*ldaA] * B[k+j*ldaB];
      C[i+j*ldaC] = cij;
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

  /* General non-square case with separate leading dimensions. */
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
