extern void sgemm_ (char*, char*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*); 

const char* sgemm_desc = "BLAS-backed general sgemm.";

/* This routine performs a sgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values.    
 * This function wraps a call to the BLAS-3 routine sgemm, via the standard FORTRAN interface - hence the reference semantics. */
/* Unified SGEMM API wrapping BLAS for arbitrary M,N,K */
void sgemm (int M, int N, int K, float* A, float* B, float* C)
{
  char TRANSA = 'N';
  char TRANSB = 'N';
  float ALPHA = 1.f;
  float BETA  = 1.f;
  int LDA = M;
  int LDB = K;
  int LDC = M;
  sgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
}   
