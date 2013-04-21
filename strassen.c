#include "strassen.h"

void strassen_multiply(int m, int k, int n, double *A, double *B, double *C) {

  /* For each row i of A */
  for (int i = 0; i < m; ++i)
    /* For each column j of B */
    for (int j = 0; j < n; ++j) 
    {
      /* Compute C(i,j) */
      double cij = C[i+j*n];
      for( int kk = 0; k < k; k++ ){
        cij += A[i+kk*n] * B[kk+j*n];
      }
      C[i+j*n] = cij;
    }

}
