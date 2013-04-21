#include "harness.h"
#define MIN_M 1024
#define MAX_M 1024
#define MIN_K 1024
#define MAX_K 1024
#define MIN_N 1024
#define MAX_N 1024

void multiply(int m, int k, int n, double *A, double *B, double *C);

void initialize(int m, int k, int n, double* A, double* B, double* C) {
  srand48(time(NULL));
  int i;
  for(i = 0; i < m*k; i++) A[i] = 2 * drand48() - 1;
  for(i = 0; i < k*n; i++) B[i] = 2 * drand48() - 1;
  for(i = 0; i < m*n; i++) C[i] = 2 * drand48() - 1;
}

int main(int argc, char **argv) {
  srand48(time(NULL));
  int i, m, k, n;
  for (m = MIN_M; m <= MAX_M; m*=2) {
    for (k = MIN_K; k <= MAX_K; k*=2) {
      for (n = MIN_N; n <= MAX_N; n*=2) {
        double *A = (double*) malloc(m * k * sizeof(double));
        double *B = (double*) malloc(k * n * sizeof(double));
        double *C = (double*) malloc(m * n * sizeof(double));
        initialize(m, k, n, A, B, C);
        strassen_multiply(m, k, n, A, B, C);

        // check for correctness
        // memset(C, 0, sizeof(double) * m * n); //if commented, this tests C = A*B instead of C += A*B
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans, m,n,k, -1, A,m, B,k, 1, C,m);
        for(i = 0; i < m*k; i++) A[i] = fabs( A[i] );
        for(i = 0; i < k*n; i++) B[i] = fabs( B[i] );
        for(i = 0; i < m*n; i++) C[i] = fabs( C[i] );
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans, m,n,k, -3.0*DBL_EPSILON*n, A,m, B,k, 1, C,m);
        for(i = 0; i < m*n; i++) {
          if(C[i] > 0) {
            printf("FAILURE: error in matrix multiply exceeds an acceptable margin\n");
            exit(EXIT_FAILURE);
          }
        }

        free(A);
        free(B);
        free(C);
      }
    }
  }
  printf("SUCCESS\n");
  return 0;
}
