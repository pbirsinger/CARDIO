#include "strassen.h"

// A is m * k
// B is k * n
// C is m * n

/**
FUTURE OPTIMIZATIONS:
for additions/ subtractions
use daxpy
    look at example in carma.c

for base case
    use mkl
        example in carma
    dgemm
        google and look at wiki
**/

int should_run_base_case(Problem p) {
  if (p.m == 32 || p.k == 32 || p.n == 32){
    return 1;
  }
  else{
    return 0;
  }
}

void matrix_add(int m, int n, double *A, double *B, double *C) {
    int i,j;
    for (i =0; i<m; i++){
        for (j =0; j <n; j++){
            C[i + j*m] = A[i+j*m] + B[i+j*m];
        }
    }
}

void matrix_subtract(int m, int n, double *A, double *B, double *C) {
    int i,j;
    for (i =0; i<m; i++){
        for (j =0; j <n; j++){
            C[i + j*m] = A[i+j*m] - B[i+j*m];
        }
    }
}

void naive_multiply(int m, int k, int n, double *A, double *B, double *C) {
    int i,j,kk;
    double cij;
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
          cij = 0.0;
          for( kk = 0; kk < k; kk++ ){
            cij += A[i+kk*m] * B[kk+j*k];
          }
          C[i+j*n] = cij;
        }
    }
}

void base_case (int m, int k, int n, double *A, double *B, double *C) {
  naive_multiply(m, k, n, A, B, C);
}

void strassen_multiply(int m, int k, int n, double *A, double *B, double *C) {
    Problem problem = {m, k, m, k, n, m, A, B, C};
    int i =0;
    int j =0;
    if (should_run_base_case(problem)){
        base_case(m, k, n, A, B, C);
        //switch to MKL 
    }
    else{
        //assumign m,k,n even for now
        //should i be mallocing for all of these? 
        //might not have to for T0,1,6
        //assuming stored column major

        int T_m = m/2;
        int T_k = k/2;
        int S_k = T_k;
        int S_n = n/2;

        //compute T0-6, S0-6
        double *T0 = (double*) malloc(T_m * T_k * sizeof(double));
        double *T1 = (double*) malloc(T_m * T_k * sizeof(double));
        double *T2 = (double*) malloc(T_m * T_k * sizeof(double));
        double *T3 = (double*) malloc(T_m * T_k * sizeof(double));
        double *T4 = (double*) malloc(T_m * T_k * sizeof(double));
        double *T5 = (double*) malloc(T_m * T_k * sizeof(double));
        double *T6 = (double*) malloc(T_m * T_k * sizeof(double));

        double *S0 = (double*) malloc(S_k * S_n * sizeof(double));
        double *S1 = (double*) malloc(S_k * S_n * sizeof(double));
        double *S2 = (double*) malloc(S_k * S_n * sizeof(double));
        double *S3 = (double*) malloc(S_k * S_n * sizeof(double));
        double *S4 = (double*) malloc(S_k * S_n * sizeof(double));
        double *S5 = (double*) malloc(S_k * S_n * sizeof(double));
        double *S6 = (double*) malloc(S_k * S_n * sizeof(double));


        for (i = 0; i < T_k; i++){
            memcpy(T0 +i * T_m, A + i*m, T_m * sizeof(double));
            memcpy(T1 +i * T_m, A +(i+T_k) *m, T_m * sizeof(double));

            memcpy(T2 +i * T_m, A +T_m+ i*m, T_m * sizeof(double));
            memcpy(T4 +i * T_m, A + i*m, T_m * sizeof(double));

            memcpy(T5 +i * T_m, A + (i+T_k)*m, T_m * sizeof(double));
            memcpy(T6 +i * T_m, A +T_m+ (i+T_k)*m, T_m * sizeof(double));
        }

        for (i = 0; i < S_n; i++){
            memcpy(S0 +i * S_k, B + i*k, S_k * sizeof(double));

            memcpy(S1 +i * S_k, B +S_k+ i*k, S_k * sizeof(double));
            memcpy(S2 +i * S_k, B + (i+S_n)*k, S_k * sizeof(double));
            memcpy(S3 +i * S_k, B + S_k+(i+S_n)*k, S_k * sizeof(double));
            memcpy(S4 +i * S_k, B +S_k+ (i+S_n)*k, S_k * sizeof(double));
            memcpy(S5 +i * S_k, B+ S_k+(i+S_n)*k, S_k * sizeof(double));
        }

        //use daxpy for adds and subtracts 
        //cblas_daxpy()

        matrix_subtract(T_m, T_k, T0, T2, T4);
        matrix_add(T_m, T_k, T2, T6, T2);
        matrix_subtract(T_m, T_k, T2, T1, T3);
        matrix_add(T_m, T_k, T1, T3, T5);

        matrix_subtract(S_k, S_n, S4, S2, S4);
        matrix_add(S_k, S_n, S2, S0, S2);
        matrix_subtract(S_k, S_n, S3, S2, S3);
        matrix_subtract(S_k, S_n, S3, S1, S6);

        // compute Q0-6
        double *Q0 = (double*) malloc(T_m * S_n * sizeof(double));
        double *Q1 = (double*) malloc(T_m * S_n * sizeof(double));
        double *Q2 = (double*) malloc(T_m * S_n * sizeof(double));
        double *Q3 = (double*) malloc(T_m * S_n * sizeof(double));
        double *Q4 = (double*) malloc(T_m * S_n * sizeof(double));
        double *Q5 = (double*) malloc(T_m * S_n * sizeof(double));
        double *Q6 = (double*) malloc(T_m * S_n * sizeof(double));

        strassen_multiply(T_m, T_k, S_n, T0, S0, Q0);
        strassen_multiply(T_m, T_k, S_n, T1, S1, Q1);
        strassen_multiply(T_m, T_k, S_n, T2, S2, Q2);
        strassen_multiply(T_m, T_k, S_n, T3, S3, Q3);
        strassen_multiply(T_m, T_k, S_n, T4, S4, Q4);
        strassen_multiply(T_m, T_k, S_n, T5, S5, Q5);
        strassen_multiply(T_m, T_k, S_n, T6, S6, Q6);

        //compute U1-3
        double *U1 = (double*) malloc(T_m * S_n * sizeof(double));
        double *U2 = (double*) malloc(T_m * S_n * sizeof(double));
        double *U3 = (double*) malloc(T_m * S_n * sizeof(double));

        matrix_add(T_m, S_n, Q0, Q3, U1);
        matrix_add(T_m, S_n, U1, Q4, U2);
        matrix_add(T_m, S_n, U1, Q2, U3);

        //compute C
        double *C11 = (double*) malloc(T_m * S_n * sizeof(double));
        double *C12 = (double*) malloc(T_m * S_n * sizeof(double));
        double *C21 = (double*) malloc(T_m * S_n * sizeof(double));
        double *C22 = (double*) malloc(T_m * S_n * sizeof(double));

        matrix_add(T_m, S_n, Q0, Q1, C11);
        matrix_add(T_m, S_n, U3, Q5, C12);
        matrix_subtract(T_m, S_n, U2, Q6, C21);
        matrix_add(T_m, S_n, U2, Q2, C22);

        for (i=0; i <S_n; i++){
            for (j=0; j <T_m; j++){
                C[j + i*m] = C11[j+i*T_m];
                C[j + (i+S_n)*m] = C12[j+i*T_m];
                C[j + T_m + i*m] = C21[j+i*T_m];
                C[j + T_m + (i+S_n)*m] = C22[j+i*T_m];
            }
        }

        free(T0);
        free(T1);
        free(T2);
        free(T3);
        free(T4);
        free(T5);
        free(T6);
        free(S0);
        free(S1);
        free(S2);
        free(S3);
        free(S4);
        free(S5);
        free(S6);
        free(Q0);
        free(Q1);
        free(Q2);
        free(Q3);
        free(Q4);
        free(Q5);
        free(Q6);
        free(U1);
        free(U2);
        free(U3);
        free(C11);
        free(C12);
        free(C21);
        free(C22);
    }
}




