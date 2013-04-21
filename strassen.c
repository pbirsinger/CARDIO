#include "strassen.h"

void strassen_multiply(int m, int k, int n, double *A, double *B, double *C) {
    Problem problem = {m, k, m, k, n, m, A, B, C};
    if should_run_base_case(problem){
      base_case(m, k, n, A, B, C);
    }

    else{
        //assumign m,k,n even for now
        //should i be mallocing for all of these? 
        //might not have to for T0,1,6

        int T_m = m/2;
        int T_k = k/2;
        int S_k = T_k;
        int S_n = n/2;

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

        double *Q0 = (double*) malloc(T_m * S_n * sizeof(double));
        double *Q1 = (double*) malloc(T_m * S_n * sizeof(double));
        double *Q2 = (double*) malloc(T_m * S_n * sizeof(double));
        double *Q3 = (double*) malloc(T_m * S_n * sizeof(double));
        double *Q4 = (double*) malloc(T_m * S_n * sizeof(double));
        double *Q5 = (double*) malloc(T_m * S_n * sizeof(double));
        double *Q6 = (double*) malloc(T_m * S_n * sizeof(double));

        double *U1 = (double*) malloc(T_m * S_n * sizeof(double));
        double *U2 = (double*) malloc(T_m * S_n * sizeof(double));
        double *U3 = (double*) malloc(T_m * S_n * sizeof(double));

        // compute T0 -6
        // compute S0 -6

        // compute Q0-6
        strassen_multiply(T_m, T_k, S_n, T0, S0, Q0);
        strassen_multiply(T_m, T_k, S_n, T1, S1, Q1);
        strassen_multiply(T_m, T_k, S_n, T2, S2, Q2);
        strassen_multiply(T_m, T_k, S_n, T3, S3, Q3);
        strassen_multiply(T_m, T_k, S_n, T4, S4, Q4);
        strassen_multiply(T_m, T_k, S_n, T5, S5, Q5);
        strassen_multiply(T_m, T_k, S_n, T6, S6, Q6);

        //compute U1-3
        

        for (int i=0; i<3; i++){
          //compute U1-3
        }

    }
}

void matrix_add(int m, int k, int n, double *A, double *B, double *C) {

}

void base_case (int m, int k, int n, double *A, double *B, double *C) {
  naive_multiply(m, k, n, A, B, C);
}


void naive_multiply(int m, int k, int n, double *A, double *B, double *C) {
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      double cij = C[i+j*n];
      for( int kk = 0; k < k; k++ ){
        cij += A[i+kk*n] * B[kk+j*n];
      }
      C[i+j*n] = cij;
    }
  }
}
