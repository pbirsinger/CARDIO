#include "strassen.h"

// A is m * k
// B is k * n
// C is m * n

void strassen_multiply(int m, int k, int n, double *A, double *B, double *C) {
    Problem problem = {m, k, m, k, n, m, A, B, C};
    if should_run_base_case(problem){
        base_case(m, k, n, A, B, C);
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

        //set T0, T1, T6
        for (int i=0; i < T_k; i++){
            for (int j=0; j < T_m; j++){
                T0[j + i*T_m] = A[j + i*m];
                T1[j + i*T_m] = A[j + (i+T_k)*m];

                T2[j + i*T_m] = A[j + T_m+ i*m]; // still need to add

                T4[j + i*T_m] = A[j + i*m];  //still need to subtract
                T5[j + i*T_m] = A[j + (i+T_k)*m]; //still need to add

                T6[j + i*T_m] = A[j +T_m + (i+T_k)*m];
            }
        }

        //set S0, S1, S5
        for (int i=0; i < S_n; i++){
            for (int j=0; j < S_k; j++){
                S0[j + i*S_k] = B[j + i*k];
                S1[j + i*S_k] = B[j +S_k + i*k];

                S2[j + i*S_k] = B[j +(i*S_n)*k]; //**
                S3[j + i*S_k] = B[j +S_k + (i+S_n)*m]; //**
                S4[j + i*S_k] = B[j +S_k + (i+S_n)*m]; //**

                S5[j + i*S_k] = B[j +S_k + (i+S_n)*m];
            }
        }

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

        for (int i=0; i <S_n; i++){
            for (int j=0; j <T_m; j++){
                C[j + i*m] = C11[j+i*T_m];
                C[j + (i+S_n)*m] = C12[j+i*T_m];
                C[j + T_m + i*m] = C21[j+i*T_m];
                C[j + T_m + (i+S_n)*m] = C22[j+i*T_m];
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
        free(C0);
        free(C1);
        free(C2);
        free(C3);
    }
}

int should_run_base_case(Problem p) {
  if (p.m == 1 || p.k == 1 || p.n == 1){
    return 1;
  }
  else{
    return 0;
  }
}

void matrix_add(int m, int n, double *A, double *B, double *C) {
    for (int i =0; i<m; i++){
        for (int j =0; j <n; j++){
            C[i + j*m] = A[i+j*m] + B[i+j*m];
        }
    }
}

void matrix_subtract(int m, int n, double *A, double *B, double *C) {
    for (int i =0; i<m; i++){
        for (int j =0; j <n; j++){
            C[i + j*m] = A[i+j*m] - B[i+j*m];
        }
    }
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
