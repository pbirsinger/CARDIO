#include "carma.h"
#include "sejits.h"

int dim_to_split(int m, int k, int n) {
  if (n >= k && n >= m) return SPLIT_N;
  if (m >= k && m >= n) return SPLIT_M;
  return SPLIT_K;
}

int get_num_subproblems(Problem p) {
  return 2;
}

Problem get_next_subproblem(Problem p, Problem* subproblems, int prob_num) {
  if (prob_num == 0) {
    int split_dim = dim_to_split(p.m, p.k, p.n);
    if (split_dim == SPLIT_N) {
      Problem subproblem1 = {p.M, p.K, p.m, p.k, p.n/2, p.CM, p.A, p.B, p.C};
      return subproblem1;
    } else if (split_dim == SPLIT_M) {
      Problem subproblem1 = {p.M, p.K, p.m/2, p.k, p.n, p.CM, p.A, p.B, p.C};
      return subproblem1;
    } else { // SPLIT_K
      Problem subproblem1 = {p.M, p.K, p.m, p.k/2, p.n, p.CM, p.A, p.B, p.C};
      return subproblem1;
    }

  } else {
    int split_dim = dim_to_split(p.m, p.k, p.n);
    if (split_dim == SPLIT_N) {
      double *B2 = p.B + p.n/2 * p.K;
      Problem subproblem2 = {p.M, p.K, p.m, p.k, p.n/2, p.CM, p.A, B2, p.C + p.n/2 * p.CM};
      return subproblem2;
    } else if (split_dim == SPLIT_M) {
      double *A2 = p.A + p.m/2;
      Problem subproblem2 = {p.M, p.K, p.m/2, p.k, p.n, p.CM, A2, p.B, p.C + p.m/2};
      return subproblem2;
    } else { // SPLIT_K
      double *A2 = p.A + p.k/2 * p.M;
      double *B2 = p.B + p.k/2;
      double *Q1 = (double*) malloc(p.m * p.n * sizeof(double));
      Problem subproblem2 = {p.M, p.K, p.m, p.k/2, p.n, p.m, A2, B2, Q1};
      return subproblem2;
    }
  }
}

Result merge(Result* results) {
  Result r1 = results[0];
  Result r2 = results[1];

  if (r1.C + r1.n * r1.CM == r2.C) { // split n
    Result r = {r1.m, 2*r1.n, r1.CM, r1.C};
    return r;

  } else if (r1.C + r1.m == r2.C) { // split m
    Result r = {2*r1.m, r1.n, r1.CM, r1.C};
    return r;

  } else { // split k
    int x;
    for (x = 0; x < r1.n; x++) {
      cblas_daxpy(r2.m, 1, r2.C + r2.m * x, 1, r1.C + r1.CM * x, 1);
    }
    free(r2.C);
    Result r = {r1.m, r1.n, r1.CM, r1.C};
    return r;
  }
}

int should_run_base_case(Problem problem, int depth) {
  if (depth >= MAX_DEPTH) {
    return 1;
  } else {
    return 0;
  }
}

Result base_case(Problem p) {
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,p.m,p.n,p.k,1,p.A,p.M,p.B,p.K,0,p.C,p.CM);
  Result r = {p.m, p.n, p.CM, p.C};
  return r;
}

void multiply(int m, int k, int n, double *A, double *B, double *C) {
  Problem problem = {m, k, m, k, n, m, A, B, C};
  solve(problem, 0);
}
