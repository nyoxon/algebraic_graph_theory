#include "eig.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>

static void rand_unit_vector(double* x, int n) {
  double s = 0.0;

  for (int i = 0; i < n; i++) {
    x[i] = ((double) rand() / RAND_MAX) - 0.5;
    s += x[i] * x[i];
  }

  s = sqrt(s);

  if (s == 0.0) {
    x[0] = 1.0;
    s = 1.0;
  }

  for (int i = 0; i < n; i++) {
    x[i] /= s;
  }
}

static double dot(const double* a, const double* b, int n) {
  double s = 0.0;

  for (int i = 0; i < n; i++) {
    s += a[i] * b[i];
  }

  return s;
}

static double norm2(const double* a, int n) {
  return sqrt(dot(a, a, n));
}

static void matvec_int
(
  const int** M,
  const double* x,
  double* y,
  int n
)
{
  for (int i = 0; i < n; i++) {
    double acc = 0.0;

    const int* row = (const int*)M[i];

    for (int j = 0; j < n; j++) {
      acc += (double) row[j] * x[j];
    }

    y[i] = acc;
  }
}

static double rayleigh_int
(
  const int** M,
  const double* x,
  int n,
  double* tmp
)
{
  matvec_int(M, x, tmp, n);
  double num = dot(x, tmp, n);
  double den = dot(x, x, n);

  return num / den;
}

double power_method_int
(
    const int** M,
  int n,
  int maxit,
  double tol,
  double* out_vec  
)
{
  double* x = malloc(n * sizeof(double));
  double* y = malloc(n * sizeof(double));

  rand_unit_vector(x, n);

  double lambda = 0.0, lambda_old = 0.0;

  for (int it = 0; it < maxit; it++) {
    matvec_int(M, x, y, n);

    double ny = norm2(y, n);

    if (ny == 0.0) {
      rand_unit_vector(x, n);
      continue;
    }

    for (int i = 0; i < n; i++) {
      x[i] = y[i] / ny;
    }

    lambda_old = lambda;
    lambda = rayleigh_int(M, x, n, y);

    double denom = fmax(1.0, fabs(lambda));

    if (fabs(lambda - lambda_old) / denom < tol) {
      break;
    }
  }

  for (int i = 0; i < n; i++) {
    out_vec[i] = x[i];
  }

  free(x);
  free(y);

  return lambda;
}

// Calcula norma do subvetor (start, n) de x;
static double vnorm2(const double* x, int start, int n) {
  double s = 0.0;

  for (int i = start; i < n; i++) {
    s += x[i] * x[i];
  }

  return sqrt(s);
}


// Aplica a reflexão de householder à esquerda de A

static void householder_left(double **A, int n, int k, double *u){
    // constrói u a partir da coluna k
    double sigma = vnorm2(&A[k][k], 0, n-k);
    if (sigma==0.0) { for (int i=0;i<n;i++) u[i]=0.0; return; }

    // vetor v = coluna k (i>=k)
    for (int i=0;i<k;i++) u[i]=0.0;
    for (int i=k;i<n;i++) u[i]=A[i][k];

    double sign = (u[k]>=0.0)? 1.0 : -1.0;
    u[k] += sign*sigma;
    // normaliza u
    double nu = 0.0; for (int i=k;i<n;i++) nu += u[i]*u[i];
    nu = sqrt(nu); if (nu==0.0){ for (int i=0;i<n;i++) u[i]=0.0; return; }
    for (int i=k;i<n;i++) u[i] /= nu;

    // A <- (I - 2 u u^T) A
    // compute w = 2 * u^T A (linhas k..n-1)
    for (int j=k;j<n;j++){
        double dot=0.0;
        for (int i=k;i<n;i++) dot += u[i]*A[i][j];
        double factor = 2.0*dot;
        for (int i=k;i<n;i++) A[i][j] -= factor*u[i];
    }
}
// Aplica a reflexão de householder à direita

static void householder_right(double **A, int n, int k, const double *u){
    // A <- A - 2*(A u) u^T
    double *Au = malloc(n*sizeof(double));
    for (int i=0;i<n;i++){
        double s=0.0;
        for (int j=k;j<n;j++) s += A[i][j]*u[j];
        Au[i]=s;
    }
    for (int i=0;i<n;i++){
        double factor = 2.0*Au[i];
        for (int j=k;j<n;j++) A[i][j] -= factor*u[j];
    }
    free(Au);
}


static void qr_householder(double **A, int n, double **Q_out){
    if (Q_out){
        for (int i=0;i<n;i++)
            for (int j=0;j<n;j++)
                Q_out[i][j] = (i==j)?1.0:0.0;
    }
    double *u = malloc(n*sizeof(double));
    for (int k=0;k<n-1;k++){
        householder_left(A, n, k, u);
        if (Q_out) householder_right(Q_out, n, k, u); // k, não 0
    }
    free(u);
}
/* Transforma A numa matriz onde os elementos da diagonal
são seus autovalores.
Este algoritmo tem complexidade O(n³) e
é razoavel para matrizes pequenas. */
void qr_iterate(double **A, int n, int maxit, double tol){
    double **Q = malloc(n*sizeof(*Q));
    double **R = malloc(n*sizeof(*R));
    for (int i=0;i<n;i++){ Q[i]=malloc(n*sizeof(**Q)); R[i]=malloc(n*sizeof(**R)); }

    for (int it=0; it<maxit; it++){
        // R <- A (vamos fatorar R = Q*R via Householder à esquerda)
        for (int i=0;i<n;i++) for (int j=0;j<n;j++) R[i][j]=A[i][j];

        // QR por Householder: R vira triangular sup; Q_out recebe Q
        qr_householder(R, n, Q);

        // A <- R * Q
        for (int i=0;i<n;i++){
            for (int j=0;j<n;j++){
                double s=0.0;
                for (int k=0;k<n;k++) s += R[i][k]*Q[k][j];
                A[i][j] = s;
            }
        }

        // reforço de simetria
        for (int i=0;i<n;i++)
            for (int j=i+1;j<n;j++){
                double s = 0.5*(A[i][j]+A[j][i]);
                A[i][j]=A[j][i]=s;
            }

        // parada: norma fora da diagonal
        double off=0.0;
        for (int i=0;i<n;i++) for (int j=0;j<n;j++)
            if (i!=j) off += A[i][j]*A[i][j];
        if (sqrt(off) < tol) break;
    }

    for (int i=0;i<n;i++){ free(Q[i]); free(R[i]); }
    free(Q); free(R);
}
