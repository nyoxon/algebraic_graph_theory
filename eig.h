#ifndef EIG_H
#define EIG_H

double power_method_int
(
  const int** M,
  int n,
  int maxit,
  double tol,
  double* out_vec
);

void qr_iterate
(
  double** A,
  int n,
  int maxit,
  double tol
);

#endif
