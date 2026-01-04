#include <stdio.h>
#include <stdlib.h>

#include "lib_poisson1D.h"          /* for band storage */
#include "lib_poisson1D_sparse.h"   /* CSR / CSC */

/* BLAS band matrix-vector product */
extern void dgbmv_(char*, int*, int*, int*, int*,
                   double*, double*, int*,
                   double*, int*, double*, double*, int*);

                   
int main(void)
{
  int n = 10;
  int ku = 1, kl = 1, kv = 0;
  int lab = ku + kl + 1;

  /* test vector */
  double *x  = malloc(n*sizeof(double));
  double *y_band = malloc(n*sizeof(double));
  double *y_csr  = malloc(n*sizeof(double));
  double *y_csc  = malloc(n*sizeof(double));

  for (int i = 0; i < n; i++)
    x[i] = 1.0;

  /* ---------- Band ---------- */
  double *AB = malloc(lab*n*sizeof(double));
  set_GB_operator_colMajor_poisson1D(AB, &lab, &n, &kv);

  for (int i = 0; i < n; i++)
    y_band[i] = 0.0;

  dgbmv_("N", &n, &n, &kl, &ku,
          &(double){1.0}, AB, &lab,
          x, &(int){1},
          &(double){0.0}, y_band, &(int){1});

  /* ---------- CSR ---------- */
  CSRMatrix Acsr;
  poisson1D_to_CSR(&Acsr, n);
  csr_matvec(&Acsr, x, y_csr);

  /* ---------- CSC ---------- */
  CSCMatrix Acsc;
  poisson1D_to_CSC(&Acsc, n);
  csc_matvec(&Acsc, x, y_csc);

  /* ---------- Compare ---------- */
  printf(" i | band | csr | csc\n");
  printf("----------------------\n");
  for (int i = 0; i < n; i++) {
    printf("%2d | %5.1f | %5.1f | %5.1f\n",
           i, y_band[i], y_csr[i], y_csc[i]);
  }

  /* free */
  free(x);
  free(y_band);
  free(y_csr);
  free(y_csc);
  free(AB);
  free_CSR(&Acsr);
  free_CSC(&Acsc);

  return 0;
}
