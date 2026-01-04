/**********************************************/
/* lib_poisson1D_richardson.c                 */
/* Richardson and preconditioned Richardson  */
/* methods for Poisson 1D                     */
/**********************************************/

#include "lib_poisson1D.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* BLAS interfaces (Fortran) */
extern double dnrm2_(int*, double*, int*);
extern void   dgbmv_(char*, int*, int*, int*, int*,
                     double*, double*, int*,
                     double*, int*, double*, double*, int*);
extern void   dcopy_(int*, double*, int*, double*, int*);
extern void   daxpy_(int*, double*, double*, int*, double*, int*);

/****************************************************/
/* ex7 : Spectral analysis of Poisson 1D operator   */
/****************************************************/

/* Compute all eigenvalues of the 1D Poisson operator */
void eig_poisson1D(double* eigval, int *la)
{
  int n = *la;

  for (int j = 1; j <= n; j++) {
    eigval[j-1] = 2.0 * (1.0 - cos(j * M_PI / (n + 1)));
  }
}

/* Return minimum eigenvalue */
double eigmin_poisson1D(int *la)
{
  int n = *la;
  double *eigval = (double*) malloc(n * sizeof(double));

  eig_poisson1D(eigval, la);

  double lmin = eigval[0];
  for (int i = 1; i < n; i++) {
    if (eigval[i] < lmin) lmin = eigval[i];
  }

  free(eigval);
  return lmin;
}

/* Return maximum eigenvalue */
double eigmax_poisson1D(int *la)
{
  int n = *la;
  double *eigval = (double*) malloc(n * sizeof(double));

  eig_poisson1D(eigval, la);

  double lmax = eigval[0];
  for (int i = 1; i < n; i++) {
    if (eigval[i] > lmax) lmax = eigval[i];
  }

  free(eigval);
  return lmax;
}

/* Optimal Richardson relaxation parameter */
double richardson_alpha_opt(int *la)
{
  double lmin = eigmin_poisson1D(la);
  double lmax = eigmax_poisson1D(la);

  return 2.0 / (lmin + lmax);
}

/****************************************************/
/* ex7 : Richardson iteration with optimal alpha    */
/****************************************************/

void richardson_alpha(double *AB, double *RHS, double *X,
                      double *alpha_rich,
                      int *lab, int *la, int *ku, int *kl,
                      double *tol, int *maxit,
                      double *resvec, int *nbite)
{
  int n = *la;
  int inc = 1;
  double alpha = *alpha_rich;

  double *AX = (double*) malloc(n * sizeof(double));
  double *R  = (double*) malloc(n * sizeof(double));

  double one = 1.0, zero = 0.0, minus = -1.0;

  double normb = dnrm2_(&n, RHS, &inc);
  if (normb == 0.0) normb = 1.0;

  *nbite = 0;

  for (int k = 0; k < *maxit; k++) {

    /* AX = A * X */
    dgbmv_("N", &n, &n, kl, ku,
            &one, AB, lab,
            X, &inc,
            &zero, AX, &inc);

    /* R = RHS - AX */
    dcopy_(&n, RHS, &inc, R, &inc);
    daxpy_(&n, &minus, AX, &inc, R, &inc);

    /* Relative residual */
    double res = dnrm2_(&n, R, &inc) / normb;
    resvec[k] = res;
    *nbite = k + 1;

    if (res < *tol) break;

    /* X = X + alpha * R */
    daxpy_(&n, &alpha, R, &inc, X, &inc);
  }

  free(AX);
  free(R);
}

/****************************************************/
/* ex8 / ex9 : Preconditioners                      */
/****************************************************/

/* Jacobi preconditioner: MB = D */
void extract_MB_jacobi_tridiag(double *AB, double *MB,
                               int *lab, int *la,
                               int *ku, int *kl, int *kv)
{
  int n = *la;

  /* Initialize MB to zero */
  for (int j = 0; j < n; j++)
    for (int i = 0; i < *lab; i++)
      MB[indexABCol(i, j, lab)] = 0.0;

  /* Copy diagonal */
  for (int j = 0; j < n; j++) {
    MB[indexABCol(*ku, j, lab)] =
        AB[indexABCol(*ku, j, lab)];
  }
}

/* Gauss–Seidel preconditioner: MB = D - E */
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB,
                                     int *lab, int *la,
                                     int *ku, int *kl, int *kv)
{
  int n = *la;

  /* Initialize MB to zero */
  for (int j = 0; j < n; j++)
    for (int i = 0; i < *lab; i++)
      MB[indexABCol(i, j, lab)] = 0.0;

  /* Copy diagonal and lower diagonal */
  for (int j = 0; j < n; j++) {
    MB[indexABCol(*ku, j, lab)] =
        AB[indexABCol(*ku, j, lab)];

    if (j < n-1) {
      MB[indexABCol(*ku+1, j, lab)] =
          AB[indexABCol(*ku+1, j, lab)];
    }
  }
}

/****************************************************/
/* ex8 / ex9 : Preconditioned Richardson             */
/****************************************************/

void richardson_MB(double *AB, double *RHS, double *X, double *MB,
                   int *lab, int *la, int *ku, int *kl,
                   double *tol, int *maxit,
                   double *resvec, int *nbite)
{
  int n = *la;
  int inc = 1;
  double one = 1.0, zero = 0.0, minus = -1.0;

  double *AX = (double*) malloc(n * sizeof(double));
  double *R  = (double*) malloc(n * sizeof(double));
  double *Z  = (double*) malloc(n * sizeof(double));

  double normb = dnrm2_(&n, RHS, &inc);
  if (normb == 0.0) normb = 1.0;

  *nbite = 0;

  for (int k = 0; k < *maxit; k++) {

    /* AX = A * X */
    dgbmv_("N", &n, &n, kl, ku,
            &one, AB, lab,
            X, &inc,
            &zero, AX, &inc);

    /* R = RHS - AX */
    dcopy_(&n, RHS, &inc, R, &inc);
    daxpy_(&n, &minus, AX, &inc, R, &inc);

    /* ---- Solve M Z = R ---- */
    for (int i = 0; i < n; i++) {

      double sum = R[i];

      /* GS only: subtract lower-diagonal contribution */
      if (i > 0 && MB[indexABCol(*ku+1, i-1, lab)] != 0.0) {
        sum -= MB[indexABCol(*ku+1, i-1, lab)] * Z[i-1];
      }

      Z[i] = sum / MB[indexABCol(*ku, i, lab)];
    }

    /* X = X + Z */
    daxpy_(&n, &one, Z, &inc, X, &inc);

    double res = dnrm2_(&n, R, &inc) / normb;
    resvec[k] = res;
    *nbite = k + 1;

    if (res < *tol) break;
  }

  free(AX);
  free(R);
  free(Z);
}
