#include <stdlib.h>
#include "lib_poisson1D_sparse.h"

/* ---------- CSR ---------- */

void poisson1D_to_CSR(CSRMatrix *A, int n)
{
  A->n = n;
  A->nnz = 3*n - 2;

  A->rowptr = malloc((n+1)*sizeof(int));
  A->colind = malloc(A->nnz*sizeof(int));
  A->val    = malloc(A->nnz*sizeof(double));

  int k = 0;
  A->rowptr[0] = 0;

  for (int i = 0; i < n; i++) {
    if (i > 0) {
      A->colind[k] = i-1;
      A->val[k++]  = -1.0;
    }
    A->colind[k] = i;
    A->val[k++]  = 2.0;
    if (i < n-1) {
      A->colind[k] = i+1;
      A->val[k++]  = -1.0;
    }
    A->rowptr[i+1] = k;
  }
}

void csr_matvec(const CSRMatrix *A, const double *x, double *y)
{
  for (int i = 0; i < A->n; i++) {
    y[i] = 0.0;
    for (int k = A->rowptr[i]; k < A->rowptr[i+1]; k++) {
      y[i] += A->val[k] * x[A->colind[k]];
    }
  }
}

void free_CSR(CSRMatrix *A)
{
  free(A->rowptr);
  free(A->colind);
  free(A->val);
}

/* ---------- CSC ---------- */

void poisson1D_to_CSC(CSCMatrix *A, int n)
{
  A->n = n;
  A->nnz = 3*n - 2;

  A->colptr = malloc((n+1)*sizeof(int));
  A->rowind = malloc(A->nnz*sizeof(int));
  A->val    = malloc(A->nnz*sizeof(double));

  int k = 0;
  A->colptr[0] = 0;

  for (int j = 0; j < n; j++) {
    if (j > 0) {
      A->rowind[k] = j-1;
      A->val[k++]  = -1.0;
    }
    A->rowind[k] = j;
    A->val[k++]  = 2.0;
    if (j < n-1) {
      A->rowind[k] = j+1;
      A->val[k++]  = -1.0;
    }
    A->colptr[j+1] = k;
  }
}

void csc_matvec(const CSCMatrix *A, const double *x, double *y)
{
  for (int i = 0; i < A->n; i++)
    y[i] = 0.0;

  for (int j = 0; j < A->n; j++) {
    for (int k = A->colptr[j]; k < A->colptr[j+1]; k++) {
      y[A->rowind[k]] += A->val[k] * x[j];
    }
  }
}

void free_CSC(CSCMatrix *A)
{
  free(A->colptr);
  free(A->rowind);
  free(A->val);
}
