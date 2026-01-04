#ifndef LIB_POISSON1D_SPARSE_H
#define LIB_POISSON1D_SPARSE_H

/* ---------- CSR ---------- */
typedef struct {
  int n;
  int nnz;
  int *rowptr;
  int *colind;
  double *val;
} CSRMatrix;

void poisson1D_to_CSR(CSRMatrix *A, int n);
void csr_matvec(const CSRMatrix *A, const double *x, double *y);
void free_CSR(CSRMatrix *A);

/* ---------- CSC ---------- */
typedef struct {
  int n;
  int nnz;
  int *colptr;
  int *rowind;
  double *val;
} CSCMatrix;

void poisson1D_to_CSC(CSCMatrix *A, int n);
void csc_matvec(const CSCMatrix *A, const double *x, double *y);
void free_CSC(CSCMatrix *A);

#endif
