/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int n    = *la;         
  int ldab = *lab;         

  double h    = 1.0 / (n + 1);
  double diag =  2.0 / (h * h);
  double off  = -1.0 / (h * h);

  for (int j = 0; j < n; j++) {
    for (int i = 0; i < ldab; i++) {
      AB[indexABCol(i, j, lab)] = 0.0;
    }
  }

  int row_extra  = 0; 
  (void)row_extra;   
  int row_super  = 1;  
  int row_diag   = 2;  
  int row_sub    = 3;  

  for (int j = 0; j < n; j++) {
    if (j > 0) {
      AB[indexABCol(row_super, j, lab)] = off;
    }

    AB[indexABCol(row_diag, j, lab)] = diag;

    if (j < n - 1) {
      AB[indexABCol(row_sub, j, lab)] = off;
    }
  }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int i, j;
  for (j = 0; j < *la; j++){
    for (i = 0; i < *lab; i++){
      AB[indexABCol(i, j, lab)] = 0.0;
    }
    AB[indexABCol(*kv, j, lab)] = 1.0;
  }
}


void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  double h = 1.0/((*la)+1);
  double h2 = h*h;

  for (int i = 0; i < *la; i++){
    RHS[i] = 0.0;
  }

  RHS[0]     += (*BC0)/h2;
  RHS[*la-1] += (*BC1)/h2;
}


void set_analytical_solution_DBC_1D(double* EX_SOL, double* X,
                                    int* la, double* BC0, double* BC1){
  for (int i = 0; i < *la; i++){
    EX_SOL[i] = (*BC0) + X[i]*((*BC1)-(*BC0));
  }
}


void set_grid_points_1D(double* x, int* la){
  double h = 1.0/((*la)+1);
  for (int i = 0; i < *la; i++){
    x[i] = (i+1)*h;
  }
}


double relative_forward_error(double* x, double* y, int* la){
  double num = 0.0, den = 0.0;
  for (int i = 0; i < *la; i++){
    num += (x[i]-y[i])*(x[i]-y[i]);
    den += y[i]*y[i];
  }
  return sqrt(num)/sqrt(den);
}


int indexABCol(int i, int j, int *lab){
  return j*(*lab) + i;
}


int dgbtrftridiag(int *la, int *n, int *kl, int *ku,
                  double *AB, int *lab,
                  int *ipiv, int *info)
{
  int j;
  *info = 0;

  /* Explicit GB row indices */
  int row_super = 1;
  int row_diag  = 2;
  int row_sub   = 3;

  for (j = 0; j < *la - 1; j++) {

    double pivot = AB[indexABCol(row_diag, j, lab)];
    if (pivot == 0.0) {
      *info = j + 1;
      return *info;
    }

    /* L(j+1,j) */
    AB[indexABCol(row_sub, j, lab)] /= pivot;

    /* U(j+1,j+1) */
    AB[indexABCol(row_diag, j + 1, lab)] -=
      AB[indexABCol(row_sub, j, lab)] *
      AB[indexABCol(row_super, j + 1, lab)];
  }

  return *info;
}

