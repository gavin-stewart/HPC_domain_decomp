#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#include <suitesparse/cholmod.h>
#include <suitesparse/umfpack.h>
#include <suitesparse/cs.h>

#define SUCCESS_RETURN (0)
#define FAIL_RETURN (1)

//cholmod_sparse *cholmod_sparse_from_cxsparse(cs_di *, cholmod_common *);

/*
 * Solve the system Ax = b, store the result in x.
 */
int umfpack_solve(cs_di *A, double *x, double *b) {
  if(A->m != A->n) { //If A is not square
    return FAIL_RETURN;
  }
  int n = A->n; //Equal to m by check.
  int *Ap = A->p;
  int *Ai = A->i;
  double *Ax = A->x;
  void *Symbolic = NULL;
  void *Numeric = NULL;

  int symb_status = umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, 
      (double*)NULL, (double*)NULL);

  if(!Symbolic) {
    printf("Symbolic was not assigned, status code %d\n", symb_status);
    return FAIL_RETURN;
  }

  int num_status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, 
      (double*)NULL, (double*)NULL);

  if(!Numeric) {
    printf("Numeric was not assigned, status code %d\n", num_status);
    return FAIL_RETURN;
  }

  umfpack_di_free_symbolic(&Symbolic);

  return umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, 
      (double*)NULL, (double*)NULL);

  //umfpack_di_free_numeric(&Numeric);
}

int cmp(const void *a, const void *b) {
  return (*(int*)a - *(int*)b);
}

cs_di *construct_sparse_diag(int n, int num_diags, int *diag_numbers, double **diags) {
  cs_di *A = malloc(sizeof(cs_di));
  int col, row, diag_ind;
  int *Ap, *Ai;
  double *Ax;
  int cum_elems = 0;
  double curr_val;

  Ap = calloc((n+1), sizeof(int));
  Ai = calloc(n * num_diags, sizeof(int)); //May be shrunk at the end
  Ax = calloc(n * num_diags, sizeof(double)); //May be shrunk at the end

  Ap[0] = 0;
  qsort(diag_numbers, num_diags, sizeof(int), cmp);

  for(col = 0; col < n; col++) {
    Ap[col+1] = Ap[col];
    for(diag_ind = 0; diag_ind < num_diags; diag_ind++) {
      row = col + diag_numbers[diag_ind];
      if(row < 0 || row >= n) { //This diagonal hasn't started/has ended
        continue;
      }
      curr_val = diags[diag_ind][row];
      if(curr_val != 0) {
        Ax[cum_elems] = curr_val;
        Ap[col+1] = Ap[col+1] + 1;
        Ai[cum_elems] = row;
        cum_elems++;
      }
    }
  }

  //Resize the arrays to have the correct number of elements
  Ai = realloc(Ai, cum_elems * sizeof(int));
  Ax = realloc(Ax, cum_elems * sizeof(double));

  A->n = n;
  A->m = n;
  A->p = Ap;
  A->i = Ai;
  A->x = Ax;
  A->nz = -1;  //Compressed column

  return A;
}

cs_di *construct_sqr_laplace_mat(int n, double invhsq) {
  int n_sqr = n * n;
  int row,col, ind;
  int num_diags = 5; //Quindiagonal
  int *diag_numbers = malloc(num_diags * sizeof(int));
  double **diags = malloc(num_diags * sizeof(double*));

  double *up_diag = malloc(n_sqr * sizeof(double));
  double *left_diag = malloc(n_sqr * sizeof(double));
  double *right_diag = malloc(n_sqr * sizeof(double));
  double *down_diag = malloc(n_sqr * sizeof(double));
  double *center_diag = malloc(n_sqr * sizeof(double));

  for(row = 0; row < n; row++) {
    for(col = 0; col < n; col++) {
      ind = row + n * col;
      center_diag[ind] = 4 * invhsq;
      if(col > 0) {
        down_diag[ind - n] = -invhsq;
      }
      if(col < n -1) {
        up_diag[ind + n] = -invhsq;
      }
      if(row > 0) {
        left_diag[ind - 1] = -invhsq;
      }
      if(row < n - 1) {
        right_diag[ind + 1] = -invhsq;
      }
    }
  }
  diag_numbers[0] = -n;
  diag_numbers[1] = -1;
  diag_numbers[2] = 0;
  diag_numbers[3] = 1;
  diag_numbers[4] = n;

  diags[0] = down_diag;
  diags[1] = left_diag;
  diags[2] = center_diag;
  diags[3] = right_diag;
  diags[4] = up_diag;

  return construct_sparse_diag(n_sqr, num_diags, diag_numbers, diags);
}

/*
int cholmod_factor_and_solve(cs_di *A, double *x, double *b) {
  //
  //Convert the cs_di matrix to cholmod form
  cholmod_sparse *A_chol;
  cholmod_dense *x_chol, *b_chol;
  cholmod_common c;
  cholmod_factor *L;
  int n = A->n;

  cholmod_start(&c);
  A_chol = cholmod_sparse_from_cxsparse(A, &c);
  
  //Allocate b_chol
  b_chol = cholmod_allocate_dense(n, 1, n, CHOLMOD_REAL, &c);
  memcpy(b_chol->x, b, n*sizeof(double));
  cholmod_write_dense(stdout, b_chol, "", &c);
  //Factorize A_chol
  L = cholmod_analyze(A_chol, &c);
  cholmod_factorize(A_chol, L, &c);
  //Solve for x_chol
  x_chol = cholmod_solve(CHOLMOD_A, L, b_chol, &c);
  cholmod_write_dense(stdout, x_chol, "", &c);
  //Set x to array from x_chol
  memcpy(x, x_chol->x, n*sizeof(double));
  //Free memory from L, b_chol, A_chol, x_chol (but such that 
  //A, x, and b are not affected).
  //TODO: Make sure A, x, b are not de-allocated
  cholmod_free_factor(&L, &c);
  //cholmod_free_sparse(&A_chol, &c);
  //cholmod_free_dense(&x_chol, &c);
  //cholmod_free_dense(&b_chol, &c);
  cholmod_finish(&c);
  return EXIT_SUCCESS;
}

cholmod_sparse *cholmod_sparse_from_cxsparse(cs_di *A, cholmod_common *c) {
  cholmod_sparse *B;
  B = cholmod_allocate_sparse(A->n, A->m, A->p[A->n], 1, 1, 0, CHOLMOD_REAL, c);
  B->p = A->p;
  B->i = A->i;
  B->x = A->x;
  //TEST
  cholmod_write_sparse(stdout, B, NULL, "", c);
  return B;
}*/
