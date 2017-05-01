#include <stdio.h>
#include <stdlib.h>

#include <suitesparse/umfpack.h>
#include "suitesparse/cs.h"

#define SUCCESS_RETURN (0)
#define FAIL_RETURN (1)

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
            if(row < 0 || row > n) { //This diagonal hasn't started/has ended
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

int cholmod_solve(cs_di *A, double *x, double *b) {
    //TODO: Implement
    return EXIT_FAILURE;
}
