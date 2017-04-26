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


void include_bcs(double *b, double h, double *left_bdy, double *right_bdy, 
        double *top_bdy, double *bottom_bdy) {
    //TODO: Implement
}
