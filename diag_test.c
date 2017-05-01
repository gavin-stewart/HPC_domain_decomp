/*
 *  Code to allow visual inspection of the construct_sparse_diags function.
 *  TODO: This does not currently work.  Make it work.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "local_solvers.c"
#include "suitesparse/cs.h"

int main(int argc, char** argv) {
    cs_di *A;
    int n = 5;
    int num_diags = 1;
    int diag_numbers[] = {0};
    double **diags;
    double *main_diag;// = {1,2,3,4,5};
    int i;
    main_diag = malloc(5 * sizeof(double));
    for(i = 0; i < n; i++) {
        main_diag[i] = i+1;
    }
    diags = &main_diag;

    A = construct_sparse_diag(n, num_diags, diag_numbers, diags);


    for(i = 0; i < n; i++) {
        assert(A->p[i] == i);
        assert(A->i[i] == i);
        assert(fabs(A->x[i] - (i+1)) < 1e-5);
    }
    assert(A->p[n] == n);

    printf("Diag_test: PASSED\n");
    
    return EXIT_SUCCESS;
}
