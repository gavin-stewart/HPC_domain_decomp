/*
 * Solve div grad u = sin(pi*x)sin(2pi*y) in the square
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <suitesparse/cs.h>
#include "local_solvers.c"


int main(int argc, char **argv) {
    int n = 400;
    double h = 1.0 / (n+1);
    double *f = malloc(n * n * sizeof(double));
    int i,j;
    double *u = malloc(n * n * sizeof(double));
    double *u_true = malloc(n * n * sizeof(double));
    double **diags = malloc(5 * sizeof(double*));
    double *center_diag, *vert_diag, *left_diag, *right_diag;
    int num_diags = 5;
    int *diag_ind = malloc(5 * sizeof(int));
    double max_err = 0;
    center_diag = malloc(n * n * sizeof(double));
    left_diag = malloc(n * n * sizeof(double));
    right_diag = malloc(n * n * sizeof(double));
    vert_diag = malloc(n * n * sizeof(double));

    cs_di *A;

    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            int ind = i + n * j;
            //f = sin(pi*x)sin(2pi*y)
            f[ind] = sin(M_PI * h * i) * sin(2 * h * M_PI *j); 
            u_true[ind] = -f[ind] / (5 * M_PI * M_PI);
        }
    }

    for(j = 0; j < n; j++) {
        for(i = 0; i < n; i++) {
            int ind = i + n*j;
            center_diag[ind] = -4 / (h * h);
            vert_diag[ind] = 1 / (h * h);
            if(ind > 0)
                left_diag[ind-1] = (i==0)?(0):(1 / (h*h));
            if(ind < n*n - 1)
                right_diag[ind+1] = (i==(n-1))?(0):(1 / (h*h));
        }
    }

    diag_ind[0] = -n;
    diag_ind[1] = -1;
    diag_ind[2] = 0;
    diag_ind[3] = 1;
    diag_ind[4] = n;
    diags[0] = vert_diag;
    diags[1] = left_diag;
    diags[2] = center_diag;
    diags[3] = right_diag;
    diags[4] = vert_diag;

    A = construct_sparse_diag(n*n, num_diags, diag_ind, diags);

    //cs_print(A, 0);
    
    umfpack_solve(A, u, f);
    
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            int ind = i + n * j;
            double err = fabs(u[ind] - u_true[ind]);
            if(err > max_err) {
                max_err = err;
            }
        }
    }
    printf("The maximum error was %f\n", max_err);
    return EXIT_SUCCESS;
}
