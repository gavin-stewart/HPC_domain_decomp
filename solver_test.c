#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "suitesparse/cs.h"
#include "local_solvers.c"

/*
 * Tests the umfpack_solve function against an example listed in the UMFPACK 
 * quickstart guide.  The result should be [1,2,3,4,5].
 */
int main(int argc, char **argv) {
    // Data from UMFPACk_QuickStart
    int n = 5;
    int Ap[] = {0,2,5,9,10,12};
    int Ai[] = {0,1,0,2,4,1,2,3,4,2,1,4};
    double Ax[] = {2.,3.,3.,-1.,4.,4.,-3.,1.,2.,2.,6.,1.};
    double b[] = {8.0,45.0,-3.0,3.0,19.0};
    double x[5];
    int i;
    cs_di A;

    //Define the parameters of A
    A.nzmax = 12;
    A.m = n;
    A.n = n;
    A.p = Ap;
    A.i = Ai;
    A.x = Ax;
    A.nz = -1;

    umfpack_solve(&A, (double*)x, (double*)b);
    for(i = 0; i < n; i++) {
        assert(fabs(x[i] - i - 1) < 1e-5);
    }

    printf("Local solver test: PASSED\n");

    return EXIT_SUCCESS;
}
