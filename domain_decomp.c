#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "local_solvers.c"
#include "residue.c"

#define TOL 1e-4
#define DEBUG
int main(int argc, char **argv) {
  int iter, max_iter, num_proc, ppe, lN, N, mpi_rank, i;
  double initial_resid_norm, resid_norm, h, invhsq;
  double *lresidual;
  double *lu;
  double *f;
  double *lcorrection;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Check for the right number of arguments
  if(argc != 3) {
    if(mpi_rank == 0) {
      fprintf(stderr, "Usage: %s <# local grid pts> <max iterations>\n", argv[0]);
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  
  lN = atoi(argv[1]);
  max_iter = atoi(argv[2]);

  lu = calloc(lN * lN, sizeof(double));
  f = calloc(lN * lN, sizeof(double));
  lcorrection = malloc(lN * lN * sizeof(double));

  //TODO: Set f to something 'interesting'
  for(i = 0; i < lN * lN; i++) {
    f[i] = 1;
  }
  // Compute the number of processors per edge, check that the number of 
  // threads is a square.
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  ppe = (int) floor(sqrt(num_proc));
  if(ppe * ppe != num_proc) {
    if(mpi_rank == 0) 
      fprintf(stderr, "Error: number of processors must be a square!\n");
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  
  N = ppe * lN;
  h = 1.0 / (N + 1);
  invhsq = 1.0 / (h * h);

  // Form the Laplace matrix
  cs_di *laplace_mat = construct_sqr_laplace_mat(lN, invhsq);


  // Compute the initial residual
  lresidual = residue(lu, f, lN, num_proc, invhsq, mpi_rank);
  initial_resid_norm = compute_residual(lresidual, lN, invhsq);
  resid_norm = initial_resid_norm;
  #ifdef DEBUG
  /*printf("Local residual: \n");
  for(i = 0; i < lN; i++) {
    int j;
    for(j = 0;  j < lN; j++) {
      printf("%3f ", lresidual[j + lN * i]);
    }
    printf("\n");
  }*/
  if(mpi_rank == 0)
    printf("Initial residual: %f\n", initial_resid_norm);
  #endif
  for(iter = 0; iter < max_iter && resid_norm > initial_resid_norm * TOL; iter++) {
    // Solve the local system -(Laplace)correction = residual
    umfpack_solve(laplace_mat, lcorrection, lresidual);
    // Add the correction to lu
    // TODO: Could move to function for clarity
    for(i = 0; i < lN * lN; i++) {
      lu[i] += lcorrection[i];
    }
    // Compute the new lresidual, resid_norm
    free(lresidual);
    lresidual = residue(lu, f, lN, num_proc, invhsq,  mpi_rank);
    /*int turn;
    for(turn = 0; turn < num_proc; turn++) {
      if(mpi_rank == turn) {
        printf("Residual for rank %d:\n", mpi_rank);
        for(i = 0; i < lN; i++) {
          int j;
          for(j = 0;  j < lN; j++) {
            printf("%3f ", lresidual[j + lN * i]);
          }
          printf("\n");
        }
        printf("\n\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }*/
    resid_norm = compute_residual(lresidual, lN, invhsq);
    #ifdef DEBUG
    if(mpi_rank == 0)
      printf("Iteration [%d]: residual norm: %f\n", iter, resid_norm);
    #endif
  }
  /*int turn;
  for(turn = 0; turn < num_proc; turn++) {
    if(mpi_rank==turn) {
      int row, col;
      printf("Residual for processor %d:\n\n", mpi_rank);
      for(row = 0; row < lN; row++) {
        for(col = 0;  col < lN; col++) {
          printf("%4f ", lu[row * lN + col]);
        }
        printf("\n");
      }
      printf("\n\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }*/
  free(lresidual);
  MPI_Finalize();
}
