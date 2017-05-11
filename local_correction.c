#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <suitesparse/cs.h>
#include "local_solvers.c"

#define DEBUG

void set_up_indices(int*, int*, int*);
double *extract_rect_subregion(double*, int, int, int, int, int, int);
void include_rect_subregion(double*, double*, int, int, int, int, int, int);

int compute_local_correction(double *lu, double *lresid, int lN, int overlap, int ppe, cs_di **laplace_mat) {
  double *lcorrection;
  double *resid_w_overlap;
  double h = 1.0 / (lN * ppe + 1);
  double invhsq = 1.0 / (h * h);
  int nx = lN;
  int ny = lN;
  int mpi_rank, p_row, p_col;
  int i, j;
  /*
   * These two pointers are indexed as
   * 5 |   6   | 7
   *---------------
   *   |       |
   * 3 |       | 4
   *   |       |
   *---------------
   * 0 |   1   | 2
   * where the center region is the local region and the numbered regions are
   * the nonlocal regions for the residual with overlap.
   */
  double **send_subregions = calloc(8, sizeof(double*));
  double **recv_subregions = calloc(8, sizeof(double*));
  int num_comms = 0;
  MPI_Request *send_request = calloc(8, sizeof(MPI_Request));
  MPI_Request *recv_request = calloc(8, sizeof(MPI_Request));
  int horiz_overlap[] = {0,lN,0};
  int vert_overlap[] = {0,lN,0};
  int *horiz_start_ind = malloc(3 * sizeof(int));
  int *horiz_end_ind = malloc(3 * sizeof(int));
  int *vert_start_ind = malloc(3 * sizeof(int));
  int *vert_end_ind = malloc(3 * sizeof(int));

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  p_row = mpi_rank / ppe;
  p_col = mpi_rank % ppe;

  // Determine how large the array for residuals (with overlap) should be
  if(p_col > 0) {
    horiz_overlap[0] = overlap;
    nx += overlap;
  }
  if(p_col < ppe - 1) {
    horiz_overlap[2] = overlap;
    nx += overlap;
  }
  if(p_row > 0) {
    vert_overlap[0] = overlap;
    ny += overlap;
  }
  if(p_row < ppe - 1) {
    vert_overlap[2] = overlap;
    ny += overlap;
  }

  //Set up start, end indices
  set_up_indices(horiz_overlap, horiz_start_ind, horiz_end_ind);
  set_up_indices(vert_overlap, vert_start_ind, vert_end_ind);

  num_comms = 0;
  //Set up send, recv subregions and begin communication.
  for(j = 0; j < 3; j++) {
    for(i = 0; i < 3; i++) {
      int index = 3 * j + i;
      if (index == 4) { //Skip the center cell.
        continue;
      } else if(index > 4) {
        index--;
      }
      if(vert_overlap[j] > 0 && horiz_overlap[i] > 0) { 
        //We have a processor to send to
        int comm_rank = (p_row + j - 1) * ppe + (p_col + i - 1);
        int data_size = vert_overlap[j] * horiz_overlap[i];
        recv_subregions[index] = calloc(data_size, sizeof(double));
        MPI_Irecv(recv_subregions[index], data_size, MPI_DOUBLE, comm_rank,
            987, MPI_COMM_WORLD, &recv_request[num_comms]);
        send_subregions[index] = extract_rect_subregion(lresid, lN, lN, 
            vert_start_ind[j], vert_end_ind[j], 
            horiz_start_ind[i], horiz_end_ind[i]);
        MPI_Isend(send_subregions[index], data_size, MPI_DOUBLE, comm_rank,
            987, MPI_COMM_WORLD, &send_request[num_comms]);
        num_comms++;
      }
    }
  }

  if(*laplace_mat == NULL) { 
    //Initialize the Laplace matrix here to be used in the future
    *laplace_mat = construct_laplace_mat(nx, ny, invhsq);
  }

  // Allocate the arrays for the overlapping grid
  lcorrection = malloc(nx * ny * sizeof(double));
  resid_w_overlap = malloc(nx * ny * sizeof(double));

  //Include local data in the residual
  include_rect_subregion(resid_w_overlap, lresid, nx, ny, vert_start_ind[1],
      vert_end_ind[1], horiz_start_ind[1], horiz_end_ind[1]);

  // Wait for communications to conclude
  MPI_Waitall(num_comms, recv_request, MPI_STATUSES_IGNORE);

  // Include the non-local residuals
  int turn;
  MPI_Barrier(MPI_COMM_WORLD);
  for(turn = 0; turn < ppe * ppe; turn++) {
    if(mpi_rank == turn) {
      for(j = 0; j < 3; j++) {
        for(i = 0; i < 3; i++) {
          int index = 3 * j + i;
          if (index == 4) { //Skip the center cell.
            continue;
          } else if(index > 4) {
            index--;
          }
          if(vert_overlap[j] > 0 && horiz_overlap[i] > 0) { 
            //We received some data; incorporate it
            include_rect_subregion(resid_w_overlap, recv_subregions[index], 
                nx, ny, 
                vert_start_ind[j], vert_end_ind[j],
                horiz_start_ind[i], horiz_end_ind[i]);
            printf("Proc %d including r%d: %d -- %d x %d -- %d\n", 
                mpi_rank, index,
                vert_start_ind[j], vert_end_ind[j],
                horiz_start_ind[i], horiz_end_ind[i]);
#ifdef DEBUG
            printf("Residual for mpi_rank %d\n", mpi_rank);
            int r,c;
            for(r = ny - 1; r >= 0; r--) {
              for(c = 0; c < nx; c++) {
                printf("%f ", resid_w_overlap[r * nx + c]);
              }
              printf("\n");
            }
#endif
          }
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  //See what the residual looks like

  // Compute the local correction
  umfpack_solve(*laplace_mat, lcorrection, resid_w_overlap);
  // Add the local correction to the local solution
  for(j = 0; j < lN; j++) {
    //TODO: Switch to memcpy?
    for(i = 0; i < lN; i++) {
      int local_ind = j * lN + i;
      int overlap_ind = (j + vert_overlap[0]) * nx + (i + horiz_overlap[0]);
      lu[local_ind] += lcorrection[overlap_ind];
    }
  }
  // Free allocated memory (except for send-related memory)
  free(lcorrection);
  free(resid_w_overlap);
  for(i = 0; i < 8; i++) {
    if(recv_subregions[i])
      free(recv_subregions[i]);
  }
  free(recv_subregions);
  free(horiz_start_ind);
  free(horiz_end_ind);
  free(vert_start_ind);
  free(vert_end_ind);
  free(recv_request);
  // Wait for all sends to go through
  MPI_Waitall(num_comms, send_request, MPI_STATUSES_IGNORE);
  // Free send-related memory
  free(send_request);
  for(i = 0; i < 8; i++) {
    if(send_subregions[i])
      free(send_subregions[i]);
  }
  free(send_subregions);

}

void set_up_indices(int *edge_sizes, int *start_inds, int *end_inds) {
  int i;
  start_inds[0] = 0;
  end_inds[0] = edge_sizes[0];
  for(i = 1; i < 3; i++) {
    start_inds[i] = start_inds[i-1] + edge_sizes[i-1];
    end_inds[i] = end_inds[i-1] + edge_sizes[i];
  }
}

double *extract_rect_subregion(double *region, int nx, int ny, int row_start, int row_end, int col_start, int col_end) {
  int row, col;
  int sub_ny = row_end - row_start;
  int sub_nx = col_end - col_start;
  double *subregion = malloc(sizeof(double) * sub_nx * sub_ny);
  for(row = row_start; row < row_end; row++) {
    int sub_row = row - row_start;
    for(col = col_start; col < col_end; col++) {
      int sub_col = col - col_start;
      int sub_ind = sub_row * sub_nx + sub_col;
      int ind = row * nx + col;
      subregion[sub_ind] = region[ind];
    }
  }
  return subregion;
}

void include_rect_subregion(double *region, double *subregion, int nx, int ny, 
    int row_start, int row_end, int col_start, int col_end) {
  int sub_ny = row_end - row_start;
  int sub_nx = col_end - col_start;
  int row, sub_row, col, sub_col;

  for(row = row_start; row < row_end; row++) {
    sub_row = row - row_start;
    for(col = col_start; col < col_end; col++) {
      int ind = row * nx + col;
      sub_col = col - col_start;
      int sub_ind = sub_row * sub_nx + sub_col;
      region[ind] = subregion[sub_ind];
    }
  }
}
