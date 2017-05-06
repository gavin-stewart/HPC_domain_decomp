#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "util.h"
#include <string.h>

/* compuate global residual, assuming ghost values are updated */
double compute_residual(double *lu, int lN, double invhsq)//lu is the r and invshq is the
{
    int i;
    double  gres = 0.0, lres = 0.0;
    
    for (i = 0; i <(lN)*(lN); ++i) {
            lres=lres+pow(lu[i],2);
        }
                          
        /* use allreduce for convenience; a reduce would also be sufficient */
        MPI_Allreduce(&lres, &gres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return sqrt(gres);
        }

double * residue(double* lu, double* f, int lN, int p,double invhsq,int mpirank)
{
    int i;//N is total no of rows and colums
    MPI_Status status, status1;

    
    /* Allocation of vectors, including left/upper and right/lower ghost points */
    double * lunew = (double *) calloc( lN * lN, sizeof(double));
    double * ta = (double *) calloc(lN, sizeof(double));
    double * ba = (double *) calloc(lN, sizeof(double));
    double * la = (double *) calloc(lN, sizeof(double));
    double * ra = (double *) calloc(lN, sizeof(double));


    
    double hsq = 1./invhsq;
    
    /* communicate ghost values */
    int np=sqrt(p);
    int mr=mpirank/(np);
    int mc=mpirank-mr*(np);
   
    if (mr <(np-1)) {
        int srrank=(mr+1)*np+mc;
        /* If not a top row process, send/recv bdry values to the top */
        for(i=0;i<(lN);i++)
        {
            MPI_Send(&(lu[(lN-1)*lN+i]), 1, MPI_DOUBLE, srrank, 124, MPI_COMM_WORLD);
            MPI_Recv(&(ta[i]), 1, MPI_DOUBLE, srrank, 123, MPI_COMM_WORLD, &status);
        }
    }
    if (mr > 0) {
        int srrank=(mr-1)*np+mc;
        /* If not a bottom row process, send/recv bdry values to the bottom */
        for(i=0;i<(lN);i++)
        {
            MPI_Send(&(lu[i]), 1, MPI_DOUBLE, srrank, 123, MPI_COMM_WORLD);
            MPI_Recv(&(ba[i]), 1, MPI_DOUBLE, srrank, 124, MPI_COMM_WORLD, &status);
            
        }
    }
    if (mc <(np-1)) {
        int srrank=mr*np+mc+1;
        /* If not a right-most column process, send/recv bdry values to the right */
        for(i=0;i<(lN);i++)
        {
            MPI_Send(&(lu[i*(lN)+lN-1]), 1, MPI_DOUBLE, srrank, 128, MPI_COMM_WORLD);
            MPI_Recv(&(ra[i]), 1, MPI_DOUBLE, srrank, 129, MPI_COMM_WORLD, &status);
            
        }
    }
    if (mc >0) {
        int srrank=mr*np+mc-1;
        /* If not a left-most column process, send/recv bdry values to the left */
        for(i=0;i<(lN);i++)
        {   
            MPI_Send(&(lu[i*(lN)]), 1, MPI_DOUBLE, srrank, 129, MPI_COMM_WORLD);
            MPI_Recv(&(la[i]), 1, MPI_DOUBLE, srrank, 128, MPI_COMM_WORLD, &status);
            
        }
    }
 
        /* residue calculation */
        for (i = 0; i < (lN*lN); ++i) {
            int r=i/(lN);
            int c=i-r*(lN);
            if(r==0||r==(lN-1)||c==0||c==(lN-1))
            {
                if(r==0&&c!=0&&c!=(lN-1))
                {
                    lunew[i]=((ba[c]+lu[(r)*(lN)+c-1]+lu[(r+1)*(lN)+c]+lu[(r)*(lN)+c+1]-4*lu[i])/hsq)+f[i];
                }
                if(r==(lN-1)&&c!=0&&c!=(lN-1))
                {
                    lunew[i]=((lu[(r-1)*(lN)+c]+lu[(r)*(lN)+c-1]+ta[c]+lu[(r)*(lN)+c+1]-4*lu[i])/hsq)+f[i];
                }
                if(c==0&&r!=0&&r!=(lN-1))
                {
                    lunew[i]=((lu[(r-1)*(lN)+c]+la[r]+lu[(r+1)*(lN)+c]+lu[(r)*(lN)+c+1]-4*lu[i])/hsq)+f[i];
                }
                if(c==(lN-1)&&r!=0&&r!=(lN-1))
                {
                    lunew[i]=((lu[(r-1)*(lN)+c]+lu[(r)*(lN)+c-1]+lu[(r+1)*(lN)+c]+ra[r]-4*lu[i])/hsq)+f[i];
                }
                if(r==0&&c==0)
                    lunew[i]=((ba[c]+la[r]+lu[(r+1)*(lN)+c]+lu[(r)*(lN)+c+1]-4*lu[i])/hsq)+f[i];
                if(r==0&&c==(lN-1))
                    lunew[i]=((ba[c]+lu[(r)*(lN)+c-1]+lu[(r+1)*(lN)+c]+ra[r]-4*lu[i])/hsq)+f[i];
                if(r==(lN-1)&&c==0)
                    lunew[i]=((lu[(r-1)*(lN)+c]+la[r]+ta[c]+lu[(r)*(lN)+c+1]-4*lu[i])/hsq)+f[i];
                if(r==(lN-1)&&c==(lN-1))
                    lunew[i]=((lu[(r-1)*(lN)+c]+lu[(r)*(lN)+c-1]+ta[c]+ra[r]-4*lu[i])/hsq)+f[i];
            }
            else{
                lunew[i]=((lu[(r-1)*(lN)+c]+lu[(r)*(lN)+c-1]+lu[(r+1)*(lN)+c]+lu[(r)*(lN)+c+1]-4*lu[i])/hsq)+f[i];
                //printf("%d\n",lunew[i]);
            }
            
            
        }
        

                         

    free(ta);
    free(ba);
    free(la);
    free(ra);
    return lunew;
}
