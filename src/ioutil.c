#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <assert.h>
#include "ioutil.h"
#include "iocomm.h"


void mpi_printf( const char * format, ... )
{
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    if (rank==0)
    {
        va_list args;
        va_start (args, format);
        vprintf (format, args);
        va_end (args);
    }
}

int file_exist(char fname[])
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    int fexist;
    FILE *fptr;
    if (rank == 0)
    {
        fexist = ( ((fptr=fopen(fname,"rb")) != NULL) ? 1:0 );
        if (fexist == 1) { fclose(fptr); }
    }
    assert(MPI_Bcast(&fexist,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);

    return fexist;
}

void write_iarray(const char *fname, const char *mode, int iarray[], int size)
{
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    if (rank==0)
    {
        FILE *fptr=fopen(fname,mode);
        fwrite(iarray,sizeof(int),size,fptr);
        fclose(fptr);
    }
}

void read_iarray(const char *fname, const char *mode, int iarray[], int size)
{
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);
    if (rank==0)
    {
        FILE *fptr=fopen(fname,mode);
        assert(fread(iarray,sizeof(int),size,fptr)==size);
        fclose(fptr);
    }
    assert(MPI_Bcast(iarray,size,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
}

void write_darray(const char *fname, const char *mode, double darray[], int size)
{
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    if (rank==0)
    {
        FILE *fptr=fopen(fname,mode);
        fwrite(darray,sizeof(double),size,fptr);
        fclose(fptr);
    }
}

void read_darray(const char *fname, const char *mode, double darray[], int size)
{
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);
    if (rank==0)
    { 
        FILE *fptr=fopen(fname,mode);
        assert(fread(darray,sizeof(double),size,fptr)==size);
        fclose(fptr);
    }
    assert(MPI_Bcast(darray,size,MPI_DOUBLE,0,MPI_COMM_WORLD)==MPI_SUCCESS);
}

void write_solution(const char *fname, double *U_hist, void *iocomm_, int ihist0, int ihist1)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    IOcomm *iocomm        = (IOcomm*)iocomm_;
    int nDof              = iocomm->nDof;
    int *nDof_array       = iocomm->nDof_array;
    int *iDof_displ_array = iocomm->iDof_displ_array;
    int nDof_global       = iocomm->nDof_global;
    int *universalIdx     = iocomm->universalIdx;

    double *U_universal     ; if (rank==0) { U_universal     =(double*) malloc((ihist1-ihist0)*nDof_global*sizeof(double)); } else { U_universal     =(double*) malloc(0); }
    double *U_universal_temp; if (rank==0) { U_universal_temp=(double*) malloc(                nDof_global*sizeof(double)); } else { U_universal_temp=(double*) malloc(0); }

    for (int ihist=ihist0; ihist<ihist1; ihist++)
    {
        assert(MPI_Gatherv(U_hist+nDof*ihist,nDof,MPI_DOUBLE,U_universal_temp,nDof_array,iDof_displ_array,MPI_DOUBLE,0,MPI_COMM_WORLD)==MPI_SUCCESS);
        if (rank==0)
        {
            for (int i=0; i<nDof_global; i++)
            {
                U_universal[(ihist-ihist0)*nDof_global+universalIdx[i]]=U_universal_temp[i];
            }
        }
    }
    if (rank==0)
    {
        FILE *fptr=fopen(fname,"wb");
        fwrite(U_universal,sizeof(double),(ihist1-ihist0)*nDof_global,fptr);
        fclose(fptr);
    }

    free(U_universal);
    free(U_universal_temp);
}

void read_solution(const char *fname, double *U_hist, void *iocomm_, int ihist0, int ihist1)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    IOcomm *iocomm        = (IOcomm*)iocomm_;
    int nDof              = iocomm->nDof;
    int *nDof_array       = iocomm->nDof_array;
    int *iDof_displ_array = iocomm->iDof_displ_array;
    int nDof_global       = iocomm->nDof_global;
    int *universalIdx     = iocomm->universalIdx;

    double *U_universal     ; if (rank==0) { U_universal     =(double*) malloc((ihist1-ihist0)*nDof_global*sizeof(double)); } else { U_universal     =(double*) malloc(0); }
    double *U_universal_temp; if (rank==0) { U_universal_temp=(double*) malloc(                nDof_global*sizeof(double)); } else { U_universal_temp=(double*) malloc(0); }

    if (rank==0)
    { 
        FILE *fptr=fopen(fname,"rb");
        assert(fread(U_universal,sizeof(double),(ihist1-ihist0)*nDof_global,fptr)==(ihist1-ihist0)*nDof_global);
        fclose(fptr);
    }
    for (int ihist=ihist0; ihist<ihist1; ihist++)
    {
        if (rank==0)
        { 
            for (int i=0; i<nDof_global; i++) { U_universal_temp[i]=U_universal[(ihist-ihist0)*nDof_global+universalIdx[i]]; }
        }
        assert(MPI_Scatterv(U_universal_temp,nDof_array,iDof_displ_array,MPI_DOUBLE,U_hist+nDof*ihist,nDof,MPI_DOUBLE,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    }

    free(U_universal);
    free(U_universal_temp);
}












