#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "nurbs.h"
#include "mathutil.h"
void nurbs_genuniform(void *nurbs_, int mref, int porder, int bc_periodic_x, int bc_periodic_y, int bc_periodic_z)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    int ndim=3;

    int *bc_periodic=(int*)malloc(ndim*sizeof(int));
    bc_periodic[0]=bc_periodic_x;
    bc_periodic[1]=bc_periodic_y;
    bc_periodic[2]=bc_periodic_z;

    int *nelem_global=(int*)malloc(ndim*sizeof(int));
    double **knotVector_global=(double **)malloc(ndim*sizeof(double*));

    for (int idim=0; idim<ndim; idim++) { nelem_global[idim]=int_pow(2,mref); }
    if (rank==0) { for (int idim=0; idim<ndim; idim++) { knotVector_global[idim]=(double*) malloc((nelem_global[idim]+1+2*porder)*sizeof(double)); } }
    else         { for (int idim=0; idim<ndim; idim++) { knotVector_global[idim]=(double*) malloc(0*sizeof(double)); } }

    if (rank==0)
    {
        double L[]={1.0,1.0,1.0};
        double l[ndim]; for (int idim=0; idim<ndim; idim++) { l[idim]=L[idim]/nelem_global[idim]; }
        for (int idim=0; idim<ndim; idim++)
        {
            int iknot=0;
            if (bc_periodic[idim]==0)
            {
                while (iknot < porder+1)                      { knotVector_global[idim][iknot++] = 0.; }
                while (iknot < nelem_global[idim]+porder)     { knotVector_global[idim][iknot]   = (iknot-porder)*l[idim]+0*sin(iknot)*0.25*l[idim]; iknot++; }
                while (iknot < nelem_global[idim]+1+2*porder) { knotVector_global[idim][iknot++] = L[idim]; }
            }
            else
            {
                while (iknot < nelem_global[idim]+1+2*porder) { knotVector_global[idim][iknot] = (iknot-porder)*l[idim]; iknot++; }
            }
        }
    }

    Nurbs *nurbs=(Nurbs*)nurbs_;
    nurbs->porder            = porder;
    nurbs->bc_periodic       = bc_periodic;
    nurbs->nelem_global      = nelem_global;
    nurbs->knotVector_global = knotVector_global;
}
