#include "petscsnes.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "iga.h"
void free_IGA(void *iga_)
{
    IGA *iga = (IGA*)iga_;
   
    Nurbs *nurbs = iga->nurbs;
    free(nurbs[0].nelem);
    free(nurbs[0].nbasis);
int ndim=3;
    for (int idim=0; idim<ndim; idim++) { free(nurbs[0].knotVector[idim]); } // allocated in "set_mpi_part.c"
    free(nurbs[0].knotVector);
    free(nurbs[0].sendPart);                                                 //
    free(nurbs[0].recvPart);                                                 //

    BC *bc = iga->bc;
    free(bc[0].bc_type);
    free(bc[0].par_periodic);
    free(bc[0].par_dirichlet);
    free(bc[0].par_neumann);


    free(iga->selfIdx);                                                  // allocated in "set_mpi_comm.c"
    free(iga->sendIdx);                                                  //
    free(iga->sendPtr);                                                  //
    free(iga->recvIdx);                                                  //
    free(iga->recvPtr);                                                  //
    free(iga->globalIdx);

    free(iga->U_hist);
    
}
