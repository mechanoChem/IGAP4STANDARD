#include "petscsnes.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "snesiga.h"
#include "iocomm.h"

void iga_setmpicomm(int *nselfIdx_,int **selfIdx_,int *nsendIdx_,int **sendIdx_,int **sendPtr_,int *nrecvIdx_,int **recvIdx_,int **recvPtr_,int nsendPart,int *sendBlock,int nrecvPart,int *recvBlock,int ndof,int porder,int nelem[],int nbasis[]);
void iga_setnz(int *d_nz, int *o_nz, int nbasis[], int ndof, int porder, int nsendPart, int *sendBlock, int nrecvPart, int *recvBlock);

void snesiga_init(void* app_, void* phys_, void* nurbs_, void* bc_, void* snesassmbl_, void* iocomm_)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);


    Phys *phys    = (Phys*)phys_;
    int ndof      = phys[0].ndof;
    int torder    = phys[0].torder;

    
    Nurbs *nurbs  = (Nurbs*)nurbs_;
    int porder    = nurbs[0].porder;
    int *nelem    = nurbs[0].nelem;
    int *nbasis   = nurbs[0].nbasis;
    int nsendPart = nurbs[0].nsendPart;
    int *sendPart = nurbs[0].sendPart;
    int nrecvPart = nurbs[0].nrecvPart;
    int *recvPart = nurbs[0].recvPart;

    BC *bc        = (BC*)bc_;

    IGA *app=(IGA*)app_;

    app->phys  = phys;
    app->nurbs = nurbs;
    app->bc    = bc;

    // ----
    int nselfIdx,*selfIdx,nsendIdx,*sendIdx,*sendPtr,nrecvIdx,*recvIdx,*recvPtr;
    iga_setmpicomm(&nselfIdx,&selfIdx,&nsendIdx,&sendIdx,&sendPtr,&nrecvIdx,&recvIdx,&recvPtr,nsendPart,nurbs[0].sendBlock,nrecvPart,nurbs[0].recvBlock,ndof,porder,nelem,nbasis);
    app->nselfIdx = nselfIdx;
    app->selfIdx  = selfIdx;
    app->nsendIdx = nsendIdx;
    app->sendIdx  = sendIdx;
    app->sendPtr  = sendPtr;
    app->nrecvIdx = nrecvIdx;
    app->recvIdx  = recvIdx;
    app->recvPtr  = recvPtr;

    // ---- I/O
    IOcomm *iocomm = (IOcomm*)iocomm_;
    int nDof=ndof*(nbasis[0]*nbasis[1]*nbasis[2]);
    int *nDof_array;          if (rank==0) { nDof_array      =(int*) malloc(nproc*sizeof(int)); }     else { nDof_array      =(int*) malloc(0); }
    int *iDof_displ_array;    if (rank==0) { iDof_displ_array=(int*) malloc((nproc+1)*sizeof(int)); } else { iDof_displ_array=(int*) malloc(0); }
    assert(MPI_Gather(&nDof,1,MPI_INT,nDof_array,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    if (rank==0) { iDof_displ_array[0]=0; for (int iproc=1; iproc<nproc+1; iproc++) { iDof_displ_array[iproc]=iDof_displ_array[iproc-1]+nDof_array[iproc-1]; } }
    int nDof_global;
    if (rank==0) { nDof_global=iDof_displ_array[nproc]; } assert(MPI_Bcast(&nDof_global,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    int *universalIdx_self=(int*)malloc(nDof*sizeof(int));
    int ia=0, ja=0;
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++) {
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++) {
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++) {
        for (int idof=0; idof<ndof; idof++) { universalIdx_self[ia++]=ndof*(nurbs[0].universalIdx_temp[ja])+idof; } ja++;
    }}}
    free(nurbs[0].universalIdx_temp);
    if (rank==0) { iocomm->universalIdx=(int*)malloc(nDof_global*sizeof(int)); } else { iocomm->universalIdx=(int*)malloc(0); }
    assert(MPI_Gatherv(universalIdx_self,nDof,MPI_INT,iocomm->universalIdx,nDof_array,iDof_displ_array,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    free(universalIdx_self);

    iocomm->nDof             = nDof;
    iocomm->nDof_array       = nDof_array;
    iocomm->iDof_displ_array = iDof_displ_array;
    iocomm->nDof_global      = nDof_global;

    app->U_hist=(double*)malloc((torder+1)*nDof*sizeof(double));

    // globalIdx (maps local idx to global idx) (required by SNES)
    int *globalIdx      = (int*) malloc((nselfIdx+nrecvIdx)*sizeof(int));
    MPI_Request request_send[nsendPart],request_recv[nrecvPart];
    MPI_Status status;
    ia=0;
    int ibasis; 
    int iDof_displ;
    assert(MPI_Scatter(iDof_displ_array,1,MPI_INT,&iDof_displ,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    int *globalIdx_self = (int*) malloc(nselfIdx*sizeof(int));
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
        ibasis=(nbasis[1]*nbasis[2])*ibasis_x+nbasis[2]*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { globalIdx_self[ia++]=iDof_displ+ndof*ibasis+idof; }
    }}}
    for (int i=0; i<nselfIdx; i++) { globalIdx[selfIdx[i]]=globalIdx_self[i]; }
    free(globalIdx_self);
    int *globalIdx_send=(int*) malloc(nsendIdx*sizeof(int));
    int *globalIdx_recv=(int*) malloc(nrecvIdx*sizeof(int));
    for (int i=0; i<nsendIdx; i++) { globalIdx_send[i]=globalIdx[sendIdx[i]]; }
    for (int isendPart=0; isendPart<nsendPart; isendPart++) { assert(MPI_Isend(globalIdx_send+sendPtr[isendPart],sendPtr[isendPart+1]-sendPtr[isendPart],MPI_INT,sendPart[isendPart],1,MPI_COMM_WORLD,request_send+isendPart)==MPI_SUCCESS); }
    for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) { assert(MPI_Irecv(globalIdx_recv+recvPtr[irecvPart],recvPtr[irecvPart+1]-recvPtr[irecvPart],MPI_INT,recvPart[irecvPart],1,MPI_COMM_WORLD,request_recv+irecvPart)==MPI_SUCCESS); }
    for (int isendPart=0; isendPart<nsendPart; isendPart++) assert(MPI_Wait(request_send+isendPart,&status)==MPI_SUCCESS);
    for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) assert(MPI_Wait(request_recv+irecvPart,&status)==MPI_SUCCESS);
    for (int i=0; i<nrecvIdx; i++) { globalIdx[recvIdx[i]]=globalIdx_recv[i]; }
    free(globalIdx_send);
    free(globalIdx_recv);
    app->globalIdx = globalIdx;

    // ---- To work with PETSc SNES : call in main.c snes_init(&app) and snes_solve(&app).
    SNESAssmbl *snesassmbl = (SNESAssmbl*)snesassmbl_;
    // Vec
    VecCreateMPI(MPI_COMM_WORLD,nDof,nDof_global,&snesassmbl->Residual_);
    VecDuplicate(snesassmbl->Residual_,&snesassmbl->U);
    //  Mat
    int *d_nz=(int*) malloc(nDof*sizeof(int)),*o_nz=(int*) malloc(nDof*sizeof(int));
    iga_setnz(d_nz,o_nz,nbasis,ndof,porder,nsendPart,nurbs[0].sendBlock,nrecvPart,nurbs[0].recvBlock);
    MatCreateAIJ(MPI_COMM_WORLD,nDof,nDof,nDof_global,nDof_global,0,d_nz,0,o_nz,&snesassmbl->Tangent_);
    free(nurbs[0].sendBlock);
    free(nurbs[0].recvBlock);
    free(d_nz);
    free(o_nz);
    MatSetOption(snesassmbl->Tangent_,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
    MatSetOption(snesassmbl->Tangent_,MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE);
    MatSetOption(snesassmbl->Tangent_,MAT_SUBSET_OFF_PROC_ENTRIES,PETSC_TRUE);
    //MatSetOption(Tangent_,MAT_SYMMETRIC,PETSC_TRUE);
    //MatSetOption(Tangent_,MAT_SYMMETRY_ETERNAL,PETSC_TRUE);


    // Function pointers to assembly routines.
    snesassmbl->snes_Residual=snesiga_Residual;
    snesassmbl->snes_Tangent =snesiga_Tangent;
}

void iga_setmpicomm(
int *nselfIdx_,
int **selfIdx_,
int *nsendIdx_,
int **sendIdx_,
int **sendPtr_,
int *nrecvIdx_,
int **recvIdx_,
int **recvPtr_,
int nsendPart,
int *sendBlock,
int nrecvPart,
int *recvBlock,
int ndof,
int porder,
int nelem[],
int nbasis[]
)
{
    int nbasis_y_active=nelem[1]+porder;
    int nbasis_z_active=nelem[2]+porder;
    
    // ---- 

    int ibasis,ia;
    // ---- self
    int nselfIdx=ndof*(nbasis[0]*nbasis[1]*nbasis[2]);
    int *selfIdx=(int*)malloc(nselfIdx*sizeof(int));
    ia=0;
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { selfIdx[ia++]=ndof*ibasis+idof; }
    }}}

    // ---- send
    int nsendIdx=0;
    int *sendPtr=(int*)malloc((nsendPart+1)*sizeof(int)); sendPtr[0]=nsendIdx;
    for (int i=0; i<nsendPart; i++)
    {
        switch (sendBlock[i]) {
        case 0: nsendIdx+=ndof*(nbasis[0]*nbasis[1]*porder); break;
        case 1: nsendIdx+=ndof*(nbasis[0]*nbasis[2]*porder); break;
        case 2: nsendIdx+=ndof*(nbasis[1]*nbasis[2]*porder); break;
        case 3: nsendIdx+=ndof*(nbasis[0]*porder*porder); break;
        case 4: nsendIdx+=ndof*(nbasis[1]*porder*porder); break;
        case 5: nsendIdx+=ndof*(nbasis[2]*porder*porder); break;
        case 6: nsendIdx+=ndof*(porder*porder*porder); break;
        }
        sendPtr[i+1]=nsendIdx;
    }
    int *sendIdx=(int*)malloc(nsendIdx*sizeof(int));
    ia=0;
    for (int i=0; i<nsendPart; i++)
    {
        switch (sendBlock[i]) {
        case 0: // (face0)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 1: // (face1)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 2: // (face2)
            for (int ibasis_x=0; ibasis_x<porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 3: // (edge0)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 4: // (edge1)
            for (int ibasis_x=0; ibasis_x<porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 5: // (edge2)
            for (int ibasis_x=0; ibasis_x<porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 6: // (corner)
            for (int ibasis_x=0; ibasis_x<porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        }
    }

    // ---- recv
    int nrecvIdx=0;
    int *recvPtr=(int*)malloc((nrecvPart+1)*sizeof(int)); recvPtr[0]=nrecvIdx;
    for (int i=0; i<nrecvPart; i++)
    {
        switch (recvBlock[i]) {
        case 0: nrecvIdx+=ndof*(nbasis[0]*nbasis[1]*porder); break;
        case 1: nrecvIdx+=ndof*(nbasis[0]*nbasis[2]*porder); break;
        case 2: nrecvIdx+=ndof*(nbasis[1]*nbasis[2]*porder); break;
        case 3: nrecvIdx+=ndof*(nbasis[0]*porder*porder); break;
        case 4: nrecvIdx+=ndof*(nbasis[1]*porder*porder); break;
        case 5: nrecvIdx+=ndof*(nbasis[2]*porder*porder); break;
        case 6: nrecvIdx+=ndof*(porder*porder*porder); break;
        }
        recvPtr[i+1]=nrecvIdx;
    }
    int *recvIdx=(int*)malloc(nrecvIdx*sizeof(int));
    ia=0;
    for (int i=0; i<nrecvPart; i++)
    {
        switch (recvBlock[i]) {
        case 0: // (face0)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=nbasis[2]; ibasis_z<nbasis[2]+porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 1: // (face1)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=nbasis[1]; ibasis_y<nbasis[1]+porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 2: // (face2)
            for (int ibasis_x=nbasis[0]; ibasis_x<nbasis[0]+porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 3: // (edge0)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=nbasis[1]; ibasis_y<nbasis[1]+porder  ; ibasis_y++){
            for (int ibasis_z=nbasis[2]; ibasis_z<nbasis[2]+porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 4: // (edge1)
            for (int ibasis_x=nbasis[0]; ibasis_x<nbasis[0]+porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=nbasis[2]; ibasis_z<nbasis[2]+porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 5: // (edge2)
            for (int ibasis_x=nbasis[0]; ibasis_x<nbasis[0]+porder  ; ibasis_x++){
            for (int ibasis_y=nbasis[1]; ibasis_y<nbasis[1]+porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 6: // (corner)
            for (int ibasis_x=nbasis[0]; ibasis_x<nbasis[0]+porder  ; ibasis_x++){
            for (int ibasis_y=nbasis[1]; ibasis_y<nbasis[1]+porder  ; ibasis_y++){
            for (int ibasis_z=nbasis[2]; ibasis_z<nbasis[2]+porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        }
    }
    // ----
    *nselfIdx_      = nselfIdx;
    *nsendIdx_      = nsendIdx;
    *nrecvIdx_      = nrecvIdx;

    *selfIdx_       = selfIdx;
    *sendIdx_       = sendIdx;
    *sendPtr_       = sendPtr;
    *recvIdx_       = recvIdx;
    *recvPtr_       = recvPtr;
}

void iga_setnz(int *d_nz, int *o_nz, int nbasis[], int ndof, int porder, int nsendPart, int *sendBlock, int nrecvPart, int *recvBlock)
{
    int sd[]={0,0,0}; for (int i=0; i<nsendPart; i++) { if (sendBlock[i]<3) sd[sendBlock[i]]=1; }
    int rv[]={0,0,0}; for (int i=0; i<nrecvPart; i++) { if (recvBlock[i]<3) rv[recvBlock[i]]=1; }
    int ia=0;
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++) {
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++) {
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++) {
    for (int idof=0; idof<ndof; idof++) {
        d_nz[ia]=ndof*( ( ibasis_x < porder ? ibasis_x : porder )+1+( ibasis_x > (nbasis[0]-1)-porder ? (nbasis[0]-1)-ibasis_x:porder) )
                     *( ( ibasis_y < porder ? ibasis_y : porder )+1+( ibasis_y > (nbasis[1]-1)-porder ? (nbasis[1]-1)-ibasis_y:porder) )
                     *( ( ibasis_z < porder ? ibasis_z : porder )+1+( ibasis_z > (nbasis[2]-1)-porder ? (nbasis[2]-1)-ibasis_z:porder) );
        o_nz[ia]=ndof*( ( (sd[2]==0 && ibasis_x < porder) ? ibasis_x : porder )+1+( (rv[2]==0 && ibasis_x > (nbasis[0]-1)-porder ) ? (nbasis[0]-1)-ibasis_x : porder ) )
                     *( ( (sd[1]==0 && ibasis_y < porder) ? ibasis_y : porder )+1+( (rv[1]==0 && ibasis_y > (nbasis[1]-1)-porder ) ? (nbasis[1]-1)-ibasis_y : porder ) )
                     *( ( (sd[0]==0 && ibasis_z < porder) ? ibasis_z : porder )+1+( (rv[0]==0 && ibasis_z > (nbasis[2]-1)-porder ) ? (nbasis[2]-1)-ibasis_z : porder ) )
                 -d_nz[ia];
        ia++;
    }}}}
}
// if send to/recv from self: d_nz=d_nz+o_nz, o_nz=0.
