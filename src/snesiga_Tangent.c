#include "petscsnes.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "snesiga.h"
#include "mathutil.h"
#include "quadrature.h"

#undef __FUNCT__
#define __FUNCT__ "snesiga_Tangent"

PetscErrorCode snesiga_Tangent(SNES snes, Vec U_, Mat Tangent_, Mat Pmat_, void *app_)
{
    /* ================ input parameters ================ */
    IGA *app=(IGA*)app_;
    
    Phys *phys           = app->phys;
    int ndim             = phys[0].ndim;
    int ndof             = phys[0].ndof;
    int torder           = phys[0].torder;
    double *par_mat      = phys[0].par_mat;
    tangent_ptr tangent  = phys[0].tangent;
    
    Nurbs *nurbs         = app->nurbs;
    int porder           = nurbs[0].porder;
    int *nelem           = nurbs[0].nelem;
    int *nbasis          = nurbs[0].nbasis;
    //int *nknot           = nurbs[0].nknot;
    double **knotVector  = nurbs[0].knotVector;
    int nboun            = nurbs[0].nboun;
    int *boun            = nurbs[0].boun;
    int nsendPart        = nurbs[0].nsendPart;
    int *sendPart        = nurbs[0].sendPart;
    int nrecvPart        = nurbs[0].nrecvPart;
    int *recvPart        = nurbs[0].recvPart;

    BC *bc               = app->bc;
    int *bc_type         = bc[0].bc_type;
    double *par_periodic = bc[0].par_periodic;
    //double *par_dirichlet= bc[0].par_dirichlet;
    //double *par_neumann  = bc[0].par_neumann;

    int nselfIdx         = app->nselfIdx;
    int *selfIdx         = app->selfIdx;
    int nsendIdx         = app->nsendIdx;
    int *sendIdx         = app->sendIdx;
    int *sendPtr         = app->sendPtr;
    int nrecvIdx         = app->nrecvIdx;
    int *recvIdx         = app->recvIdx;
    int *recvPtr         = app->recvPtr;
    int *globalIdx       = app->globalIdx;
    /* ================ assemble U ================ */
    MPI_Request request_send[nsendPart],request_recv[nrecvPart];
    MPI_Status status;
    double *sendbuff=(double*) malloc(nsendIdx*sizeof(double));
    double *recvbuff=(double*) malloc(nrecvIdx*sizeof(double));
    double *U = (double*) malloc((torder+1)*(nselfIdx+nrecvIdx)*sizeof(double));
    double *U_hist = app->U_hist;
    const double *U_curr; VecGetArrayRead(U_,&U_curr); for (int i=0; i<nselfIdx; i++) { U_hist[i]=U_curr[i]; } VecRestoreArrayRead(U_,&U_curr);
    for (int ihist=0; ihist<torder+1; ihist++)
    {
        for (int i=0; i<nselfIdx; i++) { U[ihist*(nselfIdx+nrecvIdx)+selfIdx[i]]=U_hist[ihist*nselfIdx+i]; }
        for (int i=0; i<nsendIdx; i++) { sendbuff[i]=U[ihist*(nselfIdx+nrecvIdx)+sendIdx[i]]; }
        for (int isendPart=0; isendPart<nsendPart; isendPart++) 
        { assert(MPI_Isend(sendbuff+sendPtr[isendPart],sendPtr[isendPart+1]-sendPtr[isendPart],MPI_DOUBLE,sendPart[isendPart],30+ihist,MPI_COMM_WORLD,request_send+isendPart)==MPI_SUCCESS);}
        for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) 
        { assert(MPI_Irecv(recvbuff+recvPtr[irecvPart],recvPtr[irecvPart+1]-recvPtr[irecvPart],MPI_DOUBLE,recvPart[irecvPart],30+ihist,MPI_COMM_WORLD,request_recv+irecvPart)==MPI_SUCCESS);}
        for (int isendPart=0; isendPart<nsendPart; isendPart++) assert(MPI_Wait(request_send+isendPart,&status)==MPI_SUCCESS);
        for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) assert(MPI_Wait(request_recv+irecvPart,&status)==MPI_SUCCESS);
        for (int i=0; i<nrecvIdx; i++) { U[ihist*(nselfIdx+nrecvIdx)+recvIdx[i]]=recvbuff[i]; }
    }
    free(sendbuff);
    free(recvbuff);
    /* ================ assemble Tangent_ ================ */
    // ---- Mesh
    int ibpe,lbpe;
    int nbpe = int_pow(porder+1,3);
    int ibasis;
    double X0[ndim],X1[ndim];
    double *kV[ndim];
    // ---- PDE
    int nddim=1+ndim+ndim*(ndim+1)/2;
    double *tangent_ = (double*) malloc(int_pow(ndof*nddim,2)*sizeof(double));
    double u[(torder+1)*nddim*ndof];
    // ---- quadrature
    double xq[ndim],xi;
    double weight;
    // ---- IGA
    double N[nbpe*nddim];
    MatZeroEntries(Pmat_);
    double val[int_pow(ndof*nbpe,2)];
    int         idx[ndof*nbpe];
    double      temp;
    int         ia,row,col;
    int nbasis_y_active=nelem[1]+porder;
    int nbasis_z_active=nelem[2]+porder;
    /* ================================================================   VOLUME INTEGRAL   ================================================================ */
    int    nquad=4;
    double cquad[nquad];
    double wquad[nquad];
    legendre_handle(cquad,wquad,nquad,0.,1.);
    // ---- loop over elements
    for (int ielem_x=0; ielem_x<nelem[0]; ielem_x++) {
    for (int ielem_y=0; ielem_y<nelem[1]; ielem_y++) {
    for (int ielem_z=0; ielem_z<nelem[2]; ielem_z++) {
    for (int i=0; i<int_pow(ndof*nbpe,2); i++) { val[i]=0.0; }
    kV[0]=knotVector[0]+ielem_x;
    kV[1]=knotVector[1]+ielem_y;
    kV[2]=knotVector[2]+ielem_z;
    // ---- loop over quadrature points
    for (int iquad_x=0; iquad_x<nquad; iquad_x++) {
    for (int iquad_y=0; iquad_y<nquad; iquad_y++) {
    for (int iquad_z=0; iquad_z<nquad; iquad_z++) {
        // ---- evaluate quadrature coord.
        X0[0]=knotVector[0][ielem_x+porder]; X1[0]=knotVector[0][ielem_x+porder+1];
        X0[1]=knotVector[1][ielem_y+porder]; X1[1]=knotVector[1][ielem_y+porder+1];
        X0[2]=knotVector[2][ielem_z+porder]; X1[2]=knotVector[2][ielem_z+porder+1];
        xi=cquad[iquad_x];xq[0]=(-xi+1.0)*X0[0]+xi*X1[0];
        xi=cquad[iquad_y];xq[1]=(-xi+1.0)*X0[1]+xi*X1[1];
        xi=cquad[iquad_z];xq[2]=(-xi+1.0)*X0[2]+xi*X1[2];
        // ---- evaluate N
        eval_Bspline(N,kV,porder,xq);
        // ---- evaluate u
        for (ia=0; ia<(torder+1)*nddim*ndof; ia++) { u[ia]=0.0; }
        for (int ibpe_x=0; ibpe_x<porder+1; ibpe_x++){
        for (int ibpe_y=0; ibpe_y<porder+1; ibpe_y++){
        for (int ibpe_z=0; ibpe_z<porder+1; ibpe_z++){
            ibpe=(porder+1)*(porder+1)*ibpe_x+(porder+1)*ibpe_y+ibpe_z;
            ibasis=(nbasis_y_active*nbasis_z_active)*(ielem_x+ibpe_x)+nbasis_z_active*(ielem_y+ibpe_y)+(ielem_z+ibpe_z);
            for (int ihist=0; ihist<torder+1; ihist++)
            {
                for (int idof=0; idof<ndof; idof++)
                {
                    ia=ndof*ibasis+idof;
                    for (int iddim=0; iddim<nddim; iddim++) { u[(nddim*ndof)*ihist+nddim*idof+iddim]+=N[nddim*ibpe+iddim]*U[(nselfIdx+nrecvIdx)*ihist+ia]; }
                }
            }
        }}}
        for (int ihist=0; ihist<torder+1; ihist++)
        {
            for (int idof=0; idof<ndof; idof++) { u[(nddim*ndof)*ihist+nddim*idof+0]+=par_periodic[ndim*idof+0]*xq[0]+par_periodic[ndim*idof+1]*xq[1]+par_periodic[ndim*idof+2]*xq[2]; }
            for (int idof=0; idof<ndof; idof++) { for (int idim=0; idim<ndim; idim++) { u[(nddim*ndof)*ihist+nddim*idof+(1+idim)]+=par_periodic[ndim*idof+idim]; } }
        }
        // ---- evaluate tangent_
        tangent(tangent_,u,par_mat);
        // ---- add to val
        weight=wquad[iquad_x]*wquad[iquad_y]*wquad[iquad_z]*(X1[0]-X0[0])*(X1[1]-X0[1])*(X1[2]-X0[2]);
        for (int ibpe_x=0; ibpe_x<porder+1; ibpe_x++){
        for (int ibpe_y=0; ibpe_y<porder+1; ibpe_y++){
        for (int ibpe_z=0; ibpe_z<porder+1; ibpe_z++){
            ibpe=(porder+1)*(porder+1)*ibpe_x+(porder+1)*ibpe_y+ibpe_z;
            for (int idof=0; idof<ndof; idof++)
            {
                row=ndof*ibpe+idof;
                for (int lbpe_x=0; lbpe_x<porder+1; lbpe_x++){
                for (int lbpe_y=0; lbpe_y<porder+1; lbpe_y++){
                for (int lbpe_z=0; lbpe_z<porder+1; lbpe_z++){
                    lbpe=(porder+1)*(porder+1)*lbpe_x+(porder+1)*lbpe_y+lbpe_z;
                    for (int ldof=0; ldof<ndof; ldof++)
                    {
                        col=ndof*lbpe+ldof;
                        temp=0.0;
                        for (int iddim=0; iddim<nddim; iddim++)
                        {
                            for (int lddim=0; lddim<nddim; lddim++)
                            {
                                temp+=N[nddim*ibpe+iddim]*tangent_[(ndof*nddim)*(nddim*idof+iddim)+(nddim*ldof+lddim)]*N[nddim*lbpe+lddim];
                            } // lddim
                        } // iddim
                        val[(ndof*nbpe)*row+col]+=weight*temp;
                    } // ldof
                }}} // lbpe
            } // idof
        }}} // ibpe
    }}} // iquad
    // ---- set idx
    for (int ibpe_x=0; ibpe_x<porder+1; ibpe_x++){
    for (int ibpe_y=0; ibpe_y<porder+1; ibpe_y++){
    for (int ibpe_z=0; ibpe_z<porder+1; ibpe_z++){
        ibpe=(porder+1)*(porder+1)*ibpe_x+(porder+1)*ibpe_y+ibpe_z;
        ibasis=(nbasis_y_active*nbasis_z_active)*(ielem_x+ibpe_x)+nbasis_z_active*(ielem_y+ibpe_y)+(ielem_z+ibpe_z);
        for (int idof=0; idof<ndof; idof++)
        {
            idx[ndof*ibpe+idof]=globalIdx[ndof*ibasis+idof];
        }
    }}}
    // ---- add to Pmat_
    MatSetValues(Pmat_,ndof*nbpe,idx,ndof*nbpe,idx,val,ADD_VALUES);
    }}} // ielem
    // ---- 
    MatAssemblyBegin(Pmat_,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Pmat_,MAT_FINAL_ASSEMBLY);

/*Mat Atrans;
MatTranspose(Pmat_,MAT_INITIAL_MATRIX,&Atrans);
MatAXPY(Atrans,-1.0,Pmat_,SAME_NONZERO_PATTERN);
double nm;
MatNorm(Atrans,NORM_INFINITY,&nm);
printf("%e\n",nm);
*/
    /* ================================================================ BOUNDARY CONDITIONS ================================================================ */
    int nbasis_x_,nbasis_y_,nbasis_z_,ibasis_z_0,ibasis_z_1;
    int ax_,ay_,az_;
    /* ----------------------------------------------------------------    Dirichlet BCs    ---------------------------------------------------------------- */
    int nrow_dirichlet=0;
    int itemp;
    for (int iboun=0; iboun<nboun; iboun++)
    {
        if      (boun[iboun]/2==0) { itemp=(2*ndof); for (int i=0; i<(2*ndof); i++) { itemp-=bc_type[boun[iboun]*(2*ndof)+i]; } nrow_dirichlet += itemp*nbasis[1]*nbasis[2]; }
        else if (boun[iboun]/2==1) { itemp=(2*ndof); for (int i=0; i<(2*ndof); i++) { itemp-=bc_type[boun[iboun]*(2*ndof)+i]; } nrow_dirichlet += itemp*nbasis[2]*nbasis[0]; }
        else                       { itemp=(2*ndof); for (int i=0; i<(2*ndof); i++) { itemp-=bc_type[boun[iboun]*(2*ndof)+i]; } nrow_dirichlet += itemp*nbasis[0]*nbasis[1]; }
    }
    int *row_dirichlet = (int*) malloc(nrow_dirichlet*sizeof(int));
    ia=0;
    for (int iboun=0; iboun<nboun; iboun++)
    {
        switch (boun[iboun]/2)
        {
        case 0:
            nbasis_x_=nbasis[1]; ax_=nbasis_z_active;
            nbasis_y_=nbasis[2]; ay_=1;
            nbasis_z_=nbasis[0]; az_=nbasis_y_active*nbasis_z_active;
            break;
        case 1:
            nbasis_x_=nbasis[0]; ax_=nbasis_y_active*nbasis_z_active;
            nbasis_y_=nbasis[2]; ay_=1;
            nbasis_z_=nbasis[1]; az_=nbasis_z_active;
            break;
        case 2:
            nbasis_x_=nbasis[0]; ax_=nbasis_y_active*nbasis_z_active;
            nbasis_y_=nbasis[1]; ay_=nbasis_z_active;
            nbasis_z_=nbasis[2]; az_=1;
            break;
        default:
            printf("error: \n");
            exit(0);
        }
        switch (boun[iboun]%2)
        {
        case 0:
            ibasis_z_0=0;
            ibasis_z_1=1;
            break;
        case 1:
            ibasis_z_0=nbasis_z_-1;
            ibasis_z_1=nbasis_z_-2;
            break;
        default: exit(EXIT_FAILURE);
        }
        // standard
        for (int idof=0; idof<ndof; idof++) {
            if (bc_type[boun[iboun]*(2*ndof)+0*ndof+idof]==0) {
                for (int ibasis_x_=0; ibasis_x_<nbasis_x_; ibasis_x_++) {
                for (int ibasis_y_=0; ibasis_y_<nbasis_y_; ibasis_y_++) {
                    ibasis=ax_*ibasis_x_+ay_*ibasis_y_+az_*ibasis_z_0;
                    row_dirichlet[ia++]=globalIdx[ndof*ibasis+idof];
        }}}}
        // high-order
        for (int idof=0; idof<ndof; idof++) {
            if (bc_type[boun[iboun]*(2*ndof)+1*ndof+idof]==0) {
                for (int ibasis_x_=0; ibasis_x_<nbasis_x_; ibasis_x_++) {
                for (int ibasis_y_=0; ibasis_y_<nbasis_y_; ibasis_y_++) {
                    ibasis=ax_*ibasis_x_+ay_*ibasis_y_+az_*ibasis_z_1;
                    row_dirichlet[ia++]=globalIdx[ndof*ibasis+idof];
        }}}}
    } // for iboun
    assert(ia==nrow_dirichlet);

    MatZeroRowsColumns(Pmat_,nrow_dirichlet,row_dirichlet,1.0,PETSC_NULL,PETSC_NULL); // symmetrized
    //MatZeroRows(Pmat_,nrow_dirichlet,row_dirichlet,1.0,PETSC_NULL,PETSC_NULL); // unsymmetrized
    free(row_dirichlet);
    
    /* ---------------------------------------------------------------- FINAL ASSEMBLY ---------------------------------------------------------------- */
    if (Tangent_ != Pmat_)
    {
        MatAssemblyBegin(Tangent_,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Tangent_,MAT_FINAL_ASSEMBLY);
    }

    free(U);
    free(tangent_);
    return 0;
}

