#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "iga.h"
void iga_dirichletinit(double *U_, void *app_)
{
    /* ================ input parameters ================ */
    IGA *app=(IGA*)app_;

    Phys *phys           = app->phys;
    int ndim             = phys[0].ndim;
    int ndof             = phys[0].ndof;
    //double *par_mat      = app->par_mat;

    Nurbs *nurbs         = app->nurbs;
    int porder           = nurbs[0].porder;
    //int *nelem           = nurbs[0].nelem;
    int *nbasis          = nurbs[0].nbasis;
    int *nknot           = nurbs[0].nknot;
    double **knotVector  = nurbs[0].knotVector;
    int nboun            = nurbs[0].nboun;
    int *boun            = nurbs[0].boun;

    BC *bc               = app->bc;
    int *bc_type         = bc[0].bc_type;
    //double *par_periodic = bc[0].par_periodic;
    double *par_dirichlet= bc[0].par_dirichlet;
    //double *par_neumann  = bc[0].par_neumann;

    /* ================================================================ BOUNDARY CONDITIONS ================================================================ */
    /* ----------------------------------------------------------------    Dirichlet BCs    ---------------------------------------------------------------- */
    int ia,ibasis;
    int ax_,ay_,az_;
    int nbasis_[ndim],ibasis_z_0,ibasis_z_1;
    int    nknot_[ndim];
    double *knotVector_[ndim];
    int nrow_dirichlet;
    double *dirichlet;
    for (int iboun=0; iboun<nboun; iboun++)
    {
        switch (boun[iboun]/2)
        {
        case 0:
            nknot_[0]=nknot[1]; knotVector_[0]=knotVector[1]; nbasis_[0]=nbasis[1]; ax_=nbasis[2];
            nknot_[1]=nknot[2]; knotVector_[1]=knotVector[2]; nbasis_[1]=nbasis[2]; ay_=1;
            nknot_[2]=nknot[0]; knotVector_[2]=knotVector[0]; nbasis_[2]=nbasis[0]; az_=nbasis[1]*nbasis[2];
            break;
        case 1:
            nknot_[0]=nknot[0]; knotVector_[0]=knotVector[0]; nbasis_[0]=nbasis[0]; ax_=nbasis[1]*nbasis[2];
            nknot_[1]=nknot[2]; knotVector_[1]=knotVector[2]; nbasis_[1]=nbasis[2]; ay_=1; 
            nknot_[2]=nknot[1]; knotVector_[2]=knotVector[1]; nbasis_[2]=nbasis[1]; az_=nbasis[2];
            break;
        case 2:
            nknot_[0]=nknot[0]; knotVector_[0]=knotVector[0]; nbasis_[0]=nbasis[0]; ax_=nbasis[1]*nbasis[2];
            nknot_[1]=nknot[1]; knotVector_[1]=knotVector[1]; nbasis_[1]=nbasis[1]; ay_=nbasis[2];
            nknot_[2]=nknot[2]; knotVector_[2]=knotVector[2]; nbasis_[2]=nbasis[2]; az_=1;
            break;
        default:
            exit(0);
        }
        switch (boun[iboun]%2)
        {
        case 0:
            ibasis_z_0=0;
            ibasis_z_1=1;
            break;
        case 1:
            ibasis_z_0=nbasis_[2]-1;
            ibasis_z_1=nbasis_[2]-2;
            break;
        default: exit(EXIT_FAILURE);
        }
        // standard
        for (int idof=0; idof<ndof; idof++) {
            if (bc_type[(2*ndof)*boun[iboun]+0*ndof+idof]==0) {
                nrow_dirichlet=nbasis_[0]*nbasis_[1];
                dirichlet     = (double*) malloc(nrow_dirichlet*sizeof(double));
                iga_dirichlet(dirichlet,par_dirichlet,boun[iboun],0,idof,knotVector_,nknot_,nbasis_,ndim,porder);
                ia=0;
                for (int ibasis_x_=0; ibasis_x_<nbasis_[0]; ibasis_x_++) {
                for (int ibasis_y_=0; ibasis_y_<nbasis_[1]; ibasis_y_++) {
                    ibasis=ax_*ibasis_x_+ay_*ibasis_y_+az_*ibasis_z_0;
                    U_[ndof*ibasis+idof]=dirichlet[ia++];
                }}
                free(dirichlet);
        }}
        // high-order
        for (int idof=0; idof<ndof; idof++) {
            if (bc_type[(2*ndof)*boun[iboun]+1*ndof+idof]==0) {
                nrow_dirichlet=nbasis_[0]*nbasis_[1];
                dirichlet     = (double*) malloc(nrow_dirichlet*sizeof(double));
                iga_dirichlet(dirichlet,par_dirichlet,boun[iboun],1,idof,knotVector_,nknot_,nbasis_,ndim,porder);
                ia=0;
                for (int ibasis_x_=0; ibasis_x_<nbasis_[0]; ibasis_x_++) {
                for (int ibasis_y_=0; ibasis_y_<nbasis_[1]; ibasis_y_++) {
                    ibasis=ax_*ibasis_x_+ay_*ibasis_y_+az_*ibasis_z_1;
                    U_[ndof*ibasis+idof]=dirichlet[ia++];
                }}
                free(dirichlet);
        }}
    } // for iboun
}

