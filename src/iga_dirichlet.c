#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "iga.h"
#include "mathutil.h"
void iga_dirichlet(double *dirichlet, double par[], int face_id, int bc_order, int idof, double *knotVector_[], int nknot_[], int nbasis_[], int ndim, int porder)
{
    // ---- current parametrization of Dirichlet B.C.:
    // u_x = a0 X + a1 Y + a2 Z              |1    |   |a0 a1 a2|
    // u_y = a3 X + a4 Y + a5 Z    -->   F = |  1  | + |a3 a4 a5|
    // u_z = a6 X + a7 Y + a8 Z              |    1|   |a6 a7 a8|


    double Z_;
    //int ibasis_z_0;
    int ibasis_z_1;
    // ---- set Z_(=const.)
    switch (face_id%2)
    {
        case 0:
            Z_=knotVector_[ndim-1][0];
            //ibasis_z_0=0;
            ibasis_z_1=1;
            break;
        case 1:
            Z_=knotVector_[ndim-1][nknot_[ndim-1]-1]-1.e-15; // need fixing for end values.
            //ibasis_z_0=nbasis_[2]-1;
            ibasis_z_1=nbasis_[2]-2;
            break;
        default: exit(EXIT_FAILURE);
    }
    // ---- 
    switch (bc_order)
    {
        case 0:
        {
            // ---- map par to the face-local coord.
            double par_[ndim];
            switch(face_id/2)
            {
                case 0: par_[0]=par[ndim*idof+1]; par_[1]=par[ndim*idof+2]; par_[2]=par[ndim*idof+0]; break; // YZ
                case 1: par_[0]=par[ndim*idof+0]; par_[1]=par[ndim*idof+2]; par_[2]=par[ndim*idof+1]; break; // XZ
                case 2: par_[0]=par[ndim*idof+0]; par_[1]=par[ndim*idof+1]; par_[2]=par[ndim*idof+2]; break; // XY
                default: exit(0);
            }
            // ---- apply affine deformation on a boundary: par_[0]*X_+par_[1]*Y_+par_[2]*Z_.
            double *dirichlet_1d[ndim-1];
            double vec_p[porder];
            double mat_p[porder*porder];
            double W_;
            for (int bdim=0; bdim<ndim-1; bdim++)
            {
                dirichlet_1d[bdim]=(double*)malloc(nbasis_[bdim]*sizeof(double));
                // 0th - (p-1)th deriv. at knotVector_[bdim][porder].
                W_=knotVector_[bdim][porder];
                for (int i=0; i<porder; i++) { vec_p[i]=0.0; } vec_p[0]=par_[bdim]*W_; vec_p[1]=par_[bdim];
                for (int i=0; i<porder; i++) { for (int j=0; j<porder; j++) { mat_p[porder*i+j]=evalN(knotVector_[bdim],nknot_[bdim],i,j,porder,W_); }}
                matinv_dot_vec(mat_p,vec_p,porder);
                for (int i=0; i<porder; i++) { dirichlet_1d[bdim][i]=vec_p[i]; }
                // 1st deriv.(=par_[bdim]=const.) at knotVector_[bdim][i].
                for (int i=porder; i<nbasis_[bdim]; i++) { dirichlet_1d[bdim][i]=dirichlet_1d[bdim][i-1]+(knotVector_[bdim][i+porder]-knotVector_[bdim][i])/((double)porder)*par_[bdim]; }
            }
            // ----
            for (int i=0; i<nbasis_[0]; i++) { for (int j=0; j<nbasis_[1]; j++) { dirichlet[nbasis_[1]*i+j]=dirichlet_1d[0][i]+dirichlet_1d[1][j]+par_[2]*Z_; }}
            // ----
            for (int bdim=0; bdim<ndim-1; bdim++) { free(dirichlet_1d[bdim]); }
            break;
        }
        case 1:
        {
            // ---- first, set u_{i,n}=0.
            iga_dirichlet(dirichlet,par,face_id,0,idof,knotVector_,nknot_,nbasis_,ndim,porder);
            // ---- then,  set u_{i,n}=some val.
            double m=0.0/evalN(knotVector_[ndim-1],nknot_[ndim-1],1,ibasis_z_1,porder,Z_);// evalN need fixing for end values...
            for (int i=0; i<nbasis_[0]; i++) { for (int j=0; j<nbasis_[1]; j++) { dirichlet[nbasis_[1]*i+j]+=m; }}
            break;
        }
    }
}
