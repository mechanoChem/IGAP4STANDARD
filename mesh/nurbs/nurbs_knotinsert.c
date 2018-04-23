#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "nurbs.h"
#include "mathutil.h"
void nurbs_knotinsert(double *U_universal, const int nelem[], int nref, int porder, int ndim, int ndof, int bc_periodic[])
{
    assert(porder==2); // currently only supports porder==2
    double vx[4],vy[4],vz[4];
    int nelem_old[ndim];
    int nbasis_old[ndim];
    int nbasis_new[ndim];
    int nbasis[ndim];
    double *Ui_new,*Ui_old;

    Ui_new=(double*) malloc((nelem[0]+porder)*(nelem[1]+porder)*(nelem[2]+porder)*sizeof(double));
    Ui_old=(double*) malloc((nelem[0]+porder)*(nelem[1]+porder)*(nelem[2]+porder)*sizeof(double));
    for (int idof=0; idof<ndof; idof++)
    {
        for (int idim=0; idim<ndim; idim++) {
            nelem_old[idim] =nelem[idim]/int_pow(2,nref);
            nbasis_old[idim]=nelem_old[idim]+porder;
            nbasis[idim]    =nelem_old[idim]+(bc_periodic[idim]==0)*porder;
        }
        for (int i=0; i<nbasis_old[0]; i++) { 
        for (int j=0; j<nbasis_old[1]; j++) { 
        for (int k=0; k<nbasis_old[2]; k++) { 
            Ui_old[nbasis_old[1]*nbasis_old[2]*i+nbasis_old[2]*j+k]=U_universal[ndof*(nbasis[1]*nbasis[2]*(i%nbasis[0])+nbasis[2]*(j%nbasis[1])+(k%nbasis[2]))+idof];
        }}}
        for (int iref=0; iref<nref; iref++)
        {
            for (int idim=0; idim<ndim; idim++) { nbasis_new[idim]=2*nelem_old[idim]+porder; }
            for (int i=0; i<nbasis_new[0]*nbasis_new[1]*nbasis_new[2]; i++) { Ui_new[i]=0; }

            for (int ix=0; ix<nbasis_old[0]; ix++) {
            for (int iy=0; iy<nbasis_old[1]; iy++) {
            for (int iz=0; iz<nbasis_old[2]; iz++) {
                switch (bc_periodic[0]) {
                case 0:
                    if (ix==0)                    { vx[0]=0    ; vx[1]=0    ; vx[2]=1    ; vx[3]=1./2.; }
                    else if (ix==1)               { vx[0]=0    ; vx[1]=1./2.; vx[2]=3./4.; vx[3]=1./4.; }
                    else if (ix<nelem_old[0])     { vx[0]=1./4.; vx[1]=3./4.; vx[2]=3./4.; vx[3]=1./4.; }
                    else if (ix==nelem_old[0])    { vx[0]=1./4.; vx[1]=3./4.; vx[2]=1./2.; vx[3]=0    ; }
                    else if (ix==nelem_old[0]+1)  { vx[0]=1./2.; vx[1]=1    ; vx[2]=0    ; vx[3]=0    ; }
                    break;
                case 1:
                    if (ix==0)                    { vx[0]=0    ; vx[1]=0    ; vx[2]=3./4.; vx[3]=1./4.; }
                    else if (ix==1)               { vx[0]=1./4.; vx[1]=3./4.; vx[2]=3./4.; vx[3]=1./4.; }
                    else if (ix<nelem_old[0])     { vx[0]=1./4.; vx[1]=3./4.; vx[2]=3./4.; vx[3]=1./4.; }
                    else if (ix==nelem_old[0])    { vx[0]=1./4.; vx[1]=3./4.; vx[2]=3./4.; vx[3]=1./4.; }
                    else if (ix==nelem_old[0]+1)  { vx[0]=1./4.; vx[1]=3./4 ; vx[2]=0    ; vx[3]=0    ; }
                    break;
                }
                switch (bc_periodic[1]) {
                case 0:
                    if (iy==0)                    { vy[0]=0    ; vy[1]=0    ; vy[2]=1    ; vy[3]=1./2.; }
                    else if (iy==1)               { vy[0]=0    ; vy[1]=1./2.; vy[2]=3./4.; vy[3]=1./4.; }
                    else if (iy<nelem_old[1])     { vy[0]=1./4.; vy[1]=3./4.; vy[2]=3./4.; vy[3]=1./4.; }
                    else if (iy==nelem_old[1])    { vy[0]=1./4.; vy[1]=3./4.; vy[2]=1./2.; vy[3]=0    ; }
                    else if (iy==nelem_old[1]+1)  { vy[0]=1./2.; vy[1]=1    ; vy[2]=0    ; vy[3]=0    ; }
                    break;
                case 1:
                    if (iy==0)                    { vy[0]=0    ; vy[1]=0    ; vy[2]=3./4.; vy[3]=1./4.; }
                    else if (iy==1)               { vy[0]=1./4.; vy[1]=3./4.; vy[2]=3./4.; vy[3]=1./4.; }
                    else if (iy<nelem_old[1])     { vy[0]=1./4.; vy[1]=3./4.; vy[2]=3./4.; vy[3]=1./4.; }
                    else if (iy==nelem_old[1])    { vy[0]=1./4.; vy[1]=3./4.; vy[2]=3./4.; vy[3]=1./4.; }
                    else if (iy==nelem_old[1]+1)  { vy[0]=1./4.; vy[1]=3./4.; vy[2]=0    ; vy[3]=0    ; }
                    break;
                }
                switch (bc_periodic[2]) {
                case 0:
                    if (iz==0)                    { vz[0]=0    ; vz[1]=0    ; vz[2]=1    ; vz[3]=1./2.; }
                    else if (iz==1)               { vz[0]=0    ; vz[1]=1./2.; vz[2]=3./4.; vz[3]=1./4.; }
                    else if (iz<nelem_old[2])     { vz[0]=1./4.; vz[1]=3./4.; vz[2]=3./4.; vz[3]=1./4.; }
                    else if (iz==nelem_old[2])    { vz[0]=1./4.; vz[1]=3./4.; vz[2]=1./2.; vz[3]=0    ; }
                    else if (iz==nelem_old[2]+1)  { vz[0]=1./2.; vz[1]=1    ; vz[2]=0    ; vz[3]=0    ; }
                    break;
                case 1:
                    if (iz==0)                    { vz[0]=0    ; vz[1]=0    ; vz[2]=3./4.; vz[3]=1./4.; }
                    else if (iz==1)               { vz[0]=1./4.; vz[1]=3./4.; vz[2]=3./4.; vz[3]=1./4.; }
                    else if (iz<nelem_old[2])     { vz[0]=1./4.; vz[1]=3./4.; vz[2]=3./4.; vz[3]=1./4.; }
                    else if (iz==nelem_old[2])    { vz[0]=1./4.; vz[1]=3./4.; vz[2]=3./4.; vz[3]=1./4.; }
                    else if (iz==nelem_old[2]+1)  { vz[0]=1./4.; vz[1]=3./4.; vz[2]=0    ; vz[3]=0    ; }
                    break;
                }
                for (int jx=0; jx<4; jx++) {
                for (int jy=0; jy<4; jy++) {
                for (int jz=0; jz<4; jz++) {
                if (vx[jx]*vy[jy]*vz[jz] != 0) 
                { Ui_new[(nbasis_new[1]*nbasis_new[2])*(2*ix-2+jx)+nbasis_new[2]*(2*iy-2+jy)+(2*iz-2+jz)]+=Ui_old[(nbasis_old[1]*nbasis_old[2])*ix+nbasis_old[2]*iy+iz]*vx[jx]*vy[jy]*vz[jz]; }
                }}}
            }}}
            for (int i=0; i<nbasis_new[0]*nbasis_new[1]*nbasis_new[2]; i++) { Ui_old[i]=Ui_new[i]; }
            for (int idim=0; idim<ndim; idim++) { nelem_old[idim]*=2; }
            for (int idim=0; idim<ndim; idim++) { nbasis_old[idim]=nelem_old[idim]+porder; }
        }
        for (int idim=0; idim<ndim; idim++) { nbasis[idim]=nelem_old[idim]+(bc_periodic[idim]==0)*porder; }
        for (int i=0; i<nbasis[0]; i++) { 
        for (int j=0; j<nbasis[1]; j++) { 
        for (int k=0; k<nbasis[2]; k++) { 
            U_universal[ndof*(nbasis[1]*nbasis[2]*i+nbasis[2]*j+k)+idof]=Ui_old[nbasis_old[1]*nbasis_old[2]*i+nbasis_old[2]*j+k];
        }}}
    } // for ndof
    free(Ui_new);
    free(Ui_old);
}
