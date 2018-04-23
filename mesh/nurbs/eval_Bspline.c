#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "nurbs.h"

#define MTOL 1.e-15

void eval_Bspline(double N[], double *kV[], int p, double xq[])
/*

requires p==2.
specificaclly for quadrature point evaluations (no end points).



kV[idim]=knotVector[idim]+ielem

Basic formulas:
                  /  1  if kV[i] =< xi < kV[i+1],
   N_{i,0}(xi) = | 
                  \  0  otherwise.

                  xi   - kV[i]                 kV[i+1+p] -   xi
   N_{i,p} =  ----------------- N_{i,p-1}  +  --------------------- N_{i+1,p-1}
               kV[i+p] - kV[i]                 kV[i+1+p] - kV[i+1]

1st-order derivative:
 dN_{i,p}             p                                 p
---------- = ----------------- N_{i,p-1}  -  --------------------- N_{i+1,p-1}
   dxi        kV[i+p] - kV[i]                 kV[i+1+p] - kV[i+1]

2nd-order derivative:
                                                                             d^2N_{i,p}
                                                                            ------------
                                                                               dxi^2 
                                                                        /                 \
                                                                       /                   \
                                                                      /                     \

                        p         dN_{i,p-1}                                                              p           dN_{i+1,p-1}
               ----------------- ---------                                       -             --------------------- ------------
                kV[i+p] - kV[i]     dxi                                                         kV[i+1+p] - kV[i+1]       dxi
                                   /   \                                                                                 /   \
                                  /     \                                                                               /     \
                                 /       \                                                                             /       \
                 p-1                             p-1                                           p-1                                   p-1
       ------------------- N_{i,p-2}  -  ------------------- N_{i+1,p-2}               ------------------- N_{i+1,p-2}  -  --------------------- N_{i+2,p-2}
        kV[i+p-1] - kV[i]                 kV[i+p] - kV[i+1]                             kV[i+p] - kV[i+1]                   kV[i+1+p] - kV[i+2]


supp(N_{i,p}) = supp(dN_{i,p}) = supp(d^2N_{i,p}) = [ kV[i] , kV[i+p+1] ]. (by induction)
Backward dependency:
N_{i,p}
N_{i,p-1} N_{i+1,p-1}
N_{i,p-2} N_{i+1,p-2} N_{i+2,p-2}
*/
{
    if (p!=2) exit(EXIT_FAILURE);
int ndim=3;
int nddim=10; 
    int ia,ibpe;
    int pplus1=p+1;
    double Np[(3*pplus1)*ndim];
    double N1[pplus1+1];
    double N0[pplus1+2]; for (int i=0; i<pplus1+2; i++) { N0[i]=0; } N0[p]=1;
    for (int idim=0; idim<ndim; idim++)
    {
        for (int i=0; i<pplus1+1; i++) { N1[i]                         =(xq[idim]-kV[idim][i])/(kV[idim][i+1]-kV[idim][i]+MTOL)*N0[i]+(kV[idim][i+2]-xq[idim])/(kV[idim][i+2]-kV[idim][i+1]+MTOL)*N0[i+1]; }
        for (int i=0; i<pplus1; i++)   { Np[(3*pplus1)*idim+i]         =(xq[idim]-kV[idim][i])/(kV[idim][i+2]-kV[idim][i]+MTOL)*N1[i]+(kV[idim][i+3]-xq[idim])/(kV[idim][i+3]-kV[idim][i+1]+MTOL)*N1[i+1]; }
        for (int i=0; i<pplus1; i++)   { Np[(3*pplus1)*idim+pplus1+i]  =                   2.0/(kV[idim][i+2]-kV[idim][i]+MTOL)*N1[i]-                     2.0/(kV[idim][i+3]-kV[idim][i+1]+MTOL)*N1[i+1]; }
        for (int i=0; i<pplus1+1; i++) { N1[i]                         =                   1.0/(kV[idim][i+1]-kV[idim][i]+MTOL)*N0[i]-                     1.0/(kV[idim][i+2]-kV[idim][i+1]+MTOL)*N0[i+1]; }
        for (int i=0; i<pplus1; i++)   { Np[(3*pplus1)*idim+2*pplus1+i]=                   2.0/(kV[idim][i+2]-kV[idim][i]+MTOL)*N1[i]-                     2.0/(kV[idim][i+3]-kV[idim][i+1]+MTOL)*N1[i+1]; }
    }

        
        for (int i=0; i<pplus1; i++){
        for (int j=0; j<pplus1; j++){
        for (int k=0; k<pplus1; k++){
            ibpe=pplus1*pplus1*i+pplus1*j+k;
            ia=nddim*ibpe;
            N[ia+0]=Np[i]         *Np[3*pplus1+j]         *Np[6*pplus1+k];
            N[ia+1]=Np[pplus1+i]  *Np[3*pplus1+j]         *Np[6*pplus1+k];
            N[ia+2]=Np[i]         *Np[3*pplus1+pplus1+j]  *Np[6*pplus1+k];
            N[ia+3]=Np[i]         *Np[3*pplus1+j]         *Np[6*pplus1+pplus1+k];
            N[ia+4]=Np[2*pplus1+i]*Np[3*pplus1+j]         *Np[6*pplus1+k];
            N[ia+5]=Np[pplus1+i]  *Np[3*pplus1+pplus1+j]  *Np[6*pplus1+k];
            N[ia+6]=Np[pplus1+i]  *Np[3*pplus1+j]         *Np[6*pplus1+pplus1+k];
            N[ia+7]=Np[i]         *Np[3*pplus1+2*pplus1+j]*Np[6*pplus1+k];
            N[ia+8]=Np[i]         *Np[3*pplus1+pplus1+j]  *Np[6*pplus1+pplus1+k];
            N[ia+9]=Np[i]         *Np[3*pplus1+j]         *Np[6*pplus1+2*pplus1+k];
        }}}
}
