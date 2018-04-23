#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "nurbs.h"
#include "mathutil.h"
double evalN(double *kV, int nk, int k, int i, int p, double xi)
/*

evaluates N_{i,p}(\xi)　using recursion.

k  : order of derivative
kV : knotVector
nk : # of knots
p  : porder
i  : ibasis
xi : coordinate of the quadrature point

*/
{

if (k==0)
{

//if (i==(nk-p-1-1) && xi==kV[nk-1]) // due to imperfection of the definition.
if (i==(nk-p-1-1) && double_abs(xi-kV[nk-1])<1.e-15 && double_abs(kV[nk-1]-kV[nk-2])<1.e-15) // temp, due to imperfection of the definition.
{ return 1.0; }
else if (p==0)
{ return ((xi>=kV[i] && kV[i+1]>xi) ? 1.0:0.0); }
else if (kV[i+p]==kV[i] && kV[i+p+1]==kV[i+1])
{ return 0.0; }
else if (kV[i+p]==kV[i])
{ return (kV[i+p+1]-xi)/(kV[i+p+1]-kV[i+1])*evalN(kV,nk,k,i+1,p-1,xi); }
else if (kV[i+p+1]==kV[i+1])
{ return (xi-kV[i])/(kV[i+p]-kV[i])*evalN(kV,nk,k,i,p-1,xi); }
else
{ return (xi-kV[i])/(kV[i+p]-kV[i])*evalN(kV,nk,k,i,p-1,xi)+(kV[i+p+1]-xi)/(kV[i+p+1]-kV[i+1])*evalN(kV,nk,k,i+1,p-1,xi); }

}
else if ( k>0 && k<=p )
{
double alpha[k+1][k+1];
int icol, irow;
alpha[0][0]=1.0;
for (irow=1; irow<k+1; irow++)
{
    alpha[irow][0]=( (kV[i+p-irow+1]==kV[i]) ? 0.0:alpha[irow-1][0]/(kV[i+p-irow+1]-kV[i]) );
    for (icol=1; icol<irow; icol++)
    {
        alpha[irow][icol]=( (kV[i+p+icol-irow+1]==kV[i+icol]) ? 0.0:(alpha[irow-1][icol]-alpha[irow-1][icol-1])/(kV[i+p+icol-irow+1]-kV[i+icol]) );
    }
    alpha[irow][irow]=( (kV[i+p+1]==kV[i+irow]) ? 0.0:-alpha[irow-1][irow-1]/(kV[i+p+1]-kV[i+irow]) );
}

double sum=0.0;
for (int j=0; j<k+1; j++){ sum+=alpha[k][j]*evalN(kV,nk,0,i+j,p-k,xi); }
return sum*factorial(p)/factorial(p-k);
}
else
{ printf("error in evalN.c: 0 =< k <= p has to be satisfied.\n"); return 0.0; }

} // evalN
