#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "mathutil.h"
void matinv_dot_vec(double mat[], double vec[], int n)
{
    matinv(mat,n);
    double vec_[n]; for (int i=0; i<n; i++) { vec_[i]=0.0; }
    int ia=0;
    for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
        vec_[i]+=mat[ia++]*vec[j];
    }}
    for (int i=0; i<n; i++) { vec[i]=vec_[i]; }
}
void matinv(double mat[], int n)
{
double mat_[n*n];
double idet;
switch(n)
{
case 2:
    idet=1./(-((*(mat+1))*(*(mat+2))) + (*(mat+0))*(*(mat+3)));
    mat_[0]=idet*((*(mat+3)));
    mat_[1]=idet*(-(*(mat+1)));
    mat_[2]=idet*(-(*(mat+2)));
    mat_[3]=idet*((*(mat+0)));
    break;
case 3:
    idet=1./(-((*(mat+2))*(*(mat+4))*(*(mat+6))) + (*(mat+1))*(*(mat+5))*(*(mat+6)) + (*(mat+2))*(*(mat+3))*(*(mat+7)) - (*(mat+0))*(*(mat+5))*(*(mat+7)) - (*(mat+1))*(*(mat+3))*(*(mat+8)) + (*(mat+0))*(*(mat+4))*(*(mat+8)));
    mat_[0]=idet*(-((*(mat+5))*(*(mat+7))) + (*(mat+4))*(*(mat+8)));
    mat_[1]=idet*((*(mat+2))*(*(mat+7)) - (*(mat+1))*(*(mat+8)));
    mat_[2]=idet*(-((*(mat+2))*(*(mat+4))) + (*(mat+1))*(*(mat+5)));
    mat_[3]=idet*((*(mat+5))*(*(mat+6)) - (*(mat+3))*(*(mat+8)));
    mat_[4]=idet*(-((*(mat+2))*(*(mat+6))) + (*(mat+0))*(*(mat+8)));
    mat_[5]=idet*((*(mat+2))*(*(mat+3)) - (*(mat+0))*(*(mat+5)));
    mat_[6]=idet*(-((*(mat+4))*(*(mat+6))) + (*(mat+3))*(*(mat+7)));
    mat_[7]=idet*((*(mat+1))*(*(mat+6)) - (*(mat+0))*(*(mat+7)));
    mat_[8]=idet*(-((*(mat+1))*(*(mat+3))) + (*(mat+0))*(*(mat+4)));
    break;
case 4:
    idet=1./((*(mat+3))*(*(mat+6))*(*(mat+9))*(*(mat+12)) - (*(mat+2))*(*(mat+7))*(*(mat+9))*(*(mat+12)) - (*(mat+3))*(*(mat+5))*(*(mat+10))*(*(mat+12)) + (*(mat+1))*(*(mat+7))*(*(mat+10))*(*(mat+12)) + (*(mat+2))*(*(mat+5))*(*(mat+11))*(*(mat+12)) - (*(mat+1))*(*(mat+6))*(*(mat+11))*(*(mat+12)) - (*(mat+3))*(*(mat+6))*(*(mat+8))*(*(mat+13)) + (*(mat+2))*(*(mat+7))*(*(mat+8))*(*(mat+13)) + (*(mat+3))*(*(mat+4))*(*(mat+10))*(*(mat+13)) - (*(mat+0))*(*(mat+7))*(*(mat+10))*(*(mat+13)) - (*(mat+2))*(*(mat+4))*(*(mat+11))*(*(mat+13)) + (*(mat+0))*(*(mat+6))*(*(mat+11))*(*(mat+13)) + (*(mat+3))*(*(mat+5))*(*(mat+8))*(*(mat+14)) - (*(mat+1))*(*(mat+7))*(*(mat+8))*(*(mat+14)) - (*(mat+3))*(*(mat+4))*(*(mat+9))*(*(mat+14)) + (*(mat+0))*(*(mat+7))*(*(mat+9))*(*(mat+14)) + (*(mat+1))*(*(mat+4))*(*(mat+11))*(*(mat+14)) - (*(mat+0))*(*(mat+5))*(*(mat+11))*(*(mat+14)) - (*(mat+2))*(*(mat+5))*(*(mat+8))*(*(mat+15)) + (*(mat+1))*(*(mat+6))*(*(mat+8))*(*(mat+15)) + (*(mat+2))*(*(mat+4))*(*(mat+9))*(*(mat+15)) - (*(mat+0))*(*(mat+6))*(*(mat+9))*(*(mat+15)) - (*(mat+1))*(*(mat+4))*(*(mat+10))*(*(mat+15)) + (*(mat+0))*(*(mat+5))*(*(mat+10))*(*(mat+15)));
    mat_[0]=idet*(-((*(mat+7))*(*(mat+10))*(*(mat+13))) + (*(mat+6))*(*(mat+11))*(*(mat+13)) + (*(mat+7))*(*(mat+9))*(*(mat+14)) - (*(mat+5))*(*(mat+11))*(*(mat+14)) - (*(mat+6))*(*(mat+9))*(*(mat+15)) + (*(mat+5))*(*(mat+10))*(*(mat+15)));
    mat_[1]=idet*((*(mat+3))*(*(mat+10))*(*(mat+13)) - (*(mat+2))*(*(mat+11))*(*(mat+13)) - (*(mat+3))*(*(mat+9))*(*(mat+14)) + (*(mat+1))*(*(mat+11))*(*(mat+14)) + (*(mat+2))*(*(mat+9))*(*(mat+15)) - (*(mat+1))*(*(mat+10))*(*(mat+15)));
    mat_[2]=idet*(-((*(mat+3))*(*(mat+6))*(*(mat+13))) + (*(mat+2))*(*(mat+7))*(*(mat+13)) + (*(mat+3))*(*(mat+5))*(*(mat+14)) - (*(mat+1))*(*(mat+7))*(*(mat+14)) - (*(mat+2))*(*(mat+5))*(*(mat+15)) + (*(mat+1))*(*(mat+6))*(*(mat+15)));
    mat_[3]=idet*((*(mat+3))*(*(mat+6))*(*(mat+9)) - (*(mat+2))*(*(mat+7))*(*(mat+9)) - (*(mat+3))*(*(mat+5))*(*(mat+10)) + (*(mat+1))*(*(mat+7))*(*(mat+10)) + (*(mat+2))*(*(mat+5))*(*(mat+11)) - (*(mat+1))*(*(mat+6))*(*(mat+11)));
    mat_[4]=idet*((*(mat+7))*(*(mat+10))*(*(mat+12)) - (*(mat+6))*(*(mat+11))*(*(mat+12)) - (*(mat+7))*(*(mat+8))*(*(mat+14)) + (*(mat+4))*(*(mat+11))*(*(mat+14)) + (*(mat+6))*(*(mat+8))*(*(mat+15)) - (*(mat+4))*(*(mat+10))*(*(mat+15)));
    mat_[5]=idet*(-((*(mat+3))*(*(mat+10))*(*(mat+12))) + (*(mat+2))*(*(mat+11))*(*(mat+12)) + (*(mat+3))*(*(mat+8))*(*(mat+14)) - (*(mat+0))*(*(mat+11))*(*(mat+14)) - (*(mat+2))*(*(mat+8))*(*(mat+15)) + (*(mat+0))*(*(mat+10))*(*(mat+15)));
    mat_[6]=idet*((*(mat+3))*(*(mat+6))*(*(mat+12)) - (*(mat+2))*(*(mat+7))*(*(mat+12)) - (*(mat+3))*(*(mat+4))*(*(mat+14)) + (*(mat+0))*(*(mat+7))*(*(mat+14)) + (*(mat+2))*(*(mat+4))*(*(mat+15)) - (*(mat+0))*(*(mat+6))*(*(mat+15)));
    mat_[7]=idet*(-((*(mat+3))*(*(mat+6))*(*(mat+8))) + (*(mat+2))*(*(mat+7))*(*(mat+8)) + (*(mat+3))*(*(mat+4))*(*(mat+10)) - (*(mat+0))*(*(mat+7))*(*(mat+10)) - (*(mat+2))*(*(mat+4))*(*(mat+11)) + (*(mat+0))*(*(mat+6))*(*(mat+11)));
    mat_[8]=idet*(-((*(mat+7))*(*(mat+9))*(*(mat+12))) + (*(mat+5))*(*(mat+11))*(*(mat+12)) + (*(mat+7))*(*(mat+8))*(*(mat+13)) - (*(mat+4))*(*(mat+11))*(*(mat+13)) - (*(mat+5))*(*(mat+8))*(*(mat+15)) + (*(mat+4))*(*(mat+9))*(*(mat+15)));
    mat_[9]=idet*((*(mat+3))*(*(mat+9))*(*(mat+12)) - (*(mat+1))*(*(mat+11))*(*(mat+12)) - (*(mat+3))*(*(mat+8))*(*(mat+13)) + (*(mat+0))*(*(mat+11))*(*(mat+13)) + (*(mat+1))*(*(mat+8))*(*(mat+15)) - (*(mat+0))*(*(mat+9))*(*(mat+15)));
    mat_[10]=idet*(-((*(mat+3))*(*(mat+5))*(*(mat+12))) + (*(mat+1))*(*(mat+7))*(*(mat+12)) + (*(mat+3))*(*(mat+4))*(*(mat+13)) - (*(mat+0))*(*(mat+7))*(*(mat+13)) - (*(mat+1))*(*(mat+4))*(*(mat+15)) + (*(mat+0))*(*(mat+5))*(*(mat+15)));
    mat_[11]=idet*((*(mat+3))*(*(mat+5))*(*(mat+8)) - (*(mat+1))*(*(mat+7))*(*(mat+8)) - (*(mat+3))*(*(mat+4))*(*(mat+9)) + (*(mat+0))*(*(mat+7))*(*(mat+9)) + (*(mat+1))*(*(mat+4))*(*(mat+11)) - (*(mat+0))*(*(mat+5))*(*(mat+11)));
    mat_[12]=idet*((*(mat+6))*(*(mat+9))*(*(mat+12)) - (*(mat+5))*(*(mat+10))*(*(mat+12)) - (*(mat+6))*(*(mat+8))*(*(mat+13)) + (*(mat+4))*(*(mat+10))*(*(mat+13)) + (*(mat+5))*(*(mat+8))*(*(mat+14)) - (*(mat+4))*(*(mat+9))*(*(mat+14)));
    mat_[13]=idet*(-((*(mat+2))*(*(mat+9))*(*(mat+12))) + (*(mat+1))*(*(mat+10))*(*(mat+12)) + (*(mat+2))*(*(mat+8))*(*(mat+13)) - (*(mat+0))*(*(mat+10))*(*(mat+13)) - (*(mat+1))*(*(mat+8))*(*(mat+14)) + (*(mat+0))*(*(mat+9))*(*(mat+14)));
    mat_[14]=idet*((*(mat+2))*(*(mat+5))*(*(mat+12)) - (*(mat+1))*(*(mat+6))*(*(mat+12)) - (*(mat+2))*(*(mat+4))*(*(mat+13)) + (*(mat+0))*(*(mat+6))*(*(mat+13)) + (*(mat+1))*(*(mat+4))*(*(mat+14)) - (*(mat+0))*(*(mat+5))*(*(mat+14)));
    mat_[15]=idet*(-((*(mat+2))*(*(mat+5))*(*(mat+8))) + (*(mat+1))*(*(mat+6))*(*(mat+8)) + (*(mat+2))*(*(mat+4))*(*(mat+9)) - (*(mat+0))*(*(mat+6))*(*(mat+9)) - (*(mat+1))*(*(mat+4))*(*(mat+10)) + (*(mat+0))*(*(mat+5))*(*(mat+10)));
    break;
default: assert( n>0 && n<5 );
} //switch
for (int i=0; i<n*n; i++) { mat[i]=mat_[i]; }
}
