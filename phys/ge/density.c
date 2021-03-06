#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "physge.h"
void ge_density(double hh[], double u[], double *par)
{
double e[39];
*(e+0)=(-1 + pow(1 + (*(u+1)),2) + pow((*(u+11)),2) + pow((*(u+21)),2))/2.;
*(e+1)=((1 + (*(u+1)))*(*(u+2)) + (*(u+11))*(1 + (*(u+12))) + (*(u+21))*(*(u+22)))/2.;
*(e+2)=((1 + (*(u+1)))*(*(u+3)) + (*(u+11))*(*(u+13)) + (*(u+21))*(1 + (*(u+23))))/2.;
*(e+3)=(-1 + pow((*(u+2)),2) + pow(1 + (*(u+12)),2) + pow((*(u+22)),2))/2.;
*(e+4)=((*(u+2))*(*(u+3)) + (1 + (*(u+12)))*(*(u+13)) + (*(u+22))*(1 + (*(u+23))))/2.;
*(e+5)=(-1 + pow((*(u+3)),2) + pow((*(u+13)),2) + pow(1 + (*(u+23)),2))/2.;
*(e+6)=(2*(1 + (*(u+1)))*(*(u+4)) + 2*(*(u+11))*(*(u+14)) + 2*(*(u+21))*(*(u+24)))/2.;
*(e+7)=(2*(1 + (*(u+1)))*(*(u+5)) + 2*(*(u+11))*(*(u+15)) + 2*(*(u+21))*(*(u+25)))/2.;
*(e+8)=(2*(1 + (*(u+1)))*(*(u+6)) + 2*(*(u+11))*(*(u+16)) + 2*(*(u+21))*(*(u+26)))/2.;
*(e+9)=((*(u+2))*(*(u+4)) + (1 + (*(u+1)))*(*(u+5)) + (1 + (*(u+12)))*(*(u+14)) + (*(u+11))*(*(u+15)) + (*(u+22))*(*(u+24)) + (*(u+21))*(*(u+25)))/2.;
*(e+10)=((*(u+2))*(*(u+5)) + (1 + (*(u+1)))*(*(u+7)) + (1 + (*(u+12)))*(*(u+15)) + (*(u+11))*(*(u+17)) + (*(u+22))*(*(u+25)) + (*(u+21))*(*(u+27)))/2.;
*(e+11)=((*(u+2))*(*(u+6)) + (1 + (*(u+1)))*(*(u+8)) + (1 + (*(u+12)))*(*(u+16)) + (*(u+11))*(*(u+18)) + (*(u+22))*(*(u+26)) + (*(u+21))*(*(u+28)))/2.;
*(e+12)=((*(u+3))*(*(u+4)) + (1 + (*(u+1)))*(*(u+6)) + (*(u+13))*(*(u+14)) + (*(u+11))*(*(u+16)) + (1 + (*(u+23)))*(*(u+24)) + (*(u+21))*(*(u+26)))/2.;
*(e+13)=((*(u+3))*(*(u+5)) + (1 + (*(u+1)))*(*(u+8)) + (*(u+13))*(*(u+15)) + (*(u+11))*(*(u+18)) + (1 + (*(u+23)))*(*(u+25)) + (*(u+21))*(*(u+28)))/2.;
*(e+14)=((*(u+3))*(*(u+6)) + (1 + (*(u+1)))*(*(u+9)) + (*(u+13))*(*(u+16)) + (*(u+11))*(*(u+19)) + (1 + (*(u+23)))*(*(u+26)) + (*(u+21))*(*(u+29)))/2.;
*(e+15)=(2*(*(u+2))*(*(u+5)) + 2*(1 + (*(u+12)))*(*(u+15)) + 2*(*(u+22))*(*(u+25)))/2.;
*(e+16)=(2*(*(u+2))*(*(u+7)) + 2*(1 + (*(u+12)))*(*(u+17)) + 2*(*(u+22))*(*(u+27)))/2.;
*(e+17)=(2*(*(u+2))*(*(u+8)) + 2*(1 + (*(u+12)))*(*(u+18)) + 2*(*(u+22))*(*(u+28)))/2.;
*(e+18)=((*(u+3))*(*(u+5)) + (*(u+2))*(*(u+6)) + (*(u+13))*(*(u+15)) + (1 + (*(u+12)))*(*(u+16)) + (1 + (*(u+23)))*(*(u+25)) + (*(u+22))*(*(u+26)))/2.;
*(e+19)=((*(u+3))*(*(u+7)) + (*(u+2))*(*(u+8)) + (*(u+13))*(*(u+17)) + (1 + (*(u+12)))*(*(u+18)) + (1 + (*(u+23)))*(*(u+27)) + (*(u+22))*(*(u+28)))/2.;
*(e+20)=((*(u+3))*(*(u+8)) + (*(u+2))*(*(u+9)) + (*(u+13))*(*(u+18)) + (1 + (*(u+12)))*(*(u+19)) + (1 + (*(u+23)))*(*(u+28)) + (*(u+22))*(*(u+29)))/2.;
*(e+21)=(2*(*(u+3))*(*(u+6)) + 2*(*(u+13))*(*(u+16)) + 2*(1 + (*(u+23)))*(*(u+26)))/2.;
*(e+22)=(2*(*(u+3))*(*(u+8)) + 2*(*(u+13))*(*(u+18)) + 2*(1 + (*(u+23)))*(*(u+28)))/2.;
*(e+23)=(2*(*(u+3))*(*(u+9)) + 2*(*(u+13))*(*(u+19)) + 2*(1 + (*(u+23)))*(*(u+29)))/2.;
*(e+24)=((*(e+0)) + (*(e+3)) + (*(e+5)))/sqrt(3);
*(e+25)=((*(e+0)) - (*(e+3)))/sqrt(2);
*(e+26)=((*(e+0)) + (*(e+3)) - 2*(*(e+5)))/sqrt(6);
*(e+27)=((*(e+6)) + (*(e+15)) + (*(e+21)))/sqrt(3);
*(e+28)=((*(e+7)) + (*(e+16)) + (*(e+22)))/sqrt(3);
*(e+29)=((*(e+8)) + (*(e+17)) + (*(e+23)))/sqrt(3);
*(e+30)=((*(e+6)) - (*(e+15)))/sqrt(2);
*(e+31)=((*(e+7)) - (*(e+16)))/sqrt(2);
*(e+32)=((*(e+8)) - (*(e+17)))/sqrt(2);
*(e+33)=((*(e+6)) + (*(e+15)) - 2*(*(e+21)))/sqrt(6);
*(e+34)=((*(e+7)) + (*(e+16)) - 2*(*(e+22)))/sqrt(6);
*(e+35)=((*(e+8)) + (*(e+17)) - 2*(*(e+23)))/sqrt(6);
*(e+36)=-((1.3333333333333333*(*(par+0)) - 1.3333333333333333*(*(par+0))*(*(par+1)))/(pow((*(par+0)),2)*(pow((*(par+0)),2) - 2*pow((*(par+0)),2)*(*(par+1)) - (*(par+0))*(1.3333333333333333*(*(par+0)) - 1.3333333333333333*(*(par+0))*(*(par+1))))));
*(e+37)=(2*(*(par+1)))/(pow((*(par+0)),2) - 2*pow((*(par+0)),2)*(*(par+1)) - (*(par+0))*(1.3333333333333333*(*(par+0)) - 1.3333333333333333*(*(par+0))*(*(par+1))));
*(e+38)=-(1/(pow((*(par+0)),2)*(pow((*(par+0)),2) - 2*pow((*(par+0)),2)*(*(par+1)) - (*(par+0))*(1.3333333333333333*(*(par+0)) - 1.3333333333333333*(*(par+0))*(*(par+1))))));
// PI
hh[0]=(*(e+26))*(-3*pow((*(e+25)),2) + pow((*(e+26)),2))*(*(e+36)) + (pow((*(e+25)),2) + pow((*(e+26)),2))*(*(e+37)) + pow(pow((*(e+25)),2) + pow((*(e+26)),2),2)*(*(e+38)) + pow((*(e+24)),2)*(*(par+2)) + (pow((*(e+1)),2) + pow((*(e+2)),2) + pow((*(e+4)),2))*(*(par+3)) + (pow((*(e+30)),2) + pow((*(e+31)),2) + pow((*(e+32)),2) + pow((*(e+33)),2) + pow((*(e+34)),2) + pow((*(e+35)),2))*pow((*(par+4)),2);
// S
hh[1]=2*(sqrt(2)*(*(e+25))*(-3*(*(e+26))*(*(e+36)) + (*(e+37)) + 2*(pow((*(e+25)),2) + pow((*(e+26)),2))*(*(e+38))) + (pow((*(e+25)),2)*(-3*(*(e+36)) + 4*(*(e+26))*(*(e+38))) + (*(e+26))*(3*(*(e+26))*(*(e+36)) + 2*(*(e+37)) + 4*pow((*(e+26)),2)*(*(e+38))))/sqrt(6) + (2*(*(e+24))*(*(par+2)))/sqrt(3));
hh[2]=4*(*(e+1))*(*(par+3));
hh[3]=4*(*(e+2))*(*(par+3));
hh[4]=2*(-(sqrt(2)*(*(e+25))*(-3*(*(e+26))*(*(e+36)) + (*(e+37)) + 2*(pow((*(e+25)),2) + pow((*(e+26)),2))*(*(e+38)))) + (pow((*(e+25)),2)*(-3*(*(e+36)) + 4*(*(e+26))*(*(e+38))) + (*(e+26))*(3*(*(e+26))*(*(e+36)) + 2*(*(e+37)) + 4*pow((*(e+26)),2)*(*(e+38))))/sqrt(6) + (2*(*(e+24))*(*(par+2)))/sqrt(3));
hh[5]=4*(*(e+4))*(*(par+3));
hh[6]=(2*(-(sqrt(2)*(pow((*(e+25)),2)*(-3*(*(e+36)) + 4*(*(e+26))*(*(e+38))) + (*(e+26))*(3*(*(e+26))*(*(e+36)) + 2*(*(e+37)) + 4*pow((*(e+26)),2)*(*(e+38))))) + 2*(*(e+24))*(*(par+2))))/sqrt(3);
// volume fraction
hh[7]=(double)(e[25]>0 && e[26]>-0.5*e[25]);
hh[8]=(double)(e[25]<0 && e[26]> 0.5*e[25]);
hh[9]=(double)(e[26]< 0.5*e[25] && e[26]<-0.5*e[25]);
// volume fraction: energy based
hh[10]=(double)(hh[0]<-0.5 && e[25]>0 && e[26]>0);
hh[11]=(double)(hh[0]<-0.5 && e[25]<0 && e[26]>0);
hh[12]=(double)(hh[0]<-0.5 && e[26]<0);

}
// boundary
void ge_densityb(double hh[], double u[], double traction, double *par, int order)
{

hh[0]=u[0]*traction;
hh[1]=0;
hh[2]=0;
hh[3]=0;
hh[4]=0;
hh[5]=0;
hh[6]=0;
hh[7]=0;
hh[8]=0;
hh[9]=0;
hh[10]=0;
hh[11]=0;
hh[12]=0;
}
