#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "physge.h"
void ge_init(void *phys_)
{
    Phys *phys=(Phys*)phys_;
    
    phys->ndim=3;
    phys->ndof=3;
    phys->torder=0;
    phys->sorder=2;
    phys->ndensity=1+6+3+3;
    phys->nfield=33;

    //                    r ,  curv,   e0 , e345 , le,   rho,  cii,        dt,    t
    double par_mat[]={ 0.25,  0.50,  500.,  250.,  0.025, 1.0,  1.0,  1.e-3*(1./2.), 0.0 };
    phys->par_mat=(double*)malloc(9*sizeof(double));
    for (int i=0; i<9; i++) phys->par_mat[i]=par_mat[i];
    
    phys->residual = ge_residual;
    phys->tangent  = ge_tangent;
    phys->field    = ge_field;
    phys->density  = ge_density;
    phys->densityb = ge_densityb;
    phys->assert_par_mat = ge_assert_par_mat;
    
    phys->assert_par_mat(par_mat);
}
