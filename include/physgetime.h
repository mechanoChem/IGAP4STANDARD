#if !defined (PHYSGETIME)
#define PHYSGETIME

#include "phys.h"

void getime_init(void *phys_);
void getime_residual(double *residual, double *u, double *par);
void getime_tangent(double *tangent, double *u, double *par);
void getime_density(double h[], double u[], double *par);
void getime_densityb(double h[], double u[], double traction, double *par, int order);
void getime_field(double *field, double *u, double *par);
void getime_assert_par_mat(double *par);

#endif
