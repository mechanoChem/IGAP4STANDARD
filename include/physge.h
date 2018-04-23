#if !defined (PHYSGE)
#define PHYSGE

#include "phys.h"

void ge_init(void *phys_);
void ge_residual(double *residual, double *u, double *par);
void ge_tangent(double *tangent, double *u, double *par);
void ge_density(double h[], double u[], double *par);
void ge_densityb(double h[], double u[], double traction, double *par, int order);
void ge_field(double *field, double *u, double *par);
void ge_assert_par_mat(double *par);

#endif
