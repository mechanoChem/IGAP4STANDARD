#include <math.h>
#include "physge.h"
#include <assert.h>
void ge_assert_par_mat(double *par)
{
assert(2*pow((*(par+0)),2) - 2*pow((*(par+0)),2)*(*(par+1))>0);
assert(1.3333333333333333*(*(par+0)) - 1.3333333333333333*(*(par+0))*(*(par+1))>0);
}
