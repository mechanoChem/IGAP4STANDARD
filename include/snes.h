#if !defined (SNESW)
#define SNESW

#include "petscsnes.h"

typedef PetscErrorCode (*snes_Residual_ptr)(SNES,Vec,Vec,void *);
typedef PetscErrorCode (*snes_Tangent_ptr)(SNES,Vec,Mat,Mat,void *);

typedef struct {
    snes_Residual_ptr snes_Residual;
    snes_Tangent_ptr  snes_Tangent;
    Vec               Residual_;
    Mat               Tangent_;

    SNES              snes;
    Vec               U;
    int               converged;    
} SNESAssmbl;
void free_SNESAssmbl(void *snesassmbl_);

void snes_init(void *snesassmbl_, void *app);
void snes_solve(void *snesassmbl_);

void initguess(double *U_hist, int nselfIdx, void *snesassmbl_);
#endif
