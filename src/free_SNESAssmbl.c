#include "petscsnes.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "snes.h"

void free_SNESAssmbl(void *snesassmbl_)
{
    SNESAssmbl *snesassmbl = (SNESAssmbl*)snesassmbl_;
    
    VecDestroy(&snesassmbl->U);
    VecDestroy(&snesassmbl->Residual_);
    MatDestroy(&snesassmbl->Tangent_);
    SNESDestroy(&snesassmbl->snes);
}
