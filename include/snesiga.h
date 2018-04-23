#if !defined (SNESIGA)
#define SNESIGA

#include "petscsnes.h"
#include "iga.h"
#include "snes.h"

PetscErrorCode snesiga_Residual(SNES snes, Vec U_, Vec Residual_, void *app_);
PetscErrorCode snesiga_Tangent(SNES snes, Vec U_, Mat Tangent_, Mat Pmat_, void *app_);

void snesiga_init(void* app_, void* phys_, void* nurbs_, void* bc_, void* snesassmbl_, void* iocomm_);

#endif
