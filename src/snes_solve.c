#include "petscsnes.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "snes.h"
#include "mathutil.h"

void snes_solve(void *snesassmbl_)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    SNESAssmbl *snesassmbl = (SNESAssmbl*)snesassmbl_;
    clock_t time0;

    KSP ksp; SNESGetKSP(snesassmbl->snes,&ksp);
    int                 ksp_niter;
    SNESConvergedReason snes_reason;

    // ---- subiterate
    double Fnorm_new,Fnorm_old=999.;
    for (int isubiter=0; isubiter<10000; isubiter++)
    {
        time0=clock();
        SNESSolve(snesassmbl->snes,NULL,snesassmbl->U);
        KSPGetIterationNumber(ksp,&ksp_niter);
        SNESGetConvergedReason(snesassmbl->snes,&snes_reason);
        SNESGetFunctionNorm(snesassmbl->snes,&Fnorm_new); if ( double_abs(Fnorm_old-Fnorm_new)/Fnorm_old < 1.e-6 ) { snes_reason=SNES_DIVERGED_LINE_SEARCH; } Fnorm_old=Fnorm_new;
        if (rank==0) { printf("%6d: Fnorm = %e, #lin.ite.= %6d, (%08.2f[min])\n",isubiter+1,Fnorm_new,ksp_niter,((float)(clock()-time0))/CLOCKS_PER_SEC/60.0); }
        if (snes_reason != SNES_DIVERGED_MAX_IT) { break; }
    }
    if (snes_reason<=0 && rank==0) { printf("%s.\n",SNESConvergedReasons[snes_reason]); }
    if (snes_reason>0) { snesassmbl->converged=1; } else { snesassmbl->converged=0; }
}
