#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "iocomm.h"

void free_IOcomm(void *iocomm_)
{
    IOcomm *iocomm = (IOcomm*)iocomm_;
   
    free(iocomm->nDof_array);
    free(iocomm->iDof_displ_array);
    free(iocomm->universalIdx);
}
