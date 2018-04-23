#if !defined (IOCOMM)
#define IOCOMM

typedef struct {

    int nDof;
    int *nDof_array;
    int *iDof_displ_array;
    int nDof_global;
    int *universalIdx;

} IOcomm;
void free_IOcomm(void *iocomm_);

#endif
