#if !defined (IGA_)
#define IGA_

#include "phys.h"
#include "nurbs.h"
#include "bc.h"

typedef struct{

    Phys  *phys;
    Nurbs *nurbs;
    BC    *bc;

    double *U_hist;
    int nselfIdx;
    int *selfIdx;
    int nsendIdx;
    int *sendIdx;
    int *sendPtr;
    int nrecvIdx;
    int *recvIdx;
    int *recvPtr;
    int *globalIdx;

} IGA;
void free_IGA(void *iga_);


void iga_intdensity(double HH_[], void *app_);
//void init_linadd(double *U_, double par_dirichlet[], double par_dirichlet_temp[], void *app_);
void iga_dirichlet(double *dirichlet, double par[], int face_id, int bc_order, int idof, double *knotVector_[], int nknot_[], int nbasis_[], int ndim, int porder);
void iga_dirichletinit(double *U_, void *app_);
void appl_init(void *app_, void *iocomm_);

// plot
void iga_plotfield(const char *fname, void *app_, int uid[], double uscale, int nppe, int figsize[], double relplot[],  double rotate[], double ranges[], double crange[], const char *sch, const char *opt);



#endif
