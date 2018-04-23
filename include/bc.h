#if !defined (BC_)
#define BC_

typedef struct{

    int *bc_type;
    double *par_periodic;
    double *par_dirichlet;
    double *par_neumann;

} BC;

void bc_init(void* bc_);
void neumann(double *neumann_, double par_neumann[], int face_id, int bc_order, int idof, double xq[]);

#endif
