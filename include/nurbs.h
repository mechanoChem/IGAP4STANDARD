#if !defined (NURBS)
#define NURBS

typedef struct {

    int porder;
    int *nelem;
    int *nbasis;
    int *nknot;
    double **knotVector;
    int nboun;
    int *boun;
    int nsendPart;
    int *sendPart;
    int nrecvPart;
    int *recvPart;
    int *sendBlock;              // MESH -> IGA (transient storage)
    int *recvBlock;              // MESH -> IGA (transient storage)
    int *universalIdx_temp;      // MESH -> IGA (transient storage)
    int *bc_periodic;            // Global MESH Info (transient storage)
    int *nelem_global;           // Global MESH Info (transient storage)
    double **knotVector_global;  // Global MESH Info (transient storage)

} Nurbs;

void nurbs_init(void *nurbs_);
void nurbs_genuniform(void *nurbs_, int mref, int porder, int bc_periodic_x, int bc_periodic_y, int bc_periodic_z);
void nurbs_knotinsert(double *U_universal, const int nelem[], int nref, int porder, int ndim, int ndof, int bc_periodic[]);

void eval_Bspline(double N[], double *kV[], int p, double xi[]);
double evalN(double *kV, int nk, int k, int i, int p, double xi);

#endif
