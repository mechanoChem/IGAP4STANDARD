#if !defined (PHYS)
#define PHYS

typedef void (*residual_ptr)(double *residual_, double *u, double *par);
typedef void (*tangent_ptr)(double *tangent_, double *u, double *par);
typedef void (*field_ptr)(double *field, double *u, double *par);
typedef void (*density_ptr)(double h[], double u[], double *par);
typedef void (*densityb_ptr)(double h[], double u[], double traction, double *par, int order);
typedef void (*assert_par_mat_ptr)(double *par);

typedef struct {
    int ndim;
    int ndof;
    int torder;
    int sorder;
    double *par_mat;
    int ndensity;
    int nfield;

    residual_ptr       residual;
    tangent_ptr        tangent;
    field_ptr          field;
    density_ptr        density;
    densityb_ptr       densityb;
    assert_par_mat_ptr assert_par_mat;
} Phys;

#endif
