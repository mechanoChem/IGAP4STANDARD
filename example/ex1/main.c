static char help[] = "gradelasttime on unitcube.\n";
#define PI    3.141592653589793
#define EXP1  2.718281828459045

#include "petscsnes.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "physgetime.h"
#include "snesiga.h"
#include "iocomm.h"
#include "ioutil.h"
#include "mathutil.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    PetscInitialize(&argc,&argv,(char*)0,help);


    Phys phys[1];
    getime_init(phys+0);

    
    Nurbs nurbs[1];
    nurbs_genuniform(nurbs+0,3,2,0,0,0);
    nurbs_init(nurbs+0);


    BC bc[1];
    bc_init(bc+0);


    IGA app;
    SNESAssmbl snesassmbl;
    IOcomm     iocomm;
    snesiga_init(&app,phys,nurbs,bc,&snesassmbl,&iocomm);


    snes_init(&snesassmbl,&app);

    double HH[phys[0].ndensity];
    double *U_hist =app.U_hist;
    int    nselfIdx=app.nselfIdx;

    for (int i=0; i<nselfIdx; i++) { U_hist[nselfIdx+i]=1.e-4*sin(i); }
    initguess(U_hist,nselfIdx,&snesassmbl);

    snes_solve(&snesassmbl);
    if (snesassmbl.converged)
    {
        for (int ihist=phys[0].torder; ihist>0; ihist--) for (int i=0; i<nselfIdx; i++) U_hist[nselfIdx*ihist+i]=U_hist[nselfIdx*(ihist-1)+i];
        iga_intdensity(HH,&app);
        mpi_printf("PI = %e.\n\n",HH[0]);
    }

    free_IGA(&app);
    free_SNESAssmbl(&snesassmbl);
    free_IOcomm(&iocomm);
    PetscFinalize();

    return 0;
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          BOUNDARY CONDITION        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void bc_init(void* bc_)
{
    BC *bc = (BC*)bc_;

    double par_periodic[] ={0.0,0.0,0.0,
                            0.0,0.0,0.0,
                            0.0,0.0,0.0 };                  // F-periodic B.C.s: x=X+(F0.X+u)
    // Dirichlet/Neumann B.C.s
    //              ux,   uy,   uz, ux_n, uy_n, uz_n,       //                     ______________________
    int bc_type[]={  0,    0,    0,    1,    1,    1,       // Face 0: YZ         /                     /|
                     0,    0,    0,    1,    1,    1,       // Face 1: YZ        /          5          / |
                     0,    0,    0,    1,    1,    1,       // Face 2: ZX       /_____________________/  |
                     0,    0,    0,    1,    1,    1,       // Face 3: ZX       |                     |  |
                     0,    0,    0,    1,    1,    1,       // Face 4: XY     0 |                     | 1|    Z
                     0,    0,    0,    1,    1,    1  };    // Face 5: XY       |          2          |  /    |  Y
                                                            //                  |                     | /     | /        0: Dirichlet
                                                            //                  | ____________________|/      |/____X    1: Neumann
    double par_dirichlet[]={0.0000,0.0000,0.0000,
                            0.0000,0.0000,0.0000,
                            0.0000,0.0000,0.0000 };         // Grad u average
    double par_neumann[]  ={0.00};

    bc->par_periodic  = (double*) malloc( 9*sizeof(double)); for (int i=0; i< 9; i++) bc->par_periodic[i]  = par_periodic[i];
    bc->bc_type       = (int*)    malloc(36*sizeof(int));    for (int i=0; i<36; i++) bc->bc_type[i]       = bc_type[i];
    bc->par_dirichlet = (double*) malloc( 9*sizeof(double)); for (int i=0; i< 9; i++) bc->par_dirichlet[i] = par_dirichlet[i];
    bc->par_neumann   = (double*) malloc( 1*sizeof(double)); for (int i=0; i< 1; i++) bc->par_neumann[i]   = par_neumann[i];
}

void neumann(double *neumann_, double par_neumann[], int face_id, int bc_order, int idof, double xq[])
{
    int ndof=3;

    //                                          0                          |                          1                               // bc_order
    //                -------------------------------------------------------------------------------------------------------------------------
    //                        0        ,        1        ,        2        |        0        ,        1        ,        2             // idof
    //                -------------------------------------------------------------------------------------------------------------------------
    double neumann[]={        0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 0
                              0        ,  par_neumann[0] ,  par_neumann[0] ,        0        ,        0        ,        0        ,    // face 1
                              0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 2
                              0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 3
                              0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 4
                              0        ,        0        ,        0        ,        0        ,        0        ,        0         };  // face 5

    *neumann_=neumann[(2*ndof)*face_id+ndof*bc_order+idof];
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             INITIAL GUESS          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void initguess(double *U_hist, int nselfIdx, void *snesassmbl_)
{
    SNESAssmbl *snesassmbl = (SNESAssmbl*)snesassmbl_;

    for (int i=0; i<nselfIdx; i++) { U_hist[i]=U_hist[nselfIdx+i]; }

    double *U_;
    VecGetArray(snesassmbl->U,&U_);
    for (int i=0; i<nselfIdx; i++) { U_[i]=U_hist[i]; }
    VecRestoreArray(snesassmbl->U,&U_);
}

