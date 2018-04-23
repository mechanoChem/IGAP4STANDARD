# IGAP4STANDARD
An Isogeometric Analysis Program for partial differential equations that involve fourth-order spatial derivatives in standard code architecture.

## Overview
IGAP4STANDARD is written in C. It requires PETSc [https://www.mcs.anl.gov/petsc/] for linear/nonlinear solvers. 
IGAP4STANDARD can be used to numerically solve (initial) boundary value problems of partial differential equations that involve fourth-order spatial derivatives; examples include (transient) gradient elasticity and the Cahn-Hilliard equation. 

## Installation on Linux

1) Install PETSc 3.7.x [https://www.mcs.anl.gov/petsc/]. Set environmental variables PETSC_DIR and PETSC_ARCH as required by PETSc.

2) Set environmental variable, IGAP4_DIR, on the command line:
> export IGAP4_DIR=/path/to/your/igap4standard/dir

3) Compile the source and make shared objects:
> cd ${IGAP4_DIR}; make

## Tests

> cd ${IGAP4_DIR}/example/ex1; make; mpiexec -n 8 ./main

## Contributors

This code has been developed at the Computational Physics Group at the University of Michigan [http://www.umich.edu/~compphys/index.html].

- Koki Sagiyama (Lead Developer)
- Krishna Garikipati

## Acknowledgements

The development of this code has been supported by the following:

- Dept of Energy (DoE and labs) : Software Center for Predictive Theory and Modeling (DE-SC0008637) 
- National Science Foundation : Integrated Computational Framework for Designing Dynamically Controlled Alloy-Oxide Heterostructures (DMR1436154)

## License
BSD 2-Clause License
