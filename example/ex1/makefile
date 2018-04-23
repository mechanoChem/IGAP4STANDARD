APPNAME       = getime

COMPILER      = ${PCC}
FLAGS         = -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -O3 -march=native -mtune=native -std=c99

include ${PETSC_DIR}/lib/petsc/conf/variables

INC           = ${IGAP4_DIR}/include
LIB           = ${IGAP4_DIR}/lib

OBJ           = $(patsubst %.c,%.o,$(wildcard *.c))

main: ${OBJ}
	${COMPILER} ${FLAGS} -L${LIB} -Wl,-rpath=${LIB} -o main ${OBJ} ${PETSC_SNES_LIB} -ligap4 -lnurbs -l${APPNAME}
	rm *.o

%.o: %.c
	${COMPILER} ${FLAGS} -I${INC} -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include -o $@ -c $<

