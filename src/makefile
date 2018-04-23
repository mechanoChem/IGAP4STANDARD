COMPILER      = ${PCC}
FLAGS         = -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -O3 -march=native -mtune=native -std=c99

include ${PETSC_DIR}/lib/petsc/conf/variables

INC           = ${IGAP4_DIR}/include
LIB           = ${IGAP4_DIR}/lib

LIBNAME       = igap4

OBJ           = $(patsubst %.c,%.o,$(wildcard *.c))

main: ${OBJ}
	${COMPILER} ${FLAGS} -shared -Wl,-soname,lib${LIBNAME}.so.1 -o lib${LIBNAME}.so.1.0 ${OBJ} -lc
	mv lib${LIBNAME}.so.1.0 ${LIB}
	ln -sf ${LIB}/lib${LIBNAME}.so.1.0 ${LIB}/lib${LIBNAME}.so.1
	ln -sf ${LIB}/lib${LIBNAME}.so.1.0 ${LIB}/lib${LIBNAME}.so
	rm *.o

%.o: %.c
	${COMPILER} ${FLAGS} -I${INC} -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include -o $@ -c $<

