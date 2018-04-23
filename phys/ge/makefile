COMPILER      = ${PCC}
FLAGS         = -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -O3 -march=native -mtune=native -std=c99

INC           = ${IGAP4_DIR}/include
LIB           = ${IGAP4_DIR}/lib

MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
APPLIBNAME    := $(notdir $(patsubst %/,%,$(dir $(MAKEFILE_PATH))))

include ${PETSC_DIR}/lib/petsc/conf/variables

OBJ           = $(patsubst %.c,%.o,$(wildcard *.c))

main: ${OBJ}
	${COMPILER} -shared -Wl,-soname,lib${APPLIBNAME}.so.1 -o lib${APPLIBNAME}.so.1.0 ${OBJ} -lc
	mv lib${APPLIBNAME}.so.1.0 ${LIB}
	ln -sf ${LIB}/lib${APPLIBNAME}.so.1.0 ${LIB}/lib${APPLIBNAME}.so.1
	ln -sf ${LIB}/lib${APPLIBNAME}.so.1.0 ${LIB}/lib${APPLIBNAME}.so
	rm *.o

%.o: %.c
	${COMPILER} ${FLAGS} -I${INC} -o $@ -c $<

