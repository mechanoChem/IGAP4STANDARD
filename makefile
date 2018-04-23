main:
	mkdir ${IGAP4_DIR}/lib
	cd ${IGAP4_DIR}/src; make
	cd ${IGAP4_DIR}/mesh/nurbs; make
	cd ${IGAP4_DIR}/phys/ge; make
	cd ${IGAP4_DIR}/phys/getime; make

