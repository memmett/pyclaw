PYCLAWMAKE = $(PYCLAW)/Makefile.common

# This makefile uses the Fortran compiler and flags from PETSc. This
# means that PETSc needs to be configured with Fortran support which
# would otherwise not be necessary. However, using arbitrary environment
# variables and hoping it all works out is way more fragile than I want
# any part of.

SHARPCLAW = $(PYCLAW)/src/pyclaw/sharpclaw
ONE_D_SHARPCLAW_SOURCES = $(SHARPCLAW)/workspace.f90 $(SHARPCLAW)/flux1.f90 $(SHARPCLAW)/weno.f90 $(SHARPCLAW)/reconstruct.f90 $(SHARPCLAW)/tfluct.f90 $(SHARPCLAW)/evec.f90

sharpclaw1.so:
	gfortran -O3 -fPIC -shared -o sharpclaw1.so -std=f2008 $(ONE_D_SHARPCLAW_SOURCES)

# FCOMPILE_SINGLE = ${FC} -c ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS}

# RP_SOURCE =  $(RIEMANN)/src/rp1_acoustics.f90

# F2PY = f2py
# PYTHON = python
# MODULE = assembly

# all: 
# 	make classic1.so
# 	make sharpclaw1.so

# ${MODULE}.so: ${MODULE}.pyx ${MODULE}impl.h ${MODULE}impl.c petscdadef.h
# 	CC='${PCC}' CFLAGS='${PCC_FLAGS}' F90='${FC_LINKER}' FFLAGS='${FFLAGS}' LDFLAGS='${FC_LINKER_FLAGS}' ${PYTHON} setup.py -q build_ext --inplace
# 	${RM} -r build ${MODULE}_wrap.c

# include $(PYCLAWMAKE)

# #include ${PETSC_DIR}/conf/variables

