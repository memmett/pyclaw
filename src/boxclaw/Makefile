BOXLIB_HOME ?= /home/memmett/projects/BoxLib

MPI  := t
COMP := gfortran

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

CC       := mpicc
FC       := mpif90
F90      := mpif90
F90FLAGS += -fPIC
CFLAGS   += -fPIC

# load the BoxLib and PyBoxLib packages
include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

include $(BOXLIB_HOME)/Src/Python/GPackage.mak
VPATH_LOCATIONS += src

all: pyfboxlib.so

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
include $(BOXLIB_HOME)/Src/Python/GMakerules.mak
