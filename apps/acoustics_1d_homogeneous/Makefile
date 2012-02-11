# This Makefile...


F90=gfortran
F90FLAGS=-fPIC -pg -g -Jbuild -Ibuild
#F90FLAGS=-fPIC -O3 -Jbuild -Ibuild

SHARPCLAW = $(PYCLAW)/src/pyclaw/sharpclaw
vpath %.f90 . $(SHARPCLAW) $(RIEMANN)/src


ONE_D_SHARPCLAW_OBJS = \
	workspace.o weno.o reconstruct.o \
	rp1_acoustics.o tfluct.o evec.o \
	flux1.o 


OBJS = $(addprefix build/, $(ONE_D_SHARPCLAW_OBJS))


libsharpclaw_.so: $(OBJS)
	$(F90) -shared -o libsharpclaw_.so $(OBJS)

build/%.o: %.f90
	@mkdir -p build
	$(F90) $(F90FLAGS) -c $(OUTPUT_OPTION) $<

clean:
	rm -rf libsharpclaw_.so build
