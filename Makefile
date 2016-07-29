MODDIR = mod
SRCDIR = src
FC = gfortran
# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS = -J $(MODDIR)
# libraries needed for linking, unused in the examples
LDFLAGS =

DEBUG = -g -fbacktrace -pedantic -fwhole-file -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -fbounds-check -pg -frecursive -fcheck=all -Wall -ffpe-trap=zero,underflow,overflow
#-g turns on debugging
#-p turns on profiling

OPTIM = -O3 -fopenmp -funroll-loops

EXE = laboetie


OBJS = 	$(SRCDIR)/module_precision_kinds.f90 \
	$(SRCDIR)/module_mathematica.f90 \
	$(SRCDIR)/module_input.f90 \
	$(SRCDIR)/module_constants.f90 \
	$(SRCDIR)/module_lbmodel.f90 \
	$(SRCDIR)/module_system.f90 \
	$(SRCDIR)/module_geometry.f90 \
	$(SRCDIR)/module_moment_propagation.f90 \
	$(SRCDIR)/module_myallocations.f90 \
	$(SRCDIR)/module_collision.f90 \
	$(SRCDIR)/module_time.f90 \
	$(SRCDIR)/backup_phi_c_plus_c_minus.f90 \
	$(SRCDIR)/advect.f90 \
	$(SRCDIR)/charges_init.f90 \
	$(SRCDIR)/charge_test.f90 \
	$(SRCDIR)/check_charge_distribution_equilibrium.f90 \
	$(SRCDIR)/drop_tracers.f90 \
	$(SRCDIR)/electrostatic_pot.f90 \
	$(SRCDIR)/equilibration.f90\
	$(SRCDIR)/module_io.f90 \
	$(SRCDIR)/init_simu.f90 \
	$(SRCDIR)/just_eq_smolu.f90 \
	$(SRCDIR)/main.f90 \
	$(SRCDIR)/poisson_nernst_planck.f90 \
	$(SRCDIR)/smolu.f90 \
	$(SRCDIR)/sor.f90 \
	$(SRCDIR)/supercell_definition.f90 \
	$(SRCDIR)/velocity_profiles.f90

# symbol '@' in front of a line makes it silent. Otherwise it is printed in terminal when called

 all: $(OBJS)
	 @mkdir -p obj mod
	 $(FC) $(FCFLAGS) $(OPTIM) -o $(EXE) $(OBJS) $(LDFLAGS)

 debug: $(OBJS)
	 @mkdir -p obj mod
	 $(FC) $(FCFLAGS) $(DEBUG) -o $(EXE) $(OBJS) $(LDFLAGS)

 clean:
	rm -vf gmon.out $(EXE) $(MODDIR)/*
