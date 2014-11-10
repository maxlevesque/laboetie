MODDIR = mod
SRCDIR = src
FC = gfortran
# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS = -J $(MODDIR)
# libraries needed for linking, unused in the examples
LDFLAGS = -lfftw3

DEBUG = -Og -g -Wall -Wextra -fimplicit-none -fbacktrace -std=f2008 -pedantic -fwhole-file -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -fbounds-check -fcheck=all -fcheck-array-temporaries -Warray-temporaries -Wconversion -pg -Wunused-parameter -Wimplicit-interface -frecursive
#-g turns on debugging
#-p turns on profiling

OPTIM = -O3 -march=native -fopenmp #-ffast-math -funroll-loops
# -fopenmp for OPENMP support

EXE = laboetie


OBJS = 	$(SRCDIR)/mod_precision_kinds.f90 \
		$(SRCDIR)/mod_input.f90 \
		$(SRCDIR)/mod_constants.f90 \
		$(SRCDIR)/mod_lbmodel.f90 \
		$(SRCDIR)/mod_system.f90 \
		$(SRCDIR)/mod_geometry.f90 \
		$(SRCDIR)/mod_moment_propagation.f90 \
		$(SRCDIR)/mod_myallocations.f90 \
		$(SRCDIR)/mod_populations.f90 \
		$(SRCDIR)/mod_supercell.f90 \
	    	$(SRCDIR)/backup_phi_c_plus_c_minus.f90 \
		$(SRCDIR)/advect.f90 \
		$(SRCDIR)/charges_init.f90 \
		$(SRCDIR)/charge_test.f90 \
		$(SRCDIR)/check_charge_distribution_equilibrium.f90 \
		$(SRCDIR)/close_simu.f90 \
		$(SRCDIR)/comp_j.f90 \
		$(SRCDIR)/comp_rho.f90 \
		$(SRCDIR)/drop_tracers.f90 \
		$(SRCDIR)/electrostatic_pot.f90 \
		$(SRCDIR)/equilibration_with_constraints.f90 \
		$(SRCDIR)/equilibration_without_constraint.f90 \
		$(SRCDIR)/init_simu.f90 \
		$(SRCDIR)/just_eq_smolu.f90 \
		$(SRCDIR)/main.f90 \
		$(SRCDIR)/module_io.f90 \
		$(SRCDIR)/poisson_nernst_planck.f90 \
		$(SRCDIR)/propagation.f90 \
		$(SRCDIR)/smolu.f90 \
		$(SRCDIR)/sor.f90 \
		$(SRCDIR)/supercell_definition.f90 \
		$(SRCDIR)/velocity_profiles.f90

# symbol '@' in front of a line makes it silent. Otherwise it is printed in terminal when called

 all: $(OBJS)
	 $(FC) $(FCFLAGS) $(OPTIM) -o $(EXE) $(OBJS) $(LDFLAGS)

 debug: $(OBJS)
	 $(FC) $(FCFLAGS) $(DEBUG) -o $(EXE) $(OBJS) $(LDFLAGS)

 clean:
	rm -vf gmon.out $(EXE) $(MODDIR)/* 
