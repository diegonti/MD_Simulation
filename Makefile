# Fortran compiler
FC=gfortran

# Code location
VPATH=src/

# Fortran- related flags
F_FLAGS=-std=f2018

# Debug related flags:
COMP_D_FLAGS=-g3 -Og -fbacktrace -Wall -Wextra -pedantic -fcheck=all,no-array-temps -ffpe-trap=invalid,zero,overflow,underflow,denormal -ffpe-summary=all -Wconversion-extra 

# Release related flags:
COMP_R_FLAGS=#-O3 -funroll-loops -ftree-vectorize -finline-functions -flto=2 -fwhole-program -ftree-vectorizer-verbose=7

# ~ LINKING ~
all: MDEMI.x
MDEMI.x:  simulation.o integrators.o potentials_module.o pbc.o readers_mod.o writers_mod.o testing.o initialization.o main.o
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) $^ -o $@


# ~ COMPILING ~
main.o: main.f90 
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) -c $^

initialization.o: initialization.f90
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) -c $^

testing.o: testing.f90
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) -c $^

readers_mod.o: readers_mod.f90
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) -c $^

writers_mod.o: writers_mod.f90
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) -c $^

pbc.o: pbc.f90
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) -c $^

potentials_module.o: potentials_module.f90 pbc.o
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) -c $^

integrators.o: integrators.f90 pbc.o potentials_module.o
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) -c $^

simulation.o: simulation.f90 pbc.o
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) -c $^


# Defined recipies:
.PHONY: clean
clean:
	rm -f *.o *.mod
	rm -f ./src/*.o ./src/*.mod

.PHONY: run_serial
run_serial:
	./MDEMI.x

