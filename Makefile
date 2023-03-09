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
MDEMI.x:  pbc.o potentials_module.o simulation.o writers_mod.o  readers_mod.o testing.o initialization.o integrators.o main.o
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) $^ -o $@


# ~ COMPILING ~
main.o: main.f90 # initialization.o
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

potentials_module.o: potentials_module.f90
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) -c $^

simulation.o: simulation.f90
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) -c $^

integrators.o: integrators.f90
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) -c $^


.PHONY: run_serial
run_serial:
	./MDEMI.x parameters.nml


.PHONY: postprocess stats plots
postprocess: stats plots
stats:
	python3 scripts/stats.py -ip $(input) $(args) #-op output_folder -s start -f finish 

plots:
	python3 scripts/visualization.py -ip $(input) $(args) #-op output_folder -s start -f finish 


# Defined recipies:
.PHONY: clean
clean:
	rm -f *.o *.mod
	rm -f ./src/*.o ./src/*.mod

.PHONY: help
help:
	@sed -n "s/^##//p" Makefile



