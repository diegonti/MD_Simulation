##
# Fortran compiler
FC=mpifort

# Code location
VPATH=src/

# Fortran- related flags
F_FLAGS=-std=f2018

# Debug related flags:
COMP_D_FLAGS=-g3 -Og -fbacktrace -Wall -Wextra -pedantic -fcheck=all,no-array-temps -ffpe-trap=invalid,zero,overflow,underflow,denormal -ffpe-summary=all -Wconversion-extra 

# Release related flags:
COMP_R_FLAGS=#-O3 -funroll-loops -ftree-vectorize -finline-functions -flto=2 -fwhole-program -ftree-vectorizer-verbose=7

# ~ LINKING ~
## make: Compiles the program.
all: MDEMI.x
MDEMI.x:  pbc.o potentials_module.o simulation.o writers_mod.o  readers_mod.o testing.o initialization.o integrators.o main.o
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) $^ -o $@
##

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


## run_parallel: Runs in parallel the simulaion executable with the parameters file.
.PHONY: run_parallel
run_parallel:
	mpirun -np $(nproc) MDEMI.x parameters.nml
##

## run_serial: Runs sequantially the simulaion executable with the parameters file.
.PHONY: run_serial
run_serial:
	./MDEMI.x parameters.nml
##

## postprocess: runs the statistics and visualization Python scripts. 
##	use with input="name_logfile.log" and optional flags with args="script flags".
## stats: runs only the statistics of the output file. Use the same arguments as postprocess.
## plots: performs only the visualization of the output file. Use the same arguments as postprocess.

.PHONY: postprocess stats plots
postprocess: stats plots
stats:
	python3 scripts/stats.py -ip $(input) $(args) #-op output_folder -s start -f finish 

plots:
	python3 scripts/visualization.py -ip $(input) $(args) #-op output_folder -s start -f finish 
##

# Defined recipies:
## clean: Cleans compilation files (*.o *.mod).
.PHONY: clean
clean:
	rm -f *.o *.mod
	rm -f ./src/*.o ./src/*.mod

## deep_clean: Cleans also output files (*.log *.xyz).
.PHONY: deep_clean
deep_clean:
	rm -f *.o *.mod
	rm -f ./src/*.o ./src/*.mod
	rm -f *.log *.xyz
##

## help: Gives information of each make command.
.PHONY: help
help:
	@sed -n "s/^##//p" Makefile
##


