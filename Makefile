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
MDEMI.x: readers_mod.o testing.o initialization.o main.o
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


# Defined recipies:
.PHONY: clean
clean:
	rm *.o *.mod

.PHONY: run_serial
run_serial:
	./MDEMI.x

