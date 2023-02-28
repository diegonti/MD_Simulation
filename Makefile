# Fortran compiler
FC=gfortran

# Code location
VPATH=src/

# Fortran- related flags
F_FLAGS=-std=f2008

# Debug related flags:
COMP_D_FLAGS=-g3 -Og -fbacktrace -Wall -Wextra -pedantic -fcheck=all,no-array-temps -ffpe-trap=invalid,zero,overflow,underflow,denormal -ffpe-summary=all -Wconversion-extra 

# Release related flags:
COMP_R_FLAGS=#-O3 -funroll-loops -ftree-vectorize -finline-functions -flto=2 -fwhole-program -ftree-vectorizer-verbose=7

# ~ LINKING ~
all: MDEMI.x
MDEMI.x: main.o
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) $^ -o $@


# ~ COMPILING ~
main.o: main.f90
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) -c $^


# Defined recipies:
.PHONY: clean
clean:
	rm *.o

.PHONY: run_serial
run_serial:
	./MDEMI.x

