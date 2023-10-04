
COMPILE_OPTIONS := -c -i8 -I"${MKLROOT}/include"
# Project name
NAME := test_mkl_rng

LINK_OPTIONS = -Wl,--start-group 
LINK_OPTIONS += ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a
LINK_OPTIONS += ${MKLROOT}/lib/intel64/libmkl_sequential.a
LINK_OPTIONS += ${MKLROOT}/lib/intel64/libmkl_core.a 
LINK_OPTIONS += -Wl,--end-group -lpthread -lm -ldl
# Configuration settings
#FC := gfortran
FC := ifort
RM := rm -f
GD := ./gen_dep.awk

OBJS = mod_random_number_generator.o mod_random_number_generator_test.o

.PHONY: all clean

.SUFFIXES: .f90 .o
# general rule
.f90.o:
	$(FC) $(COMPILE_OPTIONS) $<

all: mod_random_number_generator_test


mod_random_number_generator_test: $(OBJS)
	$(FC) mod_random_number_generator.o \
	mod_random_number_generator_test.o \
	$(LINK_OPTIONS) -o mod_random_number_generator_test


# modules
$(OBJS): %.o
	$(FC) $(COMPILE_OPTIONS) $@ $<

mod_random_number_generator_test.o: mod_random_number_generator_test.f90
	$(FC) $(COMPILE_OPTIONS) mod_random_number_generator_test.f90
.PHONY:

clean:
	$(RM) *.o *.mod mod_random_number_generator_test.out
