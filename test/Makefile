# Default compiler
FC = gfortran

# Compiler options
OPTS = -O2

# Modules to include
MODS = -I../rpe/modules

# Libraries to link 
LIBS = -lrpe -L../rpe/lib


default: all

all: test

# Main target: main executable
test: test.o ../utils.o ../params.o ../lorenz96.o ../assim.o
	$(FC) $(OPTS) -o $@ $^ $(LIBS)

# Dependencies
test.o: ../utils.o ../params.o ../lorenz96.o ../assim.o
../utils.o: ../params.o

# Build rules
%.o: %.f90
	$(FC) $(OPTS) -c $< -o $(basename $<).o $(MODS)


.PHONY: clean
clean:
	rm -f *.o *.mod main
