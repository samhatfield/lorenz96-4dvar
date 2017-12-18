# Default compiler
FC = gfortran

default: all

all: test

# Main target: main executable
test: test.o ../utils.o ../params.o ../lorenz96.o ../assim.o
	$(FC) -o $@ $^

# Dependencies
test.o: ../utils.o ../params.o ../lorenz96.o ../assim.o

# Build rules
%.o: %.f90
	$(FC) -c $< -o $(basename $<).o


.PHONY: clean
clean:
	rm -f *.o *.mod main
