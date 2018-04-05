# Default compiler
FC = gfortran

# Compiler options
OPTS = -O2

# Main target: main executable
main: main.o params.o lorenz96.o utils.o io.o assim.o cg_plus.f
	$(FC) $(OPTS) -o $@ $^

# Dependencies
main.o: params.o lorenz96.o utils.o io.o assim.o
lorenz63.o: params.o
utils.o: params.o
io.o: params.o
assim.o: lorenz96.o params.o

# Build rules
%.o: %.f90
	$(FC) $(OPTS) -c $< -o $(basename $<).o

.PHONY: clean
clean:
	rm -f *.o *.mod main
