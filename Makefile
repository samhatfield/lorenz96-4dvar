# Default compiler
FC = gfortran

# Compiler options
OPTS = -O2

# Modules to include
MODS = -Irpe/modules

# Libraries to link 
LIBS = -lrpe -Lrpe/lib

# Main target: main executable
main: main.o params.o lorenz96.o utils.o io.o assim.o minimisation.o cg_plus.f90
	$(FC) $(OPTS) -o $@ $^ $(LIBS)

# Dependencies
main.o: params.o lorenz96.o utils.o io.o assim.o minimisation.o
lorenz96.o: params.o
utils.o: params.o
io.o: params.o
assim.o: lorenz96.o params.o
minimisation.o: params.o assim.o

# Build rules
%.o: %.f90
	$(FC) $(OPTS) -c $< -o $(basename $<).o $(MODS)

.PHONY: clean
clean:
	rm -f *.o *.mod main
