SRC=md.f95
BIN=md_exec
MOD=utils.o potentials.o md_struct.o io.o

compile: $(MOD)
	gfortran $^ $(SRC) -o $(BIN)

lattice: lattice.f95
	gfortran $< -o $@

%.o: %.f95
	gfortran -c $< -o $@

clean:
	$(RM) *.o *.mod $(BIN)
