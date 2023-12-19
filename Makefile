SRC=md.f95
BIN=md_exec
MOD=utils.o potentials.o md_struct.o

compile: $(MOD)
	gfortran $^ $(SRC) -o $(BIN)

%.o: %.f95
	gfortran -c $< -o $@

clean:
	$(RM) *.o *.mod $(BIN)
