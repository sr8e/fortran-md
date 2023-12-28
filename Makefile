BIN=md_3atom md_bulk_Al md_nanowire
BIN_ALONE=lattice
MOD=utils.o potentials.o md_struct.o io.o

all: $(BIN);

$(BIN):%:%.f95 $(MOD)
	gfortran $^ -o $@

$(BIN_ALONE):%:%.f95
	gfortran $< -o $@

%.o: %.f95
	gfortran -c $< -o $@

clean:
	$(RM) *.o *.mod *.txt $(BIN) $(BIN_ALONE)

.PHONY: null clean
