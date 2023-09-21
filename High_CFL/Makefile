IFC = mpiifort
IFC2 = h5pfc

FFLAGS = -qopenmp -Ofast -g -m64 -Wall -fexceptions -fstack-protector-strong -grecord-gcc-switches -mtune=generic -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -ldl -lm 

INIT = init.f90
FINAL = Vlasov_main.f90


all:
	$(IFC2) $(FFLAGS) $(INIT) -o init.out
	$(IFC2) $(FFLAGS) $(FINAL) -o Vlasov.out
clean:
	rm file* *.out *.mod *.h5 *.o div* time.txt 
