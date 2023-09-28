# FIVES
## Code Structure
The Fully Eulerian Implicit Vlasov Solver (FIVES) is a Fortran-90 Code, parallelized with MPI-OPENMP hybrid straegy, Post-processing with python and stored data with hdf5 format. With coupling predictor-corrector method in adam form, the implicit vlasov code is implemented.
## Compliling the code
The code is separated into two parts - initialization and main execution.
First, copy the initial file.
```
cp CaseX.f90 init.f90
```
Second, the file is compliled via Makefile include in the code.
```
make
```
### Run the code
the code structure is separated into two part, initial file and main code, which is compiled in previous section and generate two out files - init.out and Vlasov.out.
The init.h5 is generated via
```
./init.out
```
with the initialization, the main code can be executed via
```
mpirun -np X ./Vlasov.out
```
which the core number is assigned in init_para.txt file. the openmp command can be also included and also be assigned in init_para.txt.
### Data Post-processing
