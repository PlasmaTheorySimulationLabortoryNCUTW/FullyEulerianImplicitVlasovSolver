# FIVES
## Code Structure
The Fully Eulerian Implicit Vlasov Solver (FIVES) is a Fortran-90 Code, parallelized with MPI-OPENMP hybrid straegy, Post-processing with python and stored data with hdf5 format. With coupling predictor-corrector method in adam form, the implicit vlasov code is implemented. Coupling with Weighted Essential Non-Oscillatory (WENO) Scheme, Central difference and fast Cubic-Spline scheme, the numerical stability can be obtained. Besides, the usage of Ampere's law instead of Poisson equation allows the parallel efficiency can be boosted up.
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
The data is stored in the directory './data' with .h5 format, the data is separated into two main branches, one is fluid data, another is full data. In fluid data, the fluid quantity, field data and reduced distribution function(from NDNV to 1D1V) is stored and in full data, the whole distribution function and field data is stored.
Reading h5 file with python, the h5py is used and user can directly read it. There is an simple script in the analyzing directory.

