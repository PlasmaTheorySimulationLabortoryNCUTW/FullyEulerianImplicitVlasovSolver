module bcstrecv
     use wtime
     use comm_global
     use MPI3D
     use Vlasov
     use GRID6D
     use coef
     use Vlasov2fluid
   contains
!===============================================================================
!===============================================================================   
   subroutine initial
     implicit double precision (A-H,O-Z)
     character*10 :: dates, times, zone
     integer, dimension(:) :: values(8)
     call Date_and_Time(dates,times,zone,values)
     call const
!
!  allocate the memory space for global array
!
     call setVlasov_global
!
!  read the data from the init.bin
!
     if (MyID == 0) call get_init_data_hdf5
     call MPI_bcast(t ,1,MPI_real8,0,MPI_COMM_WORLD,IERR)
     call MPI_bcast(it,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
     call MPI_bcast(dt,1,MPI_real8,0,MPI_COMM_WORLD,IERR)
     call MPI_bcast(uex,nuex,MPI_real8,0,MPI_COMM_WORLD,IERR)
     call MPI_bcast(uey,nuey,MPI_real8,0,MPI_COMM_WORLD,IERR)
     call MPI_bcast(uez,nuez,MPI_real8,0,MPI_COMM_WORLD,IERR)
     call MPI_bcast(uix,nuix,MPI_real8,0,MPI_COMM_WORLD,IERR)
     call MPI_bcast(uiy,nuiy,MPI_real8,0,MPI_COMM_WORLD,IERR)
     call MPI_bcast(uiz,nuiz,MPI_real8,0,MPI_COMM_WORLD,IERR)
!
     call MPI_bcast(Istop,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
     if (Istop == 999) then
     	 call MPI_FINALIZE(IERR)
     	 stop
     endif
!
!  allocate the memory space for local array
!
     call setVlasov
!-------------------------------------------------------------------------------
!  read and brocast the 6-D electron distribution data to the other node
!
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     if (MyID == 0) then
       allocate (fe_global(ncxyz*nuexyz))
     else
       allocate (fe_global(1))
     endif
!
     if (MyID == 0) call read_init_array("fe",ncxyz*nuexyz,fe_global)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call data_decomposition_6D(fe_global,f(ncfe),nuex,nuey,nuez)
     deallocate (fe_global)    
!-------------------------------------------------------------------------------
!  read and brocast the 6-D ion distribution data to the other node
!
     if (MyID == 0) then
       allocate (fi_global(ncxyz*nuixyz))
     else
       allocate (fi_global(1))
     endif
!
     if (MyID == 0) call read_init_array("fi",ncxyz*nuixyz,fi_global)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call data_decomposition_6D(fi_global,f(ncfi),nuix,nuiy,nuiz)
     deallocate (fi_global)
!-------------------------------------------------------------------------------
!  read and brocast the Bx, By, Bz, Ex, Ey, Ez data to the other node
!
     if (MyID == 0) then
       allocate (A_global(ncxyz))
     else
       allocate (A_global(1))
     endif
!
     if (MyID == 0) call read_init_array("Bx",ncxyz,A_global)   ! read Bx
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call data_decomposition(A_global,f(ncBx))
!
     if (MyID == 0) call read_init_array("By",ncxyz,A_global)   ! read By
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call data_decomposition(A_global,f(ncBy))
!
     if (MyID == 0) call read_init_array("Bz",ncxyz,A_global)   ! read Bz
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call data_decomposition(A_global,f(ncBz))
!
     if (MyID == 0) call read_init_array("Ex",ncxyz,A_global)   ! read Ex
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call data_decomposition(A_global,f(ncEx))
!
     if (MyID == 0) call read_init_array("Ey",ncxyz,A_global)   ! read Ey
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call data_decomposition(A_global,f(ncEy))
!
     if (MyID == 0) call read_init_array("Ez",ncxyz,A_global)   ! read Ez
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call data_decomposition(A_global,f(ncEz))

!-------------------------------------------------------------------------------
!  read and brocast the etaB, etafe, etafi data to the other node
!

     if (MyID == 0) call read_eta_array("etaBB",ncxyz,A_global)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call data_decomposition(A_global,etaB)
! 
     if (MyID == 0) call read_eta_array("etafe",ncxyz,A_global)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call data_decomposition(A_global,etafe)
!
     if (MyID == 0) call read_eta_array("etafi",ncxyz,A_global)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call data_decomposition(A_global,etafi)

     deallocate (A_global)
!-------------------------------------------------------------------------------
     call alloc_3D3V
!
     call data_decomposition_1D(x,x_global,ncx_mpi,ix_start)
     call data_decomposition_1D(y,y_global,ncy_mpi,iy_start)
     call data_decomposition_1D(z,z_global,ncz_mpi,iz_start)
!
     dx = x(2)-x(1)
!
     if (ncy_mpi == 1) then
       dy = 0.d0
     else
       dy = y(2)-y(1)
     endif
!
     if (ncz_mpi == 1) then
       dz = 0.d0
     else
       dz = z(2)-z(1)
     endif
!-------------------------------------------------------------------------------
!  check the CFL criterion
!
!     if (C*dt > 0.5d0*dx) then
!       if (MyID == 0) then
!         open (1,file='error.txt',status='unknown')
!         write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
!                      values(5), ':', values(6)
!         write (1,*) 'The seting in ini_para.txt have errors !!'
!         write (1,*) 'dt is too large, C*dt > dx/2'
!         write (1,*) 'Please check the setting of dt !!'
!         write (1,*) 'dt =', dt
!         write (1,*) 'Program STOP at initial !!'
!         close (1)
!         write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
!                      values(5), ':', values(6)
!         write (6,*) 'The seting in ini_para.txt have errors !!'
!         write (6,*) 'dt is too large, C*dt > dx/2'
!         write (6,*) 'Please check the setting of dt !!'
!         write (6,*) 'dt =', dt
!         write (6,*) 'Program STOP at initial !!'
!       endif
!       call MPI_FINALIZE(IERR)
!       STOP
!     endif
! 11  format (1x,i4,4(a1,i2.2))
!
115  continue
!-------------------------------------------------------------------------------
   end subroutine initial
!===============================================================================
!  read the data from the file 'init.bin'
!  the file 'init.bin' is written in fortran binary format
!

   subroutine read_init_array(dsetname,size,data_out)
       use hdf5
       IMPLICIT NONE
       integer,intent(in) :: size
       character(len=*):: dsetname
       CHARACTER(LEN=7), PARAMETER :: filename = "init.h5"
       INTEGER(HID_T) :: file_id
       INTEGER(HID_T) :: dset_id
       INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
       INTEGER     ::   error
       real*8,intent(inout),dimension(size) :: data_out
       data_dims(1) = size
       CALL h5open_f(error)
       CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
       CALL h5dopen_f(file_id, dsetname, dset_id, error)
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, data_dims, error)
       CALL h5dclose_f(dset_id, error)
       CALL h5fclose_f(file_id, error)
       CALL h5close_f(error)
   end subroutine read_init_array
   subroutine read_eta_array(dsetname,size,data_out)
       use hdf5
       IMPLICIT NONE
       integer,intent(in) :: size
       CHARACTER(LEN=6), PARAMETER :: filename = "eta.h5"
       character(len=*):: dsetname
       INTEGER(HID_T) :: file_id
       INTEGER(HID_T) :: dset_id
       INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
       INTEGER     ::   error
       real*8,intent(inout),dimension(size) :: data_out
       data_dims(1) = size
       CALL h5open_f(error)
       CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
       CALL h5dopen_f(file_id, dsetname, dset_id, error)
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, data_dims, error)
       CALL h5dclose_f(dset_id, error)
       CALL h5fclose_f(file_id, error)
       CALL h5close_f(error)
   end subroutine read_eta_array

   subroutine get_init_data_hdf5
       use hdf5
       IMPLICIT NONE
       CHARACTER(LEN=7), PARAMETER :: filename = "axis.h5"
       INTEGER(HID_T) :: file_id
       INTEGER(HID_T) :: dset_id
       INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
       INTEGER     ::   error
       real*8, allocatable, dimension(:) :: data_out
       CALL h5open_f(error)
       CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
       CALL h5dopen_f(file_id, "x", dset_id, error)
       data_dims(1) = ncx
       allocate(data_out(data_dims(1)))
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, data_dims, error)
       x_global = data_out
       deallocate(data_out)
       CALL h5dclose_f(dset_id, error)
       CALL h5dopen_f(file_id, "y", dset_id, error)
       data_dims(1) = ncy
       allocate(data_out(data_dims(1)))
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, data_dims, error)
       y_global = data_out
       deallocate(data_out)
       CALL h5dclose_f(dset_id, error)
       CALL h5dopen_f(file_id, "z", dset_id, error)
       data_dims(1) = ncz
       allocate(data_out(data_dims(1)))
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, data_dims, error)
       z_global = data_out
       deallocate(data_out)
       CALL h5dclose_f(dset_id, error)
       CALL h5dopen_f(file_id, "uex", dset_id, error)
       data_dims(1) = nuex
       allocate(data_out(data_dims(1)))
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, data_dims, error)
       uex = data_out
       deallocate(data_out)
       CALL h5dclose_f(dset_id, error)
       CALL h5dopen_f(file_id, "uey", dset_id, error)
       data_dims(1) = nuey
       allocate(data_out(data_dims(1)))
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, data_dims, error)
       uey = data_out
       deallocate(data_out)
       CALL h5dclose_f(dset_id, error)
       CALL h5dopen_f(file_id, "uez", dset_id, error)
       data_dims(1) = nuez
       allocate(data_out(data_dims(1)))
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, data_dims, error)
       uez = data_out
       deallocate(data_out)
       CALL h5dclose_f(dset_id, error)
       CALL h5dopen_f(file_id, "uix", dset_id, error)
       data_dims(1) = nuix
       allocate(data_out(data_dims(1)))
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, data_dims, error)
       uix = data_out
       deallocate(data_out)
       CALL h5dclose_f(dset_id, error)
       CALL h5dopen_f(file_id, "uiy", dset_id, error)
       data_dims(1) = nuiy
       allocate(data_out(data_dims(1)))
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, data_dims, error)
       uiy = data_out
       deallocate(data_out)
       CALL h5dclose_f(dset_id, error)
       CALL h5dopen_f(file_id, "uiz", dset_id, error)
       data_dims(1) = nuiz
       allocate(data_out(data_dims(1)))
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, data_dims, error)
       uiz = data_out
       deallocate(data_out)
       CALL h5dclose_f(dset_id, error)
       CALL h5fclose_f(file_id, error)
       CALL h5close_f(error)
       allocate(data_out(data_dims(1)))
       CALL h5open_f(error)
       CALL h5fopen_f ("init.h5", H5F_ACC_RDWR_F, file_id, error)
       CALL h5dopen_f(file_id, "axis", dset_id, error)
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, data_dims, error)
       CALL h5dclose_f(dset_id, error)
       CALL h5fclose_f(file_id, error)
       CALL h5close_f(error)
       ncx  = data_out(1)
       ncy  = data_out(2)
       ncz  = data_out(3)
       nuex = data_out(4)
       nuey = data_out(5)
       nuez = data_out(6)
       nuix = data_out(7)
       nuiy = data_out(8)
       nuiz = data_out(9)
       ami  = data_out(10)
       C    = data_out(11)
       t    = data_out(12)
       it   = data_out(13)
       dt   = data_out(14)
       deallocate(data_out)
   end subroutine get_init_data_hdf5

   subroutine get_init_data(Istop)
     implicit double precision (A-H,O-Z)
     character*10 :: dates, times, zone
     integer, dimension(:) :: values(8)
!
     call Date_and_Time(dates,times,zone,values)
!
!    open (1,file='init.bin',status='old',form='binary')
     open (1,file='init.bin',status='old',form='unformatted')
     read (1) mx, my, mz, muex, muey, muez, muix, muiy, muiz, ami1, C1, t, it, dt
     close (1)
!-------------------------------------------------------------------------------
!  check the dimension along x, y, z, direction
!
     if ((mx /= ncx) .or. (my /= ncy) .or. (mz /= ncz)) then
       open (9,file='error.txt',status='unknown')
       write (9,99) values(1), '-', values(2), '-', values(3), ' ', &
                    values(5), ':', values(6)
       write (9,*) 'The seting in MPI have errors !!'
       write (9,*) 'Please check ncx, ncy, ncz in init_para.txt and init.bin !!'
       write (9,*) 'ncx =', ncx, 'in init_para.txt'
       write (9,*) 'ncx =',  mx, 'in init.bin'
       write (9,*) 'ncy =', ncy, 'in init_para.txt'
       write (9,*) 'ncy =',  my, 'in init.bin'
       write (9,*) 'ncz =', ncz, 'in init_para.txt'
       write (9,*) 'ncz =',  mz, 'in init.bin'
       write (9,*) 'Program STOP at get_init_data !!'       
       close (9)
       Istop = 999
       return
     endif
!-------------------------------------------------------------------------------
!  check the dimension along uex, uey, uez, direction
!
     if ((muex /= nuex) .or. (muey /= nuey) .or. (muez /= nuez)) then
       open (9,file='error.txt',status='unknown')
       write (9,99) values(1), '-', values(2), '-', values(3), ' ', &
                    values(5), ':', values(6)
       write (9,*) 'The seting in MPI have errors !!'
       write (9,*) 'Please check nuex, nuey, nuez in init_para.txt and init.bin !!'
       write (9,*) 'nuex =', nuex, 'in init_para.txt'
       write (9,*) 'nuex =', muex, 'in init.bin'
       write (9,*) 'nuey =', nuey, 'in init_para.txt'
       write (9,*) 'nuey =', muey, 'in init.bin'
       write (9,*) 'nuez =', nuez, 'in init_para.txt'
       write (9,*) 'nuez =', muez, 'in init.bin'
       write (9,*) 'Program STOP at get_init_data !!'       
       close (9)
       Istop = 999
       return
     endif
!-------------------------------------------------------------------------------
!  check the dimension along uix, uiy, uiz, direction
!
     if ((muix /= nuix) .or. (muiy /= nuiy) .or. (muiz /= nuiz)) then
       open (9,file='error.txt',status='unknown')
       write (9,99) values(1), '-', values(2), '-', values(3), ' ', &
                    values(5), ':', values(6)
       write (9,*) 'The seting in MPI have errors !!'
       write (9,*) 'Please check nuix, nuiy, nuiz in init_para.txt and init.bin !!'
       write (9,*) 'nuex =', nuix, 'in init_para.txt'
       write (9,*) 'nuex =', muix, 'in init.bin'
       write (9,*) 'nuey =', nuiy, 'in init_para.txt'
       write (9,*) 'nuey =', muiy, 'in init.bin'
       write (9,*) 'nuez =', nuiz, 'in init_para.txt'
       write (9,*) 'nuez =', muiz, 'in init.bin'
       write (9,*) 'Program STOP at get_init_data !!'       
       close (9)
       Istop = 999
       return
     endif
!-------------------------------------------------------------------------------
!  check the time information t, it, and dt
!
     dth = 0.5d0*dt
     if (dabs(t-dfloat(it)*dt) > dth) then
       open (9,file='error.txt',status='unknown')
       write (9,99) values(1), '-', values(2), '-', values(3), ' ', &
                    values(5), ':', values(6)
       write (9,*) 'The seting in MPI have errors !!'
       write (9,*) 'Please check t, it, and dt in init.bin !!'
       write (9,*) 't  =',  t, ' in init.bin'
       write (9,*) 'it =', it, ' in init.bin'
       write (9,*) 'dt =', dt, ' in init.bin'
       write (9,*) 'Program STOP at get_init_data !!'       
       close (9)
       Istop = 999
       return
     endif
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  read the axis data from the file 'axis.bin'
!
!     open (3,file='axis.bin',status='old',form='binary')
     open (3,file='axis.bin',status='old',form='unformatted')
     read (3) mx, my, mz, muex, muey, muez, muix, muiy, muiz, C1
!
!  check the dimension along x, y, z, direction
!
     if ((mx /= ncx) .or. (my /= ncy) .or. (mz /= ncz)) then
       open (9,file='error.txt',status='unknown')
       write (9,99) values(1), '-', values(2), '-', values(3), ' ', &
                    values(5), ':', values(6)
       write (9,*) 'The seting in MPI have errors !!'
       write (9,*) 'Please check ncx, ncy, ncz in init_para.txt and axis.bin !!'
       write (9,*) 'ncx =', ncx, 'in init_para.txt'
       write (9,*) 'ncx =',  mx, 'in axis.bin'
       write (9,*) 'ncy =', ncy, 'in init_para.txt'
       write (9,*) 'ncy =',  my, 'in axis.bin'
       write (9,*) 'ncz =', ncz, 'in init_para.txt'
       write (9,*) 'ncz =',  mz, 'in axis.bin'
       write (9,*) 'Program STOP at get_init_data !!'       
       close (9)
       Istop = 999
       close (3)
       return
     endif
!-------------------------------------------------------------------------------
!  check the dimension along uex, uey, uez, direction
!
     if ((muex /= nuex) .or. (muey /= nuey) .or. (muez /= nuez)) then
       open (9,file='error.txt',status='unknown')
       write (9,99) values(1), '-', values(2), '-', values(3), ' ', &
                    values(5), ':', values(6)
       write (9,*) 'The seting in MPI have errors !!'
       write (9,*) 'Please check nuex, nuey, nuez in init_para.txt and init.bin !!'
       write (9,*) 'nuex =', nuex, 'in init_para.txt'
       write (9,*) 'nuex =', muex, 'in init.bin'
       write (9,*) 'nuey =', nuey, 'in init_para.txt'
       write (9,*) 'nuey =', muey, 'in init.bin'
       write (9,*) 'nuez =', nuez, 'in init_para.txt'
       write (9,*) 'nuez =', muez, 'in init.bin'
       write (9,*) 'Program STOP at get_init_data !!'       
       close (9)
       Istop = 999
       close (3)
       return
     endif
!-------------------------------------------------------------------------------
!  check the dimension along uix, uiy, uiz, direction
!
     if ((muix /= nuix) .or. (muiy /= nuiy) .or. (muiz /= nuiz)) then
       open (9,file='error.txt',status='unknown')
       write (9,99) values(1), '-', values(2), '-', values(3), ' ', &
                    values(5), ':', values(6)
       write (9,*) 'The seting in MPI have errors !!'
       write (9,*) 'Please check nuix, nuiy, nuiz in init_para.txt and init.bin !!'
       write (9,*) 'nuex =', nuix, 'in init_para.txt'
       write (9,*) 'nuex =', muix, 'in init.bin'
       write (9,*) 'nuey =', nuiy, 'in init_para.txt'
       write (9,*) 'nuey =', muiy, 'in init.bin'
       write (9,*) 'nuez =', nuiz, 'in init_para.txt'
       write (9,*) 'nuez =', muiz, 'in init.bin'
       write (9,*) 'Program STOP at get_init_data !!'       
       close (9)
       Istop = 999
       close (3)
       return
     endif
!-------------------------------------------------------------------------------
     if (nuex == 1) then
       open (1,file='error.txt',status='unknown')
       write (1,99) values(1), '-', values(2), '-', values(3), ' ', &
                    values(5), ':', values(6)
       write (1,*) 'The seting in initial condition have errors !!'
       write (1,*) 'nuex =', nuex
       write (1,*) "nuex don't be equal to 1"
       write (1,*) 'Please check nuex !!'
       write (1,*) 'Program STOP at get_init_data !!'
       close (1)
       write (6,99) values(1), '-', values(2), '-', values(3), ' ', &
                    values(5), ':', values(6)
       write (6,*) 'The seting in initial condition have errors !!'
       write (6,*) 'nuex =', nuex
       write (6,*) "nuex don't be equal to 1"
       write (6,*) 'Please check nuex !!'
       write (6,*) 'Program STOP at get_init_data !!'
       close (3)
       Istop = 999
       return
     endif
!-------------------------------------------------------------------------------
     if (nuey == 1 .and. nuez /= 1) then
       open (1,file='error.txt',status='unknown')
       write (1,99) values(1), '-', values(2), '-', values(3), ' ', &
                    values(5), ':', values(6)
       write (1,*) 'The seting in initial condition have errors !!'
       write (1,*) 'nuey =', nuey, ' nuez =', nuez
       write (1,*) 'If nuey is equal to 1, nuez must be equal to 1'
       write (1,*) 'Please check nuey and nuez !!'
       write (1,*) 'Program STOP at get_init_data !!'
       close (1)
       write (6,99) values(1), '-', values(2), '-', values(3), ' ', &
                    values(5), ':', values(6)
       write (6,*) 'The seting in initial condition have errors !!'
       write (6,*) 'nuey =', nuey, ' nuez =', nuez
       write (6,*) 'If nuey is equal to 1, nuez must be equal to 1'
       write (6,*) 'Please check nuey and nuez !!'
       write (6,*) 'Program STOP at get_init_data !!'
       close (3)
       Istop = 999
       return
     endif
!-------------------------------------------------------------------------------
     if (nuix == 1) then
       open (1,file='error.txt',status='unknown')
       write (1,99) values(1), '-', values(2), '-', values(3), ' ', &
                    values(5), ':', values(6)
       write (1,*) 'The seting in initial condition have errors !!'
       write (1,*) 'nuix =', nuix
       write (1,*) "nuix don't be equal to 1"
       write (1,*) 'Please check nuix !!'
       write (1,*) 'Program STOP at get_init_data !!'
       close (1)
       write (6,99) values(1), '-', values(2), '-', values(3), ' ', &
                    values(5), ':', values(6)
       write (6,*) 'The seting in initial condition have errors !!'
       write (6,*) 'nuix =', nuix
       write (6,*) "nuix don't be equal to 1"
       write (6,*) 'Please check nuix !!'
       write (6,*) 'Program STOP at get_init_data !!'
       close (3)
       Istop = 999
       return
     endif
!-------------------------------------------------------------------------------
     if (nuiy == 1 .and. nuiz /= 1) then
       open (1,file='error.txt',status='unknown')
       write (1,99) values(1), '-', values(2), '-', values(3), ' ', &
                    values(5), ':', values(6)
       write (1,*) 'The seting in initial condition have errors !!'
       write (1,*) 'nuiy =', nuiy, ' nuiz =', nuiz
       write (1,*) 'If nuiy is equal to 1, nuiz must be equal to 1'
       write (1,*) 'Please check nuiy and nuiz !!'
       write (1,*) 'Program STOP at get_init_data !!'
       close (1)
       write (6,99) values(1), '-', values(2), '-', values(3), ' ', &
                    values(5), ':', values(6)
       write (6,*) 'The seting in initial condition have errors !!'
       write (6,*) 'nuiy =', nuiy, ' nuiz =', nuiz
       write (6,*) 'If nuiy is equal to 1, nuiz must be equal to 1'
       write (6,*) 'Please check nuiy and nuiz !!'
       write (6,*) 'Program STOP at get_init_data !!'
       close (3)
       Istop = 999
       return
     endif
!-------------------------------------------------------------------------------
     read (3) x_global, y_global, z_global, uex, uey, uez, uix, uiy, uiz
     close (3)
!
  99 format (1x,i4,4(a1,i2.2))
   end subroutine get_init_data

!===============================================================================
   subroutine current0
     implicit double precision (A-H,O-Z)
!-------------------------------------------------------------------------------
!  get the current0 of electron
!
     call current(ncx_mpi,ncy_mpi,ncz_mpi,nuex,nuey,nuez,f(ncfe),vex,vey,vez, &
                  Buex,Cuex,huex,Buey,Cuey,huey,Buez,Cuez,huez,wkJex,wkJey,wkJez)
!-------------------------------------------------------------------------------
!  get the current0 of ion
!
     call current(ncx_mpi,ncy_mpi,ncz_mpi,nuix,nuiy,nuiz,f(ncfi),vix,viy,viz, &
                  Buix,Cuix,huix,Buiy,Cuiy,huiy,Buiz,Cuiz,huiz,wkJix,wkJiy,wkJiz)
!-------------------------------------------------------------------------------
     wkJx0 = wkJix-wkJex
     wkJy0 = wkJiy-wkJey
     wkJz0 = wkJiz-wkJez   
   end subroutine current0
!===============================================================================
   subroutine get_fluid
     implicit double precision (A-H,O-Z)
!
!  get the plasma density
!
     call density(ncx_mpi,ncy_mpi,ncz_mpi,nuex,nuey,nuez,f(ncfe), &
                  Buex,Cuex,huex,Buey,Cuey,huey,Buez,Cuez,huez,Rhoe)
     call density(ncx_mpi,ncy_mpi,ncz_mpi,nuix,nuiy,nuiz,f(ncfi), &
                  Buix,Cuix,huix,Buiy,Cuiy,huiy,Buiz,Cuiz,huiz,Rhoi)
!-------------------------------------------------------------------------------
!  get the average velocity of electron
!   
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call current(ncx_mpi,ncy_mpi,ncz_mpi,nuex,nuey,nuez,f(ncfe),uex3D,uey3D,uez3D, &
                  Buex,Cuex,huex,Buey,Cuey,huey,Buez,Cuez,huez,avguex,avguey,avguez)
!
     avguex = avguex/Rhoe
     avguey = avguey/Rhoe
     avguez = avguez/Rhoe
!
     call current(ncx_mpi,ncy_mpi,ncz_mpi,nuex,nuey,nuez,f(ncfe),vex,vey,vez, &
                  Buex,Cuex,huex,Buey,Cuey,huey,Buez,Cuez,huez,avgvex,avgvey,avgvez)
!
     avgvex = avgvex/Rhoe
     avgvey = avgvey/Rhoe
     avgvez = avgvez/Rhoe
!-------------------------------------------------------------------------------
!  get the average velocity of ion
!   
     call current(ncx_mpi,ncy_mpi,ncz_mpi,nuix,nuiy,nuiz,f(ncfi),uix3D,uiy3D,uiz3D, &
                  Buix,Cuix,huix,Buiy,Cuiy,huiy,Buiz,Cuiz,huiz,avguix,avguiy,avguiz)
!
     avguix = avguix/Rhoi
     avguiy = avguiy/Rhoi
     avguiz = avguiz/Rhoi
!
     call current(ncx_mpi,ncy_mpi,ncz_mpi,nuix,nuiy,nuiz,f(ncfi),vix,viy,viz, &
                  Buix,Cuix,huix,Buiy,Cuiy,huiy,Buiz,Cuiz,huiz,avgvix,avgviy,avgviz)
!
     avgvix = avgvix/Rhoi
     avgviy = avgviy/Rhoi
     avgviz = avgviz/Rhoi
!-------------------------------------------------------------------------------
     call pressure(ncx_mpi,ncy_mpi,ncz_mpi,nuex,nuey,nuez,f(ncfe),1.d0, &
                   Buex,Cuex,huex,Buey,Cuey,huey,Buez,Cuez,huez, &
                   uex,uey,uez,avguex,avguey,avguez, &
                   vex,vey,vez,avgvex,avgvey,avgvez, &
                   Pex,Pey,Pez)
!          
     call pressure(ncx_mpi,ncy_mpi,ncz_mpi,nuix,nuiy,nuiz,f(ncfi),ami, &
                   Buix,Cuix,huix,Buiy,Cuiy,huiy,Buiz,Cuiz,huiz, &
                   uix,uiy,uiz,avguix,avguiy,avguiz, &
                   vix,viy,viz,avgvix,avgviy,avgviz, &
                   Pix,Piy,Piz)
   end subroutine get_fluid
!===============================================================================
   subroutine data_decomposition_1D(f_recv,f_bcst,nc_mpi,i_start)
     implicit double precision (A-H,O-Z) 
     real*8, allocatable, dimension(:) :: ftemp
     dimension :: index(2), f_recv(1), f_bcst(1)
!     
     if (MyID == 0) then
!
       do i = 1, nc_mpi
         f_recv(i) = f_bcst(i)      
       enddo
!
       do id = 1, NP-1
         itag = 10*id
!
         call MPI_RECV(index(1),2,MPI_INTEGER,id,itag,MPI_WORLD,ISTATUS,IERR)
!
         is = index(1)
         NN  = index(2)
         allocate (ftemp(NN))
!
         ii = 0
         do 20 i = is, is+NN-1
           ii = ii+1
           ftemp(ii) = f_bcst(i)
    20   continue
!
         call MPI_SEND(ftemp(1),NN,MPI_REAL8,id,itag+1,MPI_WORLD,IERR)        
!
         deallocate (ftemp)
       enddo
!-------------------------------------------------------------------------------
     else
       index(1) = i_start
       index(2) = nc_mpi
       itag = 10*MyID
       call MPI_SEND(index(1),2,MPI_INTEGER,0,itag,MPI_WORLD,IERR)
       call MPI_RECV(f_recv(1),nc_mpi,MPI_REAL8,0,itag+1, &
                     MPI_WORLD,ISTATUS,IERR)
     endif  
!     
   end subroutine data_decomposition_1D
!===============================================================================
!===============================================================================
   subroutine data_decomposition_6D(f_bcst,f_recv,nux,nuy,nuz)
     implicit double precision (A-H,O-Z)
     real*8, allocatable, dimension(:) :: ftemp
     dimension :: index(6), f_recv(1), f_bcst(1)
     nuxy  = nux*nuy
     nuxyz = nux*nuy*nuz
!     
     if (MyID == 0) then
!
       ii = 0 !istart_local-1
       do 10 k = 1, ncz_mpi
       do 10 j = 1, ncy_mpi
       do 10 i = 1, ncx_mpi
!
         ijk = i+(j-1)*ncx+(k-1)*ncxy
!
         do 11 ku = 1, nuz
         do 11 ju = 1, nuy
         do 11 iu = 1, nux           
           ijku = iu+(ju-1)*nux+(ku-1)*nuxy
           ijkg = ijku+(ijk-1)*nuxyz !+istart_global-1 
           ii = ii+1 
           f_recv(ii) = f_bcst(ijkg)
  11     continue      
  10   continue       
!
       do id = 1, NP-1
         itag = 10*id
!
         call MPI_RECV(index(1),6,MPI_INTEGER,id,itag,MPI_WORLD,ISTATUS,IERR)
!
         ixs = index(1)
         iys = index(2)
         izs = index(3)
         Nx  = index(4)
         Ny  = index(5)
         Nz  = index(6)
!         Nxy = Nx*Ny
         Nxyz = Nx*Ny*Nz
         NNxyz = Nxyz*nuxyz
         allocate (ftemp(NNxyz))
!
         ii = 0
         do 20 k = izs, izs+Nz-1
         do 20 j = iys, iys+Ny-1
         do 20 i = ixs, ixs+Nx-1
!
           ijk = i+(j-1)*ncx+(k-1)*ncxy
!
           do 21 ku = 1, nuz
           do 21 ju = 1, nuy
           do 21 iu = 1, nux
             ijku = iu+(ju-1)*nux+(ku-1)*nuxy
             ijkg = ijku+(ijk-1)*nuxyz !+istart_global-1 
             ii = ii+1
             ftemp(ii) = f_bcst(ijkg)
    21     continue
    20   continue
!
         call MPI_SEND(ftemp(1),NNxyz,MPI_REAL8,id,itag+1,MPI_WORLD,IERR)        
         deallocate (ftemp)
       enddo
!-------------------------------------------------------------------------------
     else
       index(1) = ix_start
       index(2) = iy_start
       index(3) = iz_start
       index(4) = ncx_mpi
       index(5) = ncy_mpi
       index(6) = ncz_mpi
       NNxyz = ncxyz_mpi*nuxyz
       itag = 10*MyID
       call MPI_SEND(index(1),6,MPI_INTEGER,0,itag,MPI_WORLD,IERR)
       call MPI_RECV(f_recv,NNxyz,MPI_REAL8,0,itag+1, &
                     MPI_WORLD,ISTATUS,IERR)
     endif  
!     
   end subroutine data_decomposition_6D
!===============================================================================
!===============================================================================
   subroutine data_decomposition(f_bcst,f_recv)
     implicit double precision (A-H,O-Z)
     real*8, allocatable, dimension(:) :: ftemp
     dimension :: index(6), f_recv(1), f_bcst(1)
!     
     if (MyID == 0) then
!
       ii = 0 !istart_local-1
       do 10 k = 1, ncz_mpi
       do 10 j = 1, ncy_mpi
       do 10 i = 1, ncx_mpi
         ii = ii+1
         ijk = i+(j-1)*ncx+(k-1)*ncxy !+istart_global-1 
         f_recv(ii) = f_bcst(ijk)      
  10   continue
!
       do id = 1, NP-1
         itag = 10*id
!
         call MPI_RECV(index(1),6,MPI_INTEGER,id,itag,MPI_WORLD,ISTATUS,IERR)
!
         ixs = index(1)
         iys = index(2)
         izs = index(3)
         Nx  = index(4)
         Ny  = index(5)
         Nz  = index(6)
!         Nxy = Nx*Ny
         Nxyz = Nx*Ny*Nz
         allocate (ftemp(Nxyz))
!
         ii = 0
         do 20 k = izs, izs+Nz-1
         do 20 j = iys, iys+Ny-1
         do 20 i = ixs, ixs+Nx-1
           ii = ii+1
           ijk = i+(j-1)*ncx+(k-1)*ncxy !+istart_global-1
           ftemp(ii) = f_bcst(ijk)
    20   continue
!
         call MPI_SEND(ftemp(1),Nxyz,MPI_REAL8,id,itag+1,MPI_WORLD,IERR)        
         deallocate (ftemp)
       enddo
!-------------------------------------------------------------------------------
     else
       index(1) = ix_start
       index(2) = iy_start
       index(3) = iz_start
       index(4) = ncx_mpi
       index(5) = ncy_mpi
       index(6) = ncz_mpi
       itag = 10*MyID
       call MPI_SEND(index(1),6,MPI_INTEGER,0,itag,MPI_WORLD,IERR)
       call MPI_RECV(f_recv,ncxyz_mpi,MPI_REAL8,0,itag+1, &
                     MPI_WORLD,ISTATUS,IERR)
     endif  
!     
   end subroutine data_decomposition
!===============================================================================
!===============================================================================
   subroutine combine_data_6D(f_bcst,f_recv,nux,nuy,nuz)
!-------------------------------------------------------------------------------
!
!  combine the data to the global array in MYID = 0
!
!  working arrays : index(6)
!
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     real*8, allocatable, dimension(:) :: ftemp
     dimension :: index(6), f_recv(1), f_bcst(1)
     nuxy  = nux*nuy
     nuxyz = nux*nuy*nuz
!
     if (MyID == 0) then
!
       ii = 0
       do 10 k = 1, ncz_mpi
       do 10 j = 1, ncy_mpi
       do 10 i = 1, ncx_mpi
!
         ijk = i+(j-1)*ncx+(k-1)*ncxy
!         
         do 11 ku = 1, nuz
         do 11 ju = 1, nuy
         do 11 iu = 1, nux           
           ijku = iu+(ju-1)*nux+(ku-1)*nuxy
           ijkg = ijku+(ijk-1)*nuxyz !+istart_global-1 
           ii = ii+1
           f_recv(ijkg) = f_bcst(ii)
  11     continue      
  10   continue
!
       do id = 1, NP-1
         itag = 10*id
!
         call MPI_RECV(index(1),6,MPI_INTEGER,id,itag,MPI_WORLD,ISTATUS,IERR)
!
         Nx = index(4)
         Ny = index(5)
         Nz = index(6)
!         Nxy = Nx*Ny
         Nxyz = Nx*Ny*Nz
         NNxyz = Nxyz*nuxyz
         allocate (ftemp(NNxyz))
!
         call MPI_RECV(ftemp(1),NNxyz,MPI_REAL8,id,itag+1,MPI_WORLD,ISTATUS,IERR)
!
         ixs = index(1)
         iys = index(2)
         izs = index(3)
         ii = 0
         do 20 k = izs, izs+Nz-1
         do 20 j = iys, iys+Ny-1
         do 20 i = ixs, ixs+Nx-1
!
           ijk = i+(j-1)*ncx+(k-1)*ncxy !+istart_global-1
!
           do 21 ku = 1, nuz
           do 21 ju = 1, nuy
           do 21 iu = 1, nux
             ijku = iu+(ju-1)*nux+(ku-1)*nuxy
             ijkg = ijku+(ijk-1)*nuxyz !+istart_global-1 
             ii = ii+1
             f_recv(ijkg) = ftemp(ii)
    21     continue           
    20   continue
!
         deallocate (ftemp)
       enddo
!-------------------------------------------------------------------------------
     else
       index(1) = ix_start
       index(2) = iy_start
       index(3) = iz_start
       index(4) = ncx_mpi
       index(5) = ncy_mpi
       index(6) = ncz_mpi
       NNxyz = ncxyz_mpi*nuxyz
       itag = 10*MyID
       call MPI_SEND(index(1),6,MPI_INTEGER,0,itag,MPI_WORLD,IERR)
       call MPI_SEND(f_bcst,NNxyz,MPI_REAL8,0,itag+1,MPI_WORLD,IERR)
     endif
     call MPI_Barrier(MPI_COMM_WORLD,ierr) 
   end subroutine combine_data_6D
!===============================================================================
!===============================================================================
   subroutine combine_data(f_bcst,f_recv)
!-------------------------------------------------------------------------------
!
!  combine the data to the global array in MYID = 0
!
!  working arrays : index(6)
!
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     real*8, allocatable, dimension(:) :: ftemp !, fcom
     dimension :: index(6), f_recv(1), f_bcst(1)
!  
     if (MyID == 0) then
!
       ii = 0
       do 10 k = 1, ncz_mpi
       do 10 j = 1, ncy_mpi
       do 10 i = 1, ncx_mpi
         ii = ii+1
         ijk = i+(j-1)*ncx+(k-1)*ncxy !+istart_global-1
         f_recv(ijk) = f_bcst(ii)      
  10   continue
!
       do id = 1, NP-1
         itag = 10*id
!
         call MPI_RECV(index(1),6,MPI_INTEGER,id,itag,MPI_WORLD,ISTATUS,IERR)
!
         Nx = index(4)
         Ny = index(5)
         Nz = index(6)
!         Nxy = Nx*Ny
         Nxyz = Nx*Ny*Nz
         allocate (ftemp(Nxyz))
!
         call MPI_RECV(ftemp(1),Nxyz,MPI_REAL8,id,itag+1,MPI_WORLD,ISTATUS,IERR)
!
         ixs = index(1)
         iys = index(2)
         izs = index(3)
         ii = 0
         do 20 k = izs, izs+Nz-1
         do 20 j = iys, iys+Ny-1
         do 20 i = ixs, ixs+Nx-1
           ii = ii+1
           ijk = i+(j-1)*ncx+(k-1)*ncxy !+istart_global-1
           f_recv(ijk) = ftemp(ii)
    20   continue
!
         deallocate (ftemp)
       enddo
!-------------------------------------------------------------------------------
     else
       index(1) = ix_start
       index(2) = iy_start
       index(3) = iz_start
       index(4) = ncx_mpi
       index(5) = ncy_mpi
       index(6) = ncz_mpi
       itag = 10*MyID
       call MPI_SEND(index(1),6,MPI_INTEGER,0,itag,MPI_WORLD,IERR)
       call MPI_SEND(f_bcst,ncxyz_mpi,MPI_REAL8,0,itag+1,MPI_WORLD,IERR)
     endif 
   end subroutine combine_data
!===============================================================================
end module bcstrecv