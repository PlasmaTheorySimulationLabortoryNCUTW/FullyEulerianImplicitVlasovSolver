!===============================================================================
!===============================================================================
   
   subroutine write_data(filename,wdata,dset_name,arraysize,judge)
    use wtime
    use MPI3D
    use comm_global
    use Vlasov
    use bcstrecv
    use hdf5
    character(len=*):: dset_name,filename
    dimension :: wdata(1)
    integer,intent(in) :: arraysize,judge
    integer :: error,comm,info,ipp,irr
    integer(HSIZE_T) :: data_dim(1),counter(1)
    integer(HSSIZE_T) :: offset(1)
    integer(HID_T) :: file_id,filespace,dset_id,memspace,plist_id
    data_dim(1) = arraysize
    comm = MPI_COMM_WORLD
    info = MPI_INFO_NULL
    mpi_size = NP
    mpi_rank = MyID
    call h5open_f(error)
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)
    !call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,error, access_prp = plist_id)
    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error, access_prp = plist_id)
    CALL h5pclose_f(plist_id, error)
    call h5screate_simple_f(1,data_dim,filespace,error)
    call h5dcreate_f(file_id,dset_name,H5T_NATIVE_DOUBLE,filespace,dset_id,error)
    call h5sclose_f(filespace,error)
    ipp = ncxyz/mpi_size
    irr = ncxyz - ipp*mpi_size
    if(mpi_rank<irr) then
        counter(1) = (ipp+1)*(arraysize/ncxyz)
        offset(1) = (mpi_rank)*(ipp+1)*(arraysize/ncxyz)
    else
        counter(1) = (ipp)*(arraysize/ncxyz)
        offset(1) = (irr*(ipp+1)+(mpi_rank-irr)*ipp)*(arraysize/ncxyz)
    endif
    CALL h5screate_simple_f(1, counter, memspace, error)
    CALL h5dget_space_f(dset_id, filespace, error)
    CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, counter, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,wdata,data_dim,error, &
                    file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)
    call h5dclose_f(dset_id,error)
    CALL h5pclose_f(plist_id, error)
    call h5fclose_f(file_id,error)
    call h5close_f(error)
   end subroutine write_data

   subroutine write_axis_data(filename,wdata,dset_name,arraysize)
    use wtime
    use MPI3D
    use comm_global
    use Vlasov
    use bcstrecv
    use hdf5
    character(len=*):: dset_name,filename
    dimension :: wdata(1)
    integer,intent(in):: arraysize
    integer :: error
    integer(HSIZE_T) :: data_dim(1), counter(1)
    integer(HSSIZE_T) :: offset(1)
    integer(HID_T) :: file_id,filespace,dset_id
    data_dim(1)=arraysize
    call h5open_f(error)
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,error)
    call h5screate_simple_f(1,data_dim,filespace,error)
    call h5dcreate_f(file_id,dset_name,H5T_NATIVE_DOUBLE,filespace,dset_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,wdata,data_dim,error)
    CALL h5sclose_f(filespace, error)
    call h5dclose_f(dset_id,error)
    call h5fclose_f(file_id,error)
    call h5close_f(error)
   end subroutine write_axis_data








   subroutine combine_Vlasov
     use wtime
     use MPI3D
     use comm_global
     use Vlasov
     use bcstrecv
     use hdf5
     implicit double precision (A-H,O-Z)
!     character*15 :: filenamed 
     character*21 :: filenamed
     character*10 :: At
     real*8, allocatable, dimension(:) :: work
     filenamed = './data/d          .h5'
     write (At,'(i10.10)') it
     filenamed(9:18) = At
!     filenamed(2:11) = At
!-------------------------------------------------------------------------------
!  output the file list in the file 'file_list.txt'
!
     if (MyID == 0) then
       if (it .eq. 0) then
         open (13,file='file_list_d.txt',status='unknown')
         write (13,120) filenamed(8:21)
         close (13)
       else
         open (13,file='file_list_d.txt',status='unknown',position='append')
         write (13,120) filenamed(8:21)
         close (13)
       endif
     endif
 120 format (a15)
!-------------------------------------------------------------------------------
!  brocast the data to the node 0, and output data.
!-------------------------------------------------------------------------------
     if (MyID == 0) then
        allocate(work(14))
        work(1) = ncx
        work(2) = ncy
        work(3) = ncz
        work(4) = nuex
        work(5) = nuey
        work(6) = nuez
        work(7) = nuix
        work(8) = nuiy
        work(9) = nuiz
        work(10)= ami
        work(11)=C
        work(12)=t
        work(13)=it
        work(14)=dt 
        call write_axis_data(filenamed,work,"axis",14)
        deallocate(work)
     endif
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!-------------------------------------------------------------------------------
!  brocast fe to the node 0, and output data.
     call write_data(filenamed,f(ncfe:ncfe_end),"fe",ncxyz*nuexyz,2)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!-------------------------------------------------------------------------------
!  brocast fi to the node 0, and output data.   
     call write_data(filenamed,f(ncfi:ncfi_end),"fi",ncxyz*nuixyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!-------------------------------------------------------------------------------
!  brocast Bx, By, Bz, Ex, Ey, Ez to the node 0

     call write_data(filenamed,f(ncBx:ncBx_end),"Bx",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamed,f(ncBy:ncBy_end),"By",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamed,f(ncBz:ncBz_end),"Bz",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamed,f(ncEx:ncEx_end),"Ex",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamed,f(ncEy:ncEy_end),"Ey",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamed,f(ncEz:ncEz_end),"Ez",ncxyz,1)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!-------------------------------------------------------------------------------
   end subroutine combine_Vlasov

   subroutine combine_fluid
     use wtime
     use MPI3D
     use comm_global
     use Vlasov
     use bcstrecv
     use hdf5
     implicit double precision (A-H,O-Z)
!     character*15 :: filenamef 
     character*21 :: filenamef
     character*10 :: At
     real*8, allocatable, dimension(:) :: work
     filenamef = './data/f          .h5'
     write (At,'(i10.10)') it
!     filenamef(2:11) = At
	 filenamef(9:18) = At
!-------------------------------------------------------------------------------
!  output the file list in the file 'file_list_f.txt'
!
     if (MyID == 0) then
       if (it .eq. 0) then
         open (12,file='file_list_f.txt',status='unknown')
         write (12,120) filenamef(8:21)
         close (12)
       else
         open (12,file='file_list_f.txt',status='unknown',position='append')
         write (12,120) filenamef(8:21)
         close (12)
       endif
     endif
 120 format (a15)
!-------------------------------------------------------------------------------
     if (MyID == 0) then
        allocate(work(14))
        work(1) = ncx
        work(2) = ncy
        work(3) = ncz
        work(4) = nuex
        work(5) = nuey
        work(6) = nuez
        work(7) = nuix
        work(8) = nuiy
        work(9) = nuiz
        work(10)= ami
        work(11)=C
        work(12)=t
        work(13)=it
        work(14)=dt 
        call write_axis_data(filenamef,work,"axis",14)
        deallocate(work)
     endif
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!-------------------------------------------------------------------------------
     call get_distribution_1D(ncx_mpi,ncy_mpi,ncz_mpi,nuex,nuey,nuez,f(ncfe), &
                              Buex,Cuex,huex,Buey,Cuey,huey,Buez,Cuez,huez,fuex,fuey,fuez) 
!
!  output fuex(nuex,ncx,ncy,ncz) data
!
     call write_data(filenamef,fuex,"fuex",ncxyz*nuex,2)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
!  output fuey(nuey,ncx,ncy,ncz) data
!
     call write_data(filenamef,fuey,"fuey",ncxyz*nuey,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
!  output fuez(nuez,ncx,ncy,ncz) data
!
     call write_data(filenamef,fuez,"fuez",ncxyz*nuez,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!-------------------------------------------------------------------------------
     call get_distribution_1D(ncx_mpi,ncy_mpi,ncz_mpi,nuix,nuiy,nuiz,f(ncfi), &
                              Buix,Cuix,huix,Buiy,Cuiy,huiy,Buiz,Cuiz,huiz,fuix,fuiy,fuiz) 
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!  output fuix(nuix,ncx,ncy,ncz) data
!
     call write_data(filenamef,fuix,"fuix",ncxyz*nuix,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
!  output fuiy(nuiy,ncx,ncy,ncz) data
!
     call write_data(filenamef,fuiy,"fuiy",ncxyz*nuiy,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
!  output fuiz(nuiz,ncx,ncy,ncz) data
!
     call write_data(filenamef,fuiz,"fuiz",ncxyz*nuiz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!-------------------------------------------------------------------------------
!  brocast Bx, By, Bz, Ex, Ey, Ez to the node 0
!

!
     call write_data(filenamef,f(ncBx:ncBx_end),"Bx",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,f(ncBy:ncBy_end),"By",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,f(ncBz:ncBz_end),"Bz",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,f(ncEx:ncEx_end),"Ex",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,f(ncEy:ncEy_end),"Ey",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,f(ncEz:ncEz_end),"Ez",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!-------------------------------------------------------------------------------
!  brocast Rho, Vx, Vy, Vz and P to the node 0
!
     call get_fluid
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,Rhoe,"Rhoe",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,Rhoi,"Rhoi",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,avgvex,"Vex",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,avgvey,"Vey",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,avgvez,"Vez",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,avgvix,"Vix",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,avgviy,"Viy",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,avgviz,"Viz",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,Pex,"Pex",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,Pey,"Pey",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,Pez,"Pez",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,Pix,"Pix",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,Piy,"Piy",ncxyz,0)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call write_data(filenamef,Piz,"Piz",ncxyz,1)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
   end subroutine combine_fluid
!===============================================================================
