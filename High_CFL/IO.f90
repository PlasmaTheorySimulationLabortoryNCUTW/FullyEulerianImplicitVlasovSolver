!===============================================================================
!===============================================================================
   subroutine combine_Vlasov
     use wtime
     use MPI3D
     use comm_global
     use Vlasov
     use bcstrecv
     implicit double precision (A-H,O-Z)
!     character*15 :: filenamed 
     character*22 :: filenamed
     character*10 :: At
     filenamed = './data/d          .bin'
     write (At,'(i10.10)') it
     filenamed(9:18) = At
!     filenamed(2:11) = At
!-------------------------------------------------------------------------------
!  output the file list in the file 'file_list.txt'
!
     if (MyID == 0) then
       if (it .eq. 0) then
         open (13,file='file_list_d.txt',status='unknown')
         write (13,120) filenamed(8:22)
         close (13)
       else
         open (13,file='file_list_d.txt',status='unknown',position='append')
         write (13,120) filenamed(8:22)
         close (13)
       endif
     endif
 120 format (a15)
!-------------------------------------------------------------------------------
!  brocast the data to the node 0, and output data.
!-------------------------------------------------------------------------------
     if (MyID == 0) then
       open (1,file=filenamed,status='unknown',form='unformatted')
       write (1) ncx, ncy, ncz, nuex, nuey, nuez, nuix, nuiy, nuiz, ami, C, t, it, dt
     endif
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!-------------------------------------------------------------------------------
!  brocast fe to the node 0, and output data.
! 
     if (MyID == 0) then
       allocate (fe_global(ncxyz*nuexyz))
     else
       allocate (fe_global(1))
     endif
!     
     call combine_data_6D(f(ncfe),fe_global,nuex,nuey,nuez)
     if (MyID == 0) write (1) fe_global
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     deallocate (fe_global)
!-------------------------------------------------------------------------------
!  brocast fi to the node 0, and output data.
! 
     if (MyID == 0) then
       allocate (fi_global(ncxyz*nuixyz))
     else
       allocate (fi_global(1))
     endif
!     
     call combine_data_6D(f(ncfi),fi_global,nuix,nuiy,nuiz)
     if (MyID == 0) write (1) fi_global
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     deallocate (fi_global)
!-------------------------------------------------------------------------------
!  brocast Bx, By, Bz, Ex, Ey, Ez to the node 0
!
     if (MyID == 0) then
       allocate (A_global(ncxyz))
     else
       allocate (A_global(1))
     endif
!
     call combine_data(f(ncBx),A_global)
     if (MyID == 0) write (1) A_global    ! write Bx
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(f(ncBy),A_global)
     if (MyID == 0) write (1) A_global    ! write By
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(f(ncBz),A_global)
     if (MyID == 0) write (1) A_global    ! write Bz
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(f(ncEx),A_global)
     if (MyID == 0) write (1) A_global    ! write Ex
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(f(ncEy),A_global)
     if (MyID == 0) write (1) A_global    ! write Ey
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(f(ncEz),A_global)
     if (MyID == 0) write (1) A_global    ! write Ez
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     deallocate (A_global)
     close (1)
!-------------------------------------------------------------------------------
   end subroutine combine_Vlasov
!===============================================================================
   subroutine combine_fluid
     use wtime
     use MPI3D
     use comm_global
     use Vlasov
     use bcstrecv
     implicit double precision (A-H,O-Z)
!     character*15 :: filenamef 
     character*22 :: filenamef
     character*10 :: At
     filenamef = './data/f          .bin'
     write (At,'(i10.10)') it
!     filenamef(2:11) = At
	 filenamef(9:18) = At
!-------------------------------------------------------------------------------
!  output the file list in the file 'file_list_f.txt'
!
     if (MyID == 0) then
       if (it .eq. 0) then
         open (12,file='file_list_f.txt',status='unknown')
         write (12,120) filenamef(8:22)
         close (12)
       else
         open (12,file='file_list_f.txt',status='unknown',position='append')
         write (12,120) filenamef(8:22)
         close (12)
       endif
     endif
 120 format (a15)
!-------------------------------------------------------------------------------
     if (MyID == 0) then
       open (2,file=filenamef,status='unknown',form='unformatted')
       write (2) ncx, ncy, ncz, ami, C, t, it, dt
     endif
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!-------------------------------------------------------------------------------
     call get_distribution_1D(ncx_mpi,ncy_mpi,ncz_mpi,nuex,nuey,nuez,f(ncfe), &
                              Buex,Cuex,huex,Buey,Cuey,huey,Buez,Cuez,huez,fuex,fuey,fuez) 
!
!  output fuex(nuex,ncx,ncy,ncz) data
!
     if (MyID == 0) then
       allocate (fuex_global(ncxyz*nuex))
     else
       allocate (fuex_global(1))
     endif
!     
     call combine_data_6D(fuex,fuex_global,nuex,1,1)
     if (MyID == 0) write (2) fuex_global
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     deallocate (fuex_global)
!
!  output fuey(nuey,ncx,ncy,ncz) data
!
     if (MyID == 0) then
       allocate (fuey_global(ncxyz*nuey))
     else
       allocate (fuey_global(1))
     endif
!     
     call combine_data_6D(fuey,fuey_global,1,nuey,1)
     if (MyID == 0) write (2) fuey_global
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     deallocate (fuey_global)
!
!  output fuez(nuez,ncx,ncy,ncz) data
!
     if (MyID == 0) then
       allocate (fuez_global(ncxyz*nuez))
     else
       allocate (fuez_global(1))
     endif
!     
     call combine_data_6D(fuez,fuez_global,1,1,nuez)
     if (MyID == 0) write (2) fuez_global
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     deallocate (fuez_global)
!-------------------------------------------------------------------------------
     call get_distribution_1D(ncx_mpi,ncy_mpi,ncz_mpi,nuix,nuiy,nuiz,f(ncfi), &
                              Buix,Cuix,huix,Buiy,Cuiy,huiy,Buiz,Cuiz,huiz,fuix,fuiy,fuiz) 
!
!  output fuix(nuix,ncx,ncy,ncz) data
!
     if (MyID == 0) then
       allocate (fuix_global(ncxyz*nuix))
     else
       allocate (fuix_global(1))
     endif
!     
     call combine_data_6D(fuix,fuix_global,nuix,1,1)
     if (MyID == 0) write (2) fuix_global
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     deallocate (fuix_global)
!
!  output fuiy(nuiy,ncx,ncy,ncz) data
!
     if (MyID == 0) then
       allocate (fuiy_global(ncxyz*nuiy))
     else
       allocate (fuiy_global(1))
     endif
!     
     call combine_data_6D(fuiy,fuiy_global,1,nuiy,1)
     if (MyID == 0) write (2) fuiy_global
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     deallocate (fuiy_global)
!
!  output fuiz(nuiz,ncx,ncy,ncz) data
!
     if (MyID == 0) then
       allocate (fuiz_global(ncxyz*nuiz))
     else
       allocate (fuiz_global(1))
     endif
!     
     call combine_data_6D(fuiz,fuiz_global,1,1,nuiz)
     if (MyID == 0) write (2) fuiz_global
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     deallocate (fuiz_global)
!-------------------------------------------------------------------------------
!  brocast Bx, By, Bz, Ex, Ey, Ez to the node 0
!
     if (MyID == 0) then
       allocate (A_global(ncxyz))
     else
       allocate (A_global(1))
     endif
!
     call combine_data(f(ncBx),A_global)
     if (MyID == 0) write (2) A_global    ! write Bx
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(f(ncBy),A_global)
     if (MyID == 0) write (2) A_global    ! write By
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(f(ncBz),A_global)
     if (MyID == 0) write (2) A_global    ! write Bz
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(f(ncEx),A_global)
     if (MyID == 0) write (2) A_global    ! write Ex
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(f(ncEy),A_global)
     if (MyID == 0) write (2) A_global    ! write Ey
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(f(ncEz),A_global)
     if (MyID == 0) write (2) A_global    ! write Ez
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!-------------------------------------------------------------------------------
!  brocast Rho, Vx, Vy, Vz and P to the node 0
!
     call get_fluid
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(Rhoe,A_global)
     if (MyID == 0) write (2) A_global  ! write Rhoe
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(Rhoi,A_global)
     if (MyID == 0) write (2) A_global  ! write Rhoi
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(avgvex,A_global)
     if (MyID == 0) write (2) A_global  ! write avgvex
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(avgvey,A_global)
     if (MyID == 0) write (2) A_global  ! write avgvey
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(avgvez,A_global)
     if (MyID == 0) write (2) A_global  ! write avgvez
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(avgvix,A_global)
     if (MyID == 0) write (2) A_global  ! write avgvix
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(avgviy,A_global)
     if (MyID == 0) write (2) A_global  ! write avgviy
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(avgviz,A_global)
     if (MyID == 0) write (2) A_global  ! write avgviz
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(Pex,A_global)
     if (MyID == 0) write (2) A_global  ! write Pex
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(Pey,A_global)
     if (MyID == 0) write (2) A_global  ! write Pey
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(Pez,A_global)
     if (MyID == 0) write (2) A_global  ! write Pez
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(Pix,A_global)
     if (MyID == 0) write (2) A_global  ! write Pix
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(Piy,A_global)
     if (MyID == 0) write (2) A_global  ! write Piy
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     call combine_data(Piz,A_global)
     if (MyID == 0) write (2) A_global  ! write Piz
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
     deallocate (A_global)
     close (2)
   end subroutine combine_fluid
!===============================================================================
