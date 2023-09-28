!===============================================================================
!  prepare the parameters for MPI
!
module MPI3D
    use coef
    use comm_global
    use omp_lib
    include 'mpif.h'
!-------------------------------------------------------------------------------
!  NP   : the total number of nodes
!  MyID : the node ID
!  ip, iq, ir : the equivalent node ID in the x, y, z direction
!  ncx_mpi : the local array dimension in x direction
!  ncy_mpi : the local array dimension in y direction
!  ncz_mpi : the local array dimension in z direction
!  ixs : the local start point of the x-axis in this node
!  iys : the local start point of the y-axis in this node
!  izs : the local start point of the z-axis in this node
!
    integer :: NP, MyID, MPI_WORLD, Istop, IERR
    integer :: ISTATUS(MPI_STATUS_SIZE)
    integer :: ip, iq, ir
    integer :: ipm1, ipp1, ipstart, ipend
    integer :: iqm1, iqp1, iqstart, iqend
    integer :: irm1, irp1, irstart, irend
    integer :: ncx_mpi, ncy_mpi, ncz_mpi
    integer :: ncxy_mpi, ncyz_mpi, ncxyz_mpi
    integer :: ix_start, iy_start, iz_start    
  contains
!===============================================================================    
  subroutine MPI3D_init
!
!  start the MPI
!
    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, NP,  IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD, MyID,IERR)
    MPI_WORLD = MPI_COMM_WORLD
!
!  read the global parameter
!
    if (MyID == 0) call get_para
    call MPI_bcast(ami,1,MPI_real8,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(C,1,MPI_real8,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(v_b,1,MPI_real8,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(ncx,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(ncy,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(ncz,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(nuex,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(nuey,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(nuez,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(nuix,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(nuiy,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(nuiz,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(ixp,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(iyp,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(izp,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(P  ,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(Q  ,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(R  ,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(it_plot_d,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(it_plot_f,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(it_plot_s,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(it_stop,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(irel,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(iter,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(ccheck,1,MPI_character,0,MPI_COMM_WORLD,IERR)
	call MPI_bcast(atype,2,MPI_character,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(mtype,1,MPI_character,0,MPI_COMM_WORLD,IERR)
	call MPI_bcast(rtype,1,MPI_character,0,MPI_COMM_WORLD,IERR)
!    call MPI_bcast(data_format,6,MPI_character,0,MPI_COMM_WORLD,IERR)
    call MPI_bcast(iOMP,1,MPI_integer,0,MPI_COMM_WORLD,IERR)
!
!  prepare the parameter for the MPI
!
    call check_MPI
    call MPI3D_decomposition
	call check_init_5th
!
  end subroutine MPI3D_init
!===============================================================================
  subroutine MPI3D_FINALIZE
    call MPI_FINALIZE(IERR)
  end subroutine MPI3D_FINALIZE
!===============================================================================
  subroutine check_MPI
    character*10 :: dates, times, zone
    integer, dimension(:) :: values(8)
!
    call Date_and_Time(dates,times,zone,values)
!-------------------------------------------------------------------------------
    if (C <= 0.d0) then
      if (MyID == 0) then
        open (1,file='error.txt',status='unknown')
        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (1,*) 'The seting in MPI have errors !!'
        write (1,*) 'C =', C, '. Light speed C must be great than 0.'
        write (1,*) 'Please check C !!'
        write (1,*) 'Program STOP at check_MPI !!'
        close (1)
        write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (6,*) 'The seting in MPI have errors !!'
        write (6,*) 'C =', C, '. Light speed C must be great than 0.'
        write (6,*) 'Please check C !!'
        write (6,*) 'Program STOP at check_MPI !!'
      endif
      call MPI_FINALIZE(IERR)
      STOP
    endif
!-------------------------------------------------------------------------------
    if (ami <= 1.d0) then
      if (MyID == 0) then
        open (1,file='error.txt',status='unknown')
        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (1,*) 'The seting in MPI have errors !!'
        write (1,*) 'ami =', ami, '. Ion mass ratio ami must be great than 0.'
        write (1,*) 'Please check C !!'
        write (1,*) 'Program STOP at check_MPI !!'
        close (1)
        write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (6,*) 'The seting in MPI have errors !!'
        write (6,*) 'ami =', ami, '. Ion mass ratio ami must be great than 0.'
        write (6,*) 'Please check C !!'
        write (6,*) 'Program STOP at check_MPI !!'
      endif
      call MPI_FINALIZE(IERR)
      STOP
    endif
!-------------------------------------------------------------------------------
    if (NP /= P*Q*R) then
      if (MyID == 0) then
        open (1,file='error.txt',status='unknown')
        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (1,*) 'The seting in MPI have errors !!'
        write (1,*) 'Please check P, Q, R !!'
        write (1,*) 'Program STOP at check_MPI !!'
        close (1)
        write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (6,*) 'The seting in MPI have errors !!'
        write (6,*) 'Please check P, Q, R !!'
        write (6,*) 'Program STOP at check_MPI !!'
      endif
      call MPI_FINALIZE(IERR)
      STOP
    endif
!-------------------------------------------------------------------------------
    if (ncx < P) then
      if (MyID == 0) then
        open (1,file='error.txt',status='unknown')
        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (1,*) 'The seting in MPI have errors !!'
        write (1,*) 'ncx =', ncx, ' <  P =', P
        write (1,*) 'Please check P and ncx !!'
        write (1,*) 'Program STOP at check_MPI !!'
        close (1)
        write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (6,*) 'The seting in MPI have errors !!'
        write (6,*) 'ncx =', ncx, ' <  P =', P
        write (6,*) 'Please check P and ncx !!'
        write (6,*) 'Program STOP at check_MPI !!'
      endif
      call MPI_FINALIZE(IERR)
      STOP
    endif
!-------------------------------------------------------------------------------
    if (ncy < Q) then
      if (MyID == 0) then
        open (1,file='error.txt',status='unknown')
        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (1,*) 'The seting in MPI have errors !!'
        write (1,*) 'ncy =', ncy, ' <  Q =', Q
        write (1,*) 'Please check Q and ncy !!'
        write (1,*) 'Program STOP at check_MPI !!'
        close (1)
        write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (6,*) 'The seting in MPI have errors !!'
        write (6,*) 'ncy =', ncy, ' <  Q =', Q
        write (6,*) 'Please check Q and ncy !!'
        write (6,*) 'Program STOP at check_MPI !!'
      endif
      call MPI_FINALIZE(IERR)
      STOP
    endif
!-------------------------------------------------------------------------------
    if (ncy == 1 .and. ncz /= 1) then
      if (MyID == 0) then
        open (1,file='error.txt',status='unknown')
        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (1,*) 'The seting in MPI have errors !!'
        write (1,*) 'ncy =', ncy, ' ncz =', ncz
        write (1,*) 'If ncy is equal to 1, ncz must be equal to 1'
        write (1,*) 'Please check ncy and ncz !!'
        write (1,*) 'Program STOP at check_MPI !!'
        close (1)
        write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (6,*) 'The seting in MPI have errors !!'
        write (6,*) 'ncy =', ncy, ' ncz =', ncz
        write (6,*) 'If ncy is equal to 1, ncz must be equal to 1'
        write (6,*) 'Please check ncy and ncz !!'
        write (6,*) 'Program STOP at check_MPI !!'
      endif
      call MPI_FINALIZE(IERR)
      STOP
    endif
!-------------------------------------------------------------------------------
    if (ncz < R) then
      if (MyID == 0) then
        open (1,file='error.txt',status='unknown')
        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (1,*) 'The seting in MPI have errors !!'
        write (1,*) 'ncz =', ncz, ' <  R =', R
        write (1,*) 'Please check R and ncz !!'
        write (1,*) 'Program STOP at check_MPI !!'
        close (1)
        write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (6,*) 'The seting in MPI have errors !!'
        write (6,*) 'ncz =', ncz, ' <  R =', R
        write (6,*) 'Please check R and ncz !!'
        write (6,*) 'Program STOP at check_MPI !!'
      endif
      call MPI_FINALIZE(IERR)
      STOP
    endif
!-------------------------------------------------------------------------------
    if (it_plot_f < 2 .or. it_plot_d < 2) then
      if (MyID == 0) then
        open (1,file='error.txt',status='unknown')
        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (1,*) 'The seting in MPI have errors !!'
        write (1,*) 'it_plot_f =', it_plot_f, ' < 2'
        write (1,*) 'it_plot_d =', it_plot_d, ' < 2'
        write (1,*) 'Please check the setting of it_plot !!'
        write (1,*) 'Program STOP at check_MPI !!'
        close (1)
        write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (6,*) 'The seting in MPI have errors !!'
        write (6,*) 'it_plot_f =', it_plot_f, ' < 2'
        write (6,*) 'it_plot_d =', it_plot_d, ' < 2'
        write (6,*) 'Please check the setting of it_plot !!'
        write (6,*) 'Program STOP at check_MPI !!'
      endif
      call MPI_FINALIZE(IERR)
      STOP
    endif
!-------------------------------------------------------------------------------
    if (it_stop < 5) then
      if (MyID == 0) then
        open (1,file='error.txt',status='unknown')
        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (1,*) 'The seting in MPI have errors !!'
        write (1,*) 'it_stop =', it_stop, ' < 5'
        write (1,*) 'Please check the setting of it_stop !!'
        write (1,*) 'Program STOP at check_MPI !!'
        close (1)
        write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (6,*) 'The seting in MPI have errors !!'
        write (6,*) 'it_stop =', it_stop, ' < 5'
        write (6,*) 'Please check the setting of it_stop !!'
        write (6,*) 'Program STOP at check_MPI !!'
      endif
      call MPI_FINALIZE(IERR)
      STOP
    endif
!-------------------------------------------------------------------------------
    if (irel < 20) then
      if (MyID == 0) then
        open (1,file='error.txt',status='unknown')
        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (1,*) 'The seting in MPI have errors !!'
        write (1,*) 'convergence condition =', irel ,' < 20'
        write (1,*) 'Please check the setting of convergence condition !!'
        write (1,*) 'Program STOP at check_MPI !!'
        close (1)
        write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (6,*) 'The seting in MPI have errors !!'
        write (6,*) 'convergence condition =', irel ,' < 20'
        write (6,*) 'Please check the setting of convergence condition !!'
        write (6,*) 'Program STOP at check_MPI !!'
      endif
      call MPI_FINALIZE(IERR)
      STOP
    endif
!-------------------------------------------------------------------------------
    if (iter < 0) then
      if (MyID == 0) then
        open (1,file='error.txt',status='unknown')
        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (1,*) 'The seting in MPI have errors !!'
        write (1,*) 'iteration times =', iter ,' < 0'
        write (1,*) 'Please check the setting of iteration times !!'
        write (1,*) 'Program STOP at check_MPI !!'
        close (1)
        write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (6,*) 'The seting in MPI have errors !!'
        write (6,*) 'iteration times =', iter ,' < 0'
        write (6,*) 'Please check the setting of iteration times !!'
        write (6,*) 'Program STOP at check_MPI !!'
      endif
      call MPI_FINALIZE(IERR)
      STOP
    endif
!-------------------------------------------------------------------------------
    if (ccheck == 'Y' .or. ccheck == 'N') then
      goto 111
    else
      if (MyID == 0) then
        open (1,file='error.txt',status='unknown')
        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (1,*) 'The seting in MPI have errors !!'
        write (1,*) 'check convergence or not =', ccheck
        write (1,*) 'Please check the setting of it_stop !!'
        write (1,*) 'Program STOP at check_MPI !!'
        close (1)
        write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (6,*) 'The seting in MPI have errors !!'
        write (6,*) 'check convergence or not =', ccheck
        write (6,*) 'Please check the setting of it_stop !!'
        write (6,*) 'Program STOP at check_MPI !!'
      endif
      call MPI_FINALIZE(IERR)
      STOP
    endif
!
111 continue
!-------------------------------------------------------------------------------
    if (mtype == 'Y' .or. mtype == 'N') then
      goto 112
    else
      if (MyID == 0) then
        open (1,file='error.txt',status='unknown')
        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (1,*) 'The seting in initial condition have error !!'
        write (1,*) 'ion movable =', mtype
        write (1,*) 'Please check the setting of ion movable !!'
        write (1,*) 'Program STOP at check_MPI !!'
        close (1)
        write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (6,*) 'The seting in initial condition have error !!'
        write (6,*) 'ion movable =', mtype
        write (6,*) 'Please check the setting of ion movable !!'
        write (6,*) 'Program STOP at check_MPI !!'
      endif
      call MPI_FINALIZE(IERR)
      STOP
    endif
!
112 continue
!-------------------------------------------------------------------------------
!    if (data_format == 'single' .or. data_format == 'double') then
!    	goto 113
!    else
!      if (MyID == 0) then
!        open (1,file='error.txt',status='unknown')
!        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
!                     values(5), ':', values(6)
!        write (1,*) 'The seting in MPI have errors !!'
!        write (1,*) 'data_format =', data_format
!        write (1,*) 'Please check the setting of data_format !!'
!        write (1,*) 'Program STOP at check_MPI !!'
!        close (1)
!        write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
!                     values(5), ':', values(6)
!        write (6,*) 'The seting in MPI have errors !!'
!        write (6,*) 'data_format =', data_format
!        write (6,*) 'Please check the setting of data_format !!'
!        write (6,*) 'Program STOP at check_MPI !!'
!      endif
!      call MPI_FINALIZE(IERR)
!      STOP
!    endif
!
!113 continue
!-------------------------------------------------------------------------------
    CALL OMP_SET_NUM_THREADS(iOMP)
!-------------------------------------------------------------------------------
 11 format (1x,i4,4(a1,i2.2))
  end subroutine check_MPI
!===============================================================================
  subroutine check_init_5th
    character*10 :: dates, times, zone
    integer, dimension(:) :: values(8)
!
    call Date_and_Time(dates,times,zone,values)
!-------------------------------------------------------------------------------
    if (ncx_mpi < 7) then
      if (MyID == 0) then
        open (1,file='error.txt',status='unknown')
        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (1,*) 'The seting in MPI have errors !!'
        write (1,*) 'ncx_mpi =', ncx_mpi, ' <  7'
        write (1,*) 'Please check P and ncx !!'
        write (1,*) 'Program STOP at check_init_5th !!'
        close (1)
        write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (6,*) 'The seting in MPI have errors !!'
        write (6,*) 'ncx_mpi =', ncx_mpi, ' <  7'
        write (6,*) 'Please check P and ncx !!'
        write (6,*) 'Program STOP at check_init_5th !!'
      endif
      call MPI_FINALIZE(IERR)
      STOP
    endif
!
    if (ncy_mpi == 1) goto 777
    if (ncy_mpi < 7) then
      if (MyID == 0) then
        open (1,file='error.txt',status='unknown')
        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (1,*) 'The seting in MPI have errors !!'
        write (1,*) 'ncy_mpi =', ncy_mpi, ' <  7'
        write (1,*) 'Please check Q and ncy !!'
        write (1,*) 'Program STOP at check_init_5th !!'
        close (1)
        write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (6,*) 'The seting in MPI have errors !!'
        write (6,*) 'ncy_mpi =', ncy_mpi, ' <  7'
        write (6,*) 'Please check Q and ncy !!'
        write (6,*) 'Program STOP at check_init_5th !!'
      endif
      call MPI_FINALIZE(IERR)
      STOP
    endif
777 continue
!
    if (ncz_mpi == 1) goto 888
    if (ncz_mpi < 7) then
      if (MyID == 0) then
        open (1,file='error.txt',status='unknown')
        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (1,*) 'The seting in MPI have errors !!'
        write (1,*) 'ncz_mpi =', ncz_mpi, ' <  7'
        write (1,*) 'Please check R and ncz !!'
        write (1,*) 'Program STOP at check_init_5th !!'
        close (1)
        write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (6,*) 'The seting in MPI have errors !!'
        write (6,*) 'ncz_mpi =', ncz_mpi, ' <  7'
        write (6,*) 'Please check R and ncz !!'
        write (6,*) 'Program STOP at check_init_5th !!'
      endif
      call MPI_FINALIZE(IERR)
      STOP
    endif
888 continue
!-------------------------------------------------------------------------------
 11 format (1x,i4,4(a1,i2.2))
  end subroutine check_init_5th
!===============================================================================
!===============================================================================
  subroutine MPI3D_decomposition
!-------------------------------------------------------------------------------
!  calculate the node ID : ip, iq, ir
!
    ir = (MyID+1)/(P*Q)
    ipq = (MyID+1)-ir*P*Q
    if (ipq /= 0) then
      ir = ir+1
    else
      ipq = P*Q
    endif
    iq = ipq/P
    ip = ipq-iq*P
    if (ip == 0) then
      ip = P
    else
      iq = iq+1 
    endif
!-------------------------------------------------------------------------------
!  calculate the dimension of x, y, z and 
!  get the global id of start point along x, y, z-direction in this node
!
    call getnc_mpi(ncx,ncx_mpi,ip,P,ix_start)
    call getnc_mpi(ncy,ncy_mpi,iq,Q,iy_start)
    call getnc_mpi(ncz,ncz_mpi,ir,R,iz_start)
    ncxy_mpi  = ncx_mpi*ncy_mpi
    ncyz_mpi  = ncy_mpi*ncz_mpi
    ncxyz_mpi = ncx_mpi*ncy_mpi*ncz_mpi
!-------------------------------------------------------------------------------
!  get the boundary ID
!
    call transferID(ip-1,iq,ir,ipm1)
    call transferID(ip+1,iq,ir,ipp1)
    call transferID(ip,iq-1,ir,iqm1)
    call transferID(ip,iq+1,ir,iqp1)
    call transferID(ip,iq,ir-1,irm1)
    call transferID(ip,iq,ir+1,irp1)
!-------------------------------------------------------------------------------
    if (ip == 1) ipm1 = MPI_PROC_NULL
    if (ip == P) ipp1 = MPI_PROC_NULL
    if (iq == 1) iqm1 = MPI_PROC_NULL
    if (iq == Q) iqp1 = MPI_PROC_NULL
    if (ir == 1) irm1 = MPI_PROC_NULL
    if (ir == R) irp1 = MPI_PROC_NULL
!-------------------------------------------------------------------------------
    call transferID(1,iq,ir,ipstart)
    call transferID(P,iq,ir,ipend)
    call transferID(ip,1,ir,iqstart)
    call transferID(ip,Q,ir,iqend)
    call transferID(ip,iq,1,irstart)
    call transferID(ip,iq,R,irend)
!-------------------------------------------------------------------------------  
  end subroutine MPI3D_decomposition
!===============================================================================
!  get the arrays dimension in this node
  subroutine getnc_mpi(nc,nc_mpi,iD,D,istart)
    integer :: D
    nc_mpi = nc/D
    nct = nc_mpi*D
    inc = nc-nct
    if (iD <= inc) then 
      nc_mpi = nc_mpi+1
      istart = 1+(iD-1)*nc_mpi
    else
      istart = 1+inc*(nc_mpi+1)+(iD-inc-1)*nc_mpi
    endif
  end subroutine getnc_mpi
!===============================================================================
!  input IDx, IDy, IDz and output node ID
  subroutine transferID(IDx,IDy,IDz,ID)
    ID = IDx+(IDy-1)*P+(IDz-1)*P*Q-1
  end subroutine transferID
!===============================================================================
end module MPI3D
!===============================================================================
