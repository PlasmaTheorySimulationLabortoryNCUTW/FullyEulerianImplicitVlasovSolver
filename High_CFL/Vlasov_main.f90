include 'comm.f90'
include 'WENO3D_MPI.f90'
include 'vector3D.f90'
include 'GPPDE.f90'
include 'Vlasov2fluid.f90'
include 'bcstrecv.f90'
include 'IO.f90'
include 'Vlasov_func.f90'
!===============================================================================
!===============================================================================
module main_loop
     use wtime
     use Vlasov
     use GPPDE
     use bcstrecv
   contains
!-------------------------------------------------------------------------------
   subroutine PCDE4(funcfp)
     implicit double precision (A-H,O-Z)
     external :: funcfp
!
     itout_f = it+it_plot_f
     itout_d = it+it_plot_d
     itout_s = it+it_plot_s
!
     !call PCDE4th0(funcfp,f,func,func2,func1,func3,wkf1,wkf2,t,dt,ntot)
     call PCDE4th05(funcfp,f,func,func1,func2,func3,func11,func12,wkf1,wkf2,t,dt,ntot)
     func3 = func3*5
     func2 = func2*5
     func1 = func1*5
!
     it = it+15
     t = dfloat(it)*dt
     call check_Vlasov
!
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     if (MyID == 0) ttt1 = MPI_Wtime()
!
     do while (it < it_stop)
!
       !if(it == itout_s) then
          !call smooth_operator(f)
          !call PCDE4th0(funcfp,f,func,func1,func2,func3,wkf1,wkf2,t,dt,ntot)
          !it = it + 3
          !itout_s = itout_s+it_plot_s
       !else
          call PCDE4th1(funcfp,f,func,func1,func2,func3,wkf1,wkf2,t,5*dt,ntot)
          it = it+5
       !endif
       call check_Vlasov
!
       t = dfloat(it)*dt
       if (MyID == 0) print *, 'time =', t      
!
!  output fluid data
!
       if (it == itout_f) then
         itout_f = itout_f+it_plot_f
         call combine_fluid
       endif
!
!  output distribution data
!
       if (it == itout_d) then
         itout_d = itout_d+it_plot_d
         call combine_Vlasov
       endif
     enddo
!
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     if (MyID == 0) then
       ttt2 = MPI_Wtime()
       open (115,file='time.txt',status='unknown')
       write (115,*) ttt1
       write (115,*) ttt2
       write (115,*) ttt2-ttt1
       close (115)
     endif         
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!   
   end subroutine PCDE4
!-------------------------------------------------------------------------------
end module main_loop
!===============================================================================
!===============================================================================
!  the main program of the 6D Vlasov simulation code
!
program Vlasov_main
!
   use wtime
   use comm_global
   use MPI3D
   use Vlasov
   use GRID6D
   use vector3D
   use GPPDE
   use bcstrecv
   use main_loop
!-------------------------------------------------------------------------------
   implicit double precision (A-H,O-Z)
   external :: Vlasov_ESim, Vlasov_ES, Vlasov_EM
   character*10 :: dates, times, zone
   integer, dimension(:) :: values(8)
   logical :: dir_e
   inquire(file='./data/.', exist=dir_e)
!===============================================================================
!  initial process of the MPI
!
   call MPI3D_init   ! in MPI3D.f90
   if (MyID == 0 .and. (dir_e .eqv. .false.)) call system('mkdir data')
   call initial      ! in bcstrecv.f90
   call check_Vlasov ! in Vlasov_main.f90
   call current0     ! in bcstrecv.f90
!
!  main loop of 3D3V Vlasov code using 5th-order iPCCFD scheme with cubic spline
!
   select case (atype) 
   case ('ES')
     if (mtype == 'N') then ! for ion non-movable              
       call Vlasov_ESim(t,f,func,ntot)  ! in Vlasov_func.f90
       call combine_fluid  ! in IO.f90
       call combine_Vlasov ! in IO.f90
       call PCDE4(Vlasov_ESim)
     else  ! for ion movable
       call Vlasov_ES(t,f,func,ntot)  ! in Vlasov_func.f90
       call combine_fluid  ! in IO.f90
       call combine_Vlasov ! in IO.f90
       call PCDE4(Vlasov_ES)
	 endif
!
   case ('EM')
     call Vlasov_EM(t,f,func,ntot)  ! in Vlasov_func.f90
     call combine_fluid  ! in IO.f90
     call combine_Vlasov ! in IO.f90
     call PCDE4(Vlasov_EM)     
!
   case default
     if (MyID == 0) then
     call Date_and_Time(dates,times,zone,values)
	 open (1,file='error.txt',status='unknown')
       write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                    values(5), ':', values(6)
       write (1,*) 'The model type of Vlasov equation has error !!'
       write (1,*) 'type =', atype
       write (1,*) 'Please check type !!'
       write (1,*) 'Program STOP at main program !!'
       close (1)
       write (6,11) values(1), '-', values(2), '-', values(3), ' ', &
                    values(5), ':', values(6)
       write (6,*) 'The model type of Vlasov equation has error !!'
       write (6,*) 'type =', atype
       write (6,*) 'Please check type !!'
       write (6,*) 'Program STOP at main program !!'
     endif
     call MPI_FINALIZE(IERR)
     STOP
   end select 
11 format (1x,i4,4(a1,i2.2))
!===============================================================================
   call MPI3D_FINALIZE
end program Vlasov_main
!===============================================================================
!===============================================================================
subroutine check_Vlasov
   use coef
   use wtime
   use MPI3D
   use comm_global
   use Vlasov
   use GRID6D
   use vector3D
   use Vlasov2fluid
   use CubicSpline   
   implicit double precision (A-H,O-Z)
   dimension :: divB(1), divB_global(1)
   dimension :: divE(1), divE_global(1)
   character*10 :: dates, times, zone
   integer, dimension(:) :: values(8)
   call Date_and_Time(dates,times,zone,values)
!-------------------------------------------------------------------------------
!  check the div B = 0
!
   call div3d_5th(ncx_mpi,ncy_mpi,ncz_mpi,f(ncBx),f(ncBy),f(ncBz),wkA, &
                  dx,dy,dz,ixp,iyp,izp)
!
!  find the biggest divergence B value
!
   wkA = dabs(wkA)
!
   divB(1) = maxval(wkA)
!   
!  find the maximum divergence B from the node in the cluster 
!
   call MPI_ALLREDUCE(divB(1),divB_global(1),1,MPI_REAL8,MPI_MAX,MPI_WORLD,IERR)
!
   if (MyID == 0) then
   	 if (it == 0) then 
   	 	 open (1,file='divB.txt',status='unknown')
   	 else
   	   open (1,file='divB.txt',status='unknown',position='append')
   	 endif
!
     write (1,*) it*dt, divB_global(1)
!     write (1,*) 't =', it*dt, 'max(divB) =', divB_global(1)
     close (1)
   endif
!-------------------------------------------------------------------------------
!  check distribution less than 0
!
   do ii = ncfe, ncfe_end
     !if (f(ii) < 0.d0) f(ii) = 0.d0
   enddo
!
   do ii = ncfi, ncfi_end
     if (f(ii) < 0.d0) f(ii) = 0.d0
   enddo
!-------------------------------------------------------------------------------
!  check the div E = Rhoi-Rhoe
!
   call div3d_5th(ncx_mpi,ncy_mpi,ncz_mpi,f(ncEx),f(ncEy),f(ncEz),wkdE, &
                  dx,dy,dz,ixp,iyp,izp)
   call density(ncx_mpi,ncy_mpi,ncz_mpi,nuex,nuey,nuez,f(ncfe), &
                Buex,Cuex,huex,Buey,Cuey,huey,Buez,Cuez,huez,Rhoe)
   call density(ncx_mpi,ncy_mpi,ncz_mpi,nuix,nuiy,nuiz,f(ncfi), &
                Buix,Cuix,huix,Buiy,Cuiy,huiy,Buiz,Cuiz,huiz,Rhoi)
!
    if (mtype == 'N') Rhoi = 1.d0
!-------------------------------------------------------------------------------
!  find the biggest divergence E value
!
   wkA = wkdE-Rhoi+Rhoe
   wkA = dabs(wkA)
   wkdE = dabs(wkdE)
!
   divE(1) = maxval(wkA)
!   
!  find the maximum divergence B from the node in the cluster 
!
   call MPI_ALLREDUCE(divE(1),divE_global(1),1,MPI_REAL8,MPI_MAX,MPI_WORLD,IERR)
!
   if (MyID == 0) then
   	 if (it == 0) then 
   	 	 open (1,file='divE-Rhoc.txt',status='unknown')
   	 else
   	   open (1,file='divE-Rhoc.txt',status='unknown',position='append')
   	 endif
!
!     write (1,*) 't =', it*dt, 'max(divE-Rhoc) =', divE_global(1)
     write (1,*) it*dt, divE_global(1)
     close (1)
   endif    	 
!-------------------------------------------------------------------------------
end subroutine check_Vlasov
!===============================================================================
!function Gaussian(u,V,Vth)     
!  pi     = 4*datan(1.d0)
!  s2pi   = dsqrt(2.d0*pi)
  !print *, 'pi is',pi
  !print *, 'v_th',Vth
!  Gaussian = dexp(-(u-V)**2/(2*Vth**2))/Vth/s2pi
!end function Gaussian
!===============================================================================
subroutine smooth_operator(f_t)
   use coef
   use wtime
   use MPI3D
   use comm_global
   use Vlasov
   use GRID6D
   use vector3D
   use Vlasov2fluid
   use CubicSpline   
   implicit double precision (A-H,O-Z)
   dimension :: f_t(1)
   real*8, allocatable, dimension(:) :: f_tt
   N = ncx_mpi*nuex
  
  pi     = 4*datan(1.d0)
  s2pi   = dsqrt(2.d0*pi)
   !print *, N
   !print *, 'v_b is ',v_b
   allocate(f_tt(N))
   f_tt = 0.d0
   do 10 i = 1,N
      f_tt(i) = f_t(i)
   10 continue
   !print *, 'ncx_mpi is',ncx_mpi
   !print *, 'nuex is ', nuex
   do 11 k = 1, ncz_mpi
   do 11 j = 1, ncy_mpi
   do 11 i = 1, ncx_mpi
     ijk  = i+(j-1)*ncx_mpi+(k-1)*ncx_mpi*ncy_mpi
     nuexy = nuex*nuey
     nuexyz = nuex*nuey*nuez
     do 12 ku = 1, nuez
     do 12 ju = 1, nuey
     do 12 iu = 1, nuex
       ijku  = iu+(ju-1)*nuex+(ku-1)*nuexy
       ii    = ijku+nuexyz*(ijk-1)  
       temp = 0.d0
       temp2 = 0.d0
       du = uex(2)-uex(1)
       do 13 iuu = 1,nuex
         ijkuu = iuu+(ju-1)*nuex+(ku-1)*nuexy
         iii   = ijkuu+nuexyz*(ijk-1)
         G = dexp(-(uex(iu)-uex(iuu))**2/(2*v_b**2))/v_b/s2pi
         temp2 = temp2 + G*du
         !print *, 'G    is ',G
         temp = temp + f_tt(iii)*G*du
         !print*, 'temp is ',temp
       13 continue
       f_t(ii) = temp/temp2
       !print*,'temp is',temp
	12 continue
	11 continue
   !print *, 'f_tt is',f_tt(10)
   !print *, 'f_t is ', f_t(10)
    deallocate(f_tt)
end subroutine smooth_operator
