!===============================================================================
!===============================================================================
module coef

     real*8 :: relerr, rerr ! a nearly negligible number relative to 1
     real*8 :: pi     ! pi = 3.141592....
!-------------------------------------------------------------------------------
!  the setting of the predictor-corrector method
!
     integer :: irel, iter
     character(1) :: ccheck ! 'Y' : check convergence ; 'N' : doesn't check convergence
     integer :: iOMP        ! number of OpenMP threads per MPI node
   contains
!-------------------------------------------------------------------------------
     subroutine const
       implicit double precision (a-h,o-z)
       pi = 4*datan(1.d0)
       relerr = irel*epsilon(1.d0)
       rerr = 20*epsilon(1.d0)
     end subroutine const
!-------------------------------------------------------------------------------
     function DIFAB(A,B)
!
!  calculate DIFAB = A-B
!
       implicit double precision (a-h,o-z)
       DIFAB = A-B
!       return
       err  = max(dabs(A),dabs(B))*rerr
       if (dabs(DIFAB) < err) DIFAB = 0.d0
     end function DIFAB
!-------------------------------------------------------------------------------
     function sumAB(A,B)
!
!  calculate sumAB = A+B
!
       implicit double precision (a-h,o-z)
       sumAB = DIFAB(A,-B)
     end function sumAB
!------------------------------------------------------------------------------- 
     function sumABC(A,B,C)
!
!  calculate sum = A+B+C
!
       implicit double precision (a-h,o-z)
       sumABC = DIFAB(sumAB(A,B),-C)
     end function sumABC
!-------------------------------------------------------------------------------
     function sumABCD(A,B,C,D)
!
!  calculate sum = A+B+C+D
!
       implicit double precision (a-h,o-z)
       sumABCD = DIFAB(sumABC(A,B,C),-D)
     end function sumABCD
!-------------------------------------------------------------------------------
     function sumABCDE(A,B,C,D,E)
!
!  calculate sum = A+B+C+D+E
!
       implicit double precision (a-h,o-z)
       sum1 = sumABC(A,B,C)
       sum2 = sumAB(D,E)
       sumABCDE = sumAB(sum1,sum2)
     end function sumABCDE
!-------------------------------------------------------------------------------    
end module coef
!===============================================================================
!===============================================================================
module wtime
   integer :: it, it_plot_f, it_plot_d,it_plot_s, it_stop
   real*8  :: t, dt
end module wtime
!===============================================================================
!===============================================================================
module comm_global
     use wtime
     use coef
!-------------------------------------------------------------------------------
!  the dimension in real space x, y, z for global index
!
     integer :: ncx, ncy, ncz
     integer :: ncxy, ncxyz !, ncyz, ncxz
!-------------------------------------------------------------------------------
!  the dimension in electron velocity space uex, uey, uez for global index
     integer :: nuex, nuey, nuez
     integer :: nuexy, nuexyz !, nueyz, nuexz
!-------------------------------------------------------------------------------
!  the dimension in ion velocity space uix, uiy, uiz for global index
     integer :: nuix, nuiy, nuiz
     integer :: nuixy, nuixyz !, nuiyz, nuixz
!-------------------------------------------------------------------------------
!  P*Q*R = the number of nodes
     integer :: P    ! the number of nodes in x direction
     integer :: Q    ! the number of nodes in y direction
     integer :: R    ! the number of nodes in z direction
!-------------------------------------------------------------------------------
!  Boundary condition along x-, y-, z-direction
!    ixp, iyp, izp  = 1 : periodic boundary condition
!    ixp, iyp, izp \= 1 : free boundary condition
!
     integer :: ixp, iyp, izp
!-------------------------------------------------------------------------------
!  proton-to-electron mass ratio
     real*8 :: ami
!-------------------------------------------------------------------------------
!  light speed
     real*8 :: C
!-------------------------------------------------------------------------------
!  gaussian smooth width
     real*8 :: v_b
!-------------------------------------------------------------------------------
!  'ES' for electrostatic, 'EM' for electromagnetic
!   model type for the Vlasov equation 
!
     character(2) :: atype
!-------------------------------------------------------------------------------
!  'Y' for ion movable, 'N' for ion immovable
!   model type for the Vlasov equation 
!
     character(1) :: mtype
!-------------------------------------------------------------------------------
!  'Y' for relativistic, 'N' for non-relativistic
!   model type for the Vlasov equation 
!
     character(1) :: rtype
!-------------------------------------------------------------------------------
!  output data is in 'single' or 'double' precision format
!
!     character(6) :: data_format
!-------------------------------------------------------------------------------
!  'Y' : check CFL ; 'N' don't check CFL
!
!     character(1) :: CFLcheck
!-------------------------------------------------------------------------------
!  adiabatic index
!
!     real*8 :: gamma
!-------------------------------------------------------------------------------
   contains
!===============================================================================
   subroutine get_para
     implicit double precision (A-H,O-Z)
!
!  get the parameters from init_para.txt
!
     open (1,file='init_para.txt',status='old')
!     read (1,*) gamma ! adiabatic index
!
! Boundary condition along x,- y-, z-direction
!   ixp, iyp, izp  = 1 : periodic boundary condition
!   ixp, iyp, izp \= 1 : free boundary condition
!
     read (1,*) ixp ! the boundary condition for x-direction
     read (1,*) iyp ! the boundary condition for y-direction
     read (1,*) izp ! the boundary condition for z-direction
!
!  P*Q*R = the number of nodes
!
     read (1,*) P   ! P = the number of nodes in x direction
     read (1,*) Q   ! Q = the number of nodes in y direction
     read (1,*) R   ! R = the number of nodes in z direction
!
     read (1,*) it_plot_d ! the time step of the data output for the distribution
     read (1,*) it_plot_f ! the time step of the data output for the fluid data
     read (1,*) it_plot_s ! the time step of the smoothing of distribution
     read (1,*) it_stop ! the end of the iteration times in the simulation run
!
     read (1,*) irel ! convergence condition
     read (1,*) iter ! iteration times for the Predictor-Corrector method
     read (1,'(a1)') ccheck ! 'Y' : check convergence ; 'N' : doesn't check convergence
!
! the model type of the Vlasov equation 'ES' for electrostatic, 'EM' for electromagnetic
!
     read (1,'(a2)') atype
!
! the model type of the Vlasov equation 'Y' for ion movable, 'N' for ion immovable
!	 
     read (1,'(a1)') mtype
!
! the model type of the Vlasov equation 'Y' for relativistic, 'N' for non-relativistic
!	 
     read (1,'(a1)') rtype
!
!     read (1,'(a6)') data_format  ! output data is in 'single' or 'double' precision format
!     read (1,'(a1)') CFLcheck ! 'Y' : check CFL ; 'N' don't check CFL
     read (1,*) iOMP ! number of OpenMP threads per MPI node
!
     close (1)
!
!  read the mesh number along the x-, y-, z- axis, mass ratio, and light speed 
!
     open (2,file='init.bin',status='old',form='unformatted')
     read (2) ncx, ncy, ncz, nuex, nuey, nuez, nuix, nuiy, nuiz, ami, C
     close (2)
     open (3,file='gauss.bin',status='old',form='unformatted')
     read (3) v_b
     close (3)
   end subroutine get_para
!===============================================================================     
end module comm_global
!===============================================================================
!===============================================================================
!   prepare the parameters for MPI
!
include 'MPI3D.f90'
!===============================================================================
!===============================================================================
module Vlasov
     use comm_global
     use MPI3D
!-------------------------------------------------------------------------------
!  global array
!
     real*8, allocatable, dimension(:) :: fe_global, fi_global
     real*8, allocatable, dimension(:) :: A_global
     real*8, allocatable, dimension(:) :: x_global, y_global, z_global
     real*8, allocatable, dimension(:) :: fuex_global, fuey_global, fuez_global
     real*8, allocatable, dimension(:) :: fuix_global, fuiy_global, fuiz_global
!-------------------------------------------------------------------------------
! 6D array in one dimension for local memory
!
     real*8, allocatable, dimension(:) :: f, func, func1, func2, func3
     real*8, allocatable, dimension(:) :: func11, func12, func13
!
! working array for local memory
!
     real*8, allocatable, dimension(:) :: wkf1, wkf2       
     real*8, allocatable, dimension(:) :: wkdAx, wkdAy, wkdAz, wkA, wkdE
     real*8, allocatable, dimension(:) :: wkddA
     real*8, allocatable, dimension(:) :: wkEx
     real*8, allocatable, dimension(:) :: wkfex, wkfey, wkfez, wkfexv
     real*8, allocatable, dimension(:) :: wkfix, wkfiy, wkfiz
     real*8, allocatable, dimension(:) :: wkJex, wkJey, wkJez
     real*8, allocatable, dimension(:) :: wkJix, wkJiy, wkJiz
     real*8, allocatable, dimension(:) :: wkJx0, wkJy0, wkJz0
     real*8, allocatable, dimension(:) :: wkfepx , wkfepy , wkfepz
     real*8, allocatable, dimension(:) :: wkfeppx, wkfeppy, wkfeppz
     real*8, allocatable, dimension(:) :: wkfipx , wkfipy , wkfipz
     real*8, allocatable, dimension(:) :: wkfippx, wkfippy, wkfippz
!
!  local arrays for output data
!
     real*8, allocatable, dimension(:) :: Rhoe, Rhoi
     real*8, allocatable, dimension(:) :: avguex, avguey, avguez
     real*8, allocatable, dimension(:) :: avgvex, avgvey, avgvez
     real*8, allocatable, dimension(:) :: avguix, avguiy, avguiz
     real*8, allocatable, dimension(:) :: avgvix, avgviy, avgviz
     real*8, allocatable, dimension(:) :: Pex, Pey, Pez
     real*8, allocatable, dimension(:) :: Pix, Piy, Piz
     real*8, allocatable, dimension(:) :: fuex, fuey, fuez
     real*8, allocatable, dimension(:) :: fuix, fuiy, fuiz
!
! diffusion coefficient for local memory
!
     real*8, allocatable, dimension(:) :: etaB, etafe, etafi
!
! electron and ion velocity space
!
     real*8, allocatable, dimension(:) :: uex, uey, uez
     real*8, allocatable, dimension(:) :: uix, uiy, uiz
!
!  local dimension and index
!
     integer :: ntot, nofe, nofi
     integer :: ncfe, ncfe_end, ncfi, ncfi_end
     integer :: ncBx, ncBx_end, ncBy, ncBy_end, ncBz, ncBz_end
     integer :: ncEx, ncEx_end, ncEy, ncEy_end, ncEz, ncEz_end
   contains
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  allocate global arrays in rank 0
!
   subroutine setVlasov_global
     ncxyz  = ncx*ncy*ncz
     nuexy  = nuex*nuey
     nuexyz = nuex*nuey*nuez
     nuixy  = nuix*nuiy
     nuixyz = nuix*nuiy*nuiz
!
     allocate (uex(nuex), uey(nuey), uez(nuez))
     allocate (uix(nuix), uiy(nuiy), uiz(nuiz))
!
     if (MyID == 0) then
       allocate (x_global(ncx), y_global(ncy), z_global(ncz))
     else
       allocate (x_global(1), y_global(1), z_global(1))
     endif 
   end subroutine setVlasov_global
!-------------------------------------------------------------------------------
!  allocate arrays dimension used in Vlasov
!
   subroutine setVlasov
!-------------------------------------------------------------------------------
!  Let f = fe(ncfe:ncfe_end)+fi(ncfi:ncfi_end)
!         +Bx(ncBx:ncBx_end)+By(ncBy:ncBy_end)+Bz(ncBz:ncBz_end)
!         +Ex(ncEx:ncEx_end)+Ey(ncEy:ncEy_end)+Ez(ncEz:ncEz_end)
!-------------------------------------------------------------------------------
!  calculate the parameter for the local array
!
     nofe     = ncxyz_mpi*nuexyz
     nofi     = ncxyz_mpi*nuixyz
!
     ntot     = nofe+nofi+6*ncxyz_mpi
     ncfe     = 1
     ncfe_end = nofe
     ncfi     = ncfe_end+1
     ncfi_end = ncfe_end+nofi
     ncBx     = ncfi_end+1
     ncBx_end = ncfi_end+ncxyz_mpi
     ncBy     = ncBx_end+1
     ncBy_end = ncBx_end+ncxyz_mpi     
     ncBz     = ncBy_end+1
     ncBz_end = ncBy_end+ncxyz_mpi
     ncEx     = ncBz_end+1
     ncEx_end = ncBz_end+ncxyz_mpi
     ncEy     = ncEx_end+1
     ncEy_end = ncEx_end+ncxyz_mpi
     ncEz     = ncEy_end+1
     ncEz_end = ncEy_end+ncxyz_mpi
!
     allocate (f(ntot), func(ntot), func1(ntot), func2(ntot), func3(ntot)) 
     allocate (func11(ntot), func12(ntot), func13(ntot)) 
     allocate (wkf1(ntot), wkf2(ntot))  ! working array for time integration
     allocate (wkA(ncxyz_mpi), wkddA(ncxyz_mpi), wkdE(ncxyz_mpi))
     allocate (wkdAx(ncxyz_mpi), wkdAy(ncxyz_mpi), wkdAz(ncxyz_mpi)) 
     allocate (wkJex(ncxyz_mpi), wkJey(ncxyz_mpi),  wkJez(ncxyz_mpi))
     allocate (wkJix(ncxyz_mpi), wkJiy(ncxyz_mpi),  wkJiz(ncxyz_mpi))
     allocate (wkJx0(ncxyz_mpi), wkJy0(ncxyz_mpi),  wkJz0(ncxyz_mpi))
     allocate (etaB(ncxyz_mpi), etafe(ncxyz_mpi), etafi(ncxyz_mpi))
     allocate (Rhoe(ncxyz_mpi), Rhoi(ncxyz_mpi))
     allocate (avguex(ncxyz_mpi), avguey(ncxyz_mpi), avguez(ncxyz_mpi))
     allocate (avgvex(ncxyz_mpi), avgvey(ncxyz_mpi), avgvez(ncxyz_mpi))
     allocate (avguix(ncxyz_mpi), avguiy(ncxyz_mpi), avguiz(ncxyz_mpi))
     allocate (avgvix(ncxyz_mpi), avgviy(ncxyz_mpi), avgviz(ncxyz_mpi))
     allocate (Pex(ncxyz_mpi), Pey(ncxyz_mpi), Pez(ncxyz_mpi))
     allocate (Pix(ncxyz_mpi), Piy(ncxyz_mpi), Piz(ncxyz_mpi))
     allocate (fuex(nuex*ncxyz_mpi),fuey(nuey*ncxyz_mpi),fuez(nuez*ncxyz_mpi))
     allocate (fuix(nuix*ncxyz_mpi),fuiy(nuiy*ncxyz_mpi),fuiz(nuiz*ncxyz_mpi))
   end subroutine setVlasov
!===============================================================================
end module Vlasov
!===============================================================================
!===============================================================================
!  5th-order iPCCFD methods 
!    for General-Purpose PDE solver
!
include 'CenDif3D_MPI.f90'
!===============================================================================
!  calculate first order derivative and integration in velocity space 
!    using cubic spline method
!
include 'CubicSpline.f90'
!===============================================================================
module GRID6D
!
!  Arrays for spatial grids and 
!    pre-setting arrays for the General-Purpose PDE solver
!    using Finite Difference Method and Cubic spline method
!
     use MPI3D
     use comm_global
     use Vlasov
	 use CubicSpline
!
!  allocate memories
!
     real*8, allocatable, dimension(:) :: x, y, z
     real*8, allocatable, dimension(:) :: vex, vey, vez
     real*8, allocatable, dimension(:) :: vix, viy, viz
     real*8, allocatable, dimension(:) :: uex3D, uey3D, uez3D
     real*8, allocatable, dimension(:) :: uix3D, uiy3D, uiz3D
!
!  allocate memories for cubic spline
!
     real*8, allocatable, dimension(:) :: huex, Buex, Cuex
	 real*8, allocatable, dimension(:) :: huey, Buey, Cuey
	 real*8, allocatable, dimension(:) :: huez, Buez, Cuez
     real*8, allocatable, dimension(:) :: huix, Buix, Cuix
	 real*8, allocatable, dimension(:) :: huiy, Buiy, Cuiy
	 real*8, allocatable, dimension(:) :: huiz, Buiz, Cuiz
!
     real*8 :: dx, dy, dz
!
   contains
!===============================================================================
!  allocate arrays for phase space in Vlasov simulation
!-------------------------------------------------------------------------------
   subroutine alloc_3D3V
     allocate (x(ncx_mpi), y(ncy_mpi), z(ncz_mpi))
     allocate (vex(nuexyz), vey(nuexyz), vez(nuexyz))
     allocate (vix(nuixyz), viy(nuixyz), viz(nuixyz))
     allocate (uex3D(nuexyz), uey3D(nuexyz), uez3D(nuexyz))
     allocate (uix3D(nuixyz), uiy3D(nuixyz), uiz3D(nuixyz))
	 allocate (huex(nuex), Buex(nuex), Cuex(nuex))
	 allocate (huey(nuey), Buey(nuey), Cuey(nuey))
	 allocate (huez(nuez), Buez(nuez), Cuez(nuez))
	 allocate (huix(nuix), Buix(nuix), Cuix(nuix))
	 allocate (huiy(nuiy), Buiy(nuiy), Cuiy(nuiy))
	 allocate (huiz(nuiz), Buiz(nuiz), Cuiz(nuiz))
!
     if (atype == 'EM') then
       call u2v_RE
       goto 555
     endif
!
     select case (rtype)
     case ('Y')
       call u2v_RE
     case ('N')
       call u2v_NRE
     case default
       if (MyID == 0) then
         open (1,file='error.txt',status='unknown')
         write (1,*) 'The seting in initial condition has error !!'
         write (1,*) "rtpye =", rtype, "relativistic parameter must be equal to 'Y' or 'N'"
         write (1,*) 'Please check rtype !!'
         write (1,*) 'Program STOP at alloc_3D3V !!'
         close (1)
         write (6,*) 'The seting in initial condition has error !!'
         write (6,*) "rtpye =", rtype, "relativistic parameter must be equal to 'Y' or 'N'"
         write (6,*) 'Please check rtype !!'
         write (6,*) 'Program STOP at alloc_3D3V !!'
       endif
       call MPI_FINALIZE(IERR)
       STOP
     end select 
!
555  continue
     call CStable
   end subroutine alloc_3D3V
!===============================================================================
!  prepare cubic spline table
!-------------------------------------------------------------------------------
   subroutine CStable
     call geth(uex,huex,nuex)
     call FCSYP0(huex,Buex,Cuex,nuex)
!
     if (nuey == 1) then
       huey = 0.d0
       Buey = 0.d0
       Cuey = 0.d0
     else
       call geth(uey,huey,nuey)
       call FCSYP0(huey,Buey,Cuey,nuey)
     endif
!
     if (nuez == 1) then
       huez = 0.d0
       Buez = 0.d0
       Cuez = 0.d0
     else
       call geth(uez,huez,nuez)
       call FCSYP0(huez,Buez,Cuez,nuez)
     endif
!
     call geth(uix,huix,nuix)
     call FCSYP0(huix,Buix,Cuix,nuix)
!
     if (nuiy == 1) then
       huiy = 0.d0
       Buiy = 0.d0
       Cuiy = 0.d0
     else
       call geth(uiy,huiy,nuiy)
       call FCSYP0(huiy,Buiy,Cuiy,nuiy)
     endif
!
     if (nuiz == 1) then
       huiz = 0.d0
       Buiz = 0.d0
       Cuiz = 0.d0
     else
       call geth(uiz,huiz,nuiz)
       call FCSYP0(huiz,Buiz,Cuiz,nuiz)
     endif
   end subroutine CStable
!===============================================================================
!  calculate v from u
!-------------------------------------------------------------------------------
!  for non-relativistic
!
   subroutine u2v_NRE
     implicit double precision (a-h,o-z)
!
!  uex to vex, uey to vey, uez to vez
!
     do ku = 1, nuez
     do ju = 1, nuey
     do iu = 1, nuex
       ijku = iu+(ju-1)*nuex+(ku-1)*nuexy
       vex(ijku)   = uex(iu)
       vey(ijku)   = uey(ju)
       vez(ijku)   = uez(ku)
       uex3D(ijku) = uex(iu)
       uey3D(ijku) = uey(ju)
       uez3D(ijku) = uez(ku) 
     enddo
     enddo
     enddo
!
!  uix to vix, uiy to viy, uiz to viz
!
     do ku = 1, nuiz
     do ju = 1, nuiy
     do iu = 1, nuix
       ijku = iu+(ju-1)*nuix+(ku-1)*nuixy
       vix(ijku)   = uix(iu)
       viy(ijku)   = uiy(ju)
       viz(ijku)   = uiz(ku)
       uix3D(ijku) = uix(iu)
       uiy3D(ijku) = uiy(ju)
       uiz3D(ijku) = uiz(ku)
     enddo
     enddo
     enddo     
   end subroutine u2v_NRE
!-------------------------------------------------------------------------------
!  for relativistic
!
   subroutine u2v_RE
     implicit double precision (a-h,o-z)
!
!  uex to vex, uey to vey, uez to vez
!
     do ku = 1, nuez
     do ju = 1, nuey
     do iu = 1, nuex
       ijku = iu+(ju-1)*nuex+(ku-1)*nuexy
       uu = uex(iu)**2+uey(ju)**2+uez(ku)**2
       us = dsqrt(1.d0+uu/C**2)
       vex(ijku) = uex(iu)/us
       vey(ijku) = uey(ju)/us
       vez(ijku) = uez(ku)/us
       uex3D(ijku) = uex(iu)
       uey3D(ijku) = uey(ju)
       uez3D(ijku) = uez(ku) 
     enddo
     enddo
     enddo
!
!  uix to vix, uiy to viy, uiz to viz
!
     do ku = 1, nuiz
     do ju = 1, nuiy
     do iu = 1, nuix
       ijku = iu+(ju-1)*nuix+(ku-1)*nuixy
       uu = uix(iu)**2+uiy(ju)**2+uiz(ku)**2
       us = dsqrt(1.d0+uu/C**2)
       vix(ijku) = uix(iu)/us
       viy(ijku) = uiy(ju)/us
       viz(ijku) = uiz(ku)/us
       uix3D(ijku) = uix(iu)
       uiy3D(ijku) = uiy(ju)
       uiz3D(ijku) = uiz(ku)
     enddo
     enddo
     enddo     
   end subroutine u2v_RE
!===============================================================================
end module GRID6D
!===============================================================================
