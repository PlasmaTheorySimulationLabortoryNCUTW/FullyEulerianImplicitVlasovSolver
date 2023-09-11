!===============================================================================
module coef
     real(8), parameter :: pi     = 4*datan(1.d0)
     real(8), parameter :: s2pi   = dsqrt(2.d0*pi)
     real(8), parameter :: ami    = 1836.d0
     real(8), parameter :: smalle  = 1.d-5
     real(8), parameter :: smalli  = 1.d-3
     real(8), parameter :: factor = 1.d0
     real*8 :: rerr
   contains
!-------------------------------------------------------------------------------
     subroutine const
       implicit double precision (a-h,o-z)
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
end module coef
!===============================================================================
include 'CubicSpline.f90'
include 'Vlasov2fluid.f90'
!===============================================================================
module find_max
     use coef
     implicit double precision (a-h,o-z)
   contains
!-------------------------------------------------------------------------------
   function umax(du_ini,Rho,V0,Vth,small)
     du = du_ini
     uu = V0
     ff = Rho*Gaussian(uu,V0,Vth)
     do while (ff > small)
       uu = uu+du
       ff = Rho*Gaussian(uu,V0,Vth)
     enddo
     u_plus = dabs(factor*uu)
!
     du = -du_ini
     uu = V0
     ff = Rho*Gaussian(uu,V0,Vth)
     do while (ff > small)
       uu = uu+du
       ff = Rho*Gaussian(uu,V0,Vth)
     enddo
     u_minus = dabs(factor*uu)
!     print *, u_plus,u_minus
!
     umax = max(u_plus,u_minus)
   end function umax 
!-------------------------------------------------------------------------------
   function Gaussian(u,V,Vth)     
     Gaussian = dexp(-(u-V)**2/(2*Vth**2))/Vth/s2pi
   end function Gaussian
!-------------------------------------------------------------------------------   
end module find_max   
!===============================================================================
program initial_Vlasov_ES1D1V
!===============================================================================
!  initialize the initial condition of the Vlasov simulation code
!  This is the simulation for the Landau damping
!
   use coef
   use find_max
   use omp_lib
!   use integral
   use Vlasov2fluid
   use CubicSpline
   implicit double precision (a-h,o-z)
   parameter (ncx = 65, ncy = 1, ncz = 1)
   parameter (nuex = 64, nuey = 1, nuez = 1)
   parameter (nuix = 64, nuiy = 1, nuiz = 1)
   parameter (ncxy = ncx*ncy, ncyz = ncy*ncz, ncxz = ncx*ncz)
   parameter (ncxyz = ncx*ncy*ncz)
   parameter (nuexy = nuex*nuey, nueyz = nuey*nuez)
   parameter (nuexyz = nuex*nuey*nuez)
   parameter (nuixy = nuix*nuiy, nuiyz = nuiy*nuiz)
   parameter (nuixyz = nuix*nuiy*nuiz)
   parameter (nofe = nuexyz*ncxyz, nofi = nuixyz*ncxyz)
   parameter (ntot = nofe+nofi+6*ncxyz)
   real(8), allocatable, dimension(:) :: x, y, z, hx
   real(8), allocatable, dimension(:) :: uex, uey, uez
   real(8), allocatable, dimension(:) :: uix, uiy, uiz
   real(8), allocatable, dimension(:) :: fe, fi
   real(8), allocatable, dimension(:) :: Bx, By, Bz
   real(8), allocatable, dimension(:) :: Ex, Ey, Ez
   real(8), allocatable, dimension(:) :: etaB, etafe, etafi
   real(8), allocatable, dimension(:) :: wkRhoe, wkVe, wkPe, wkfe, wkJex  ! for testing
   real(8), allocatable, dimension(:) :: wkRhoi, wkVi, wkPi, wkfi, wkJix  ! for testing
   real(8), allocatable, dimension(:) :: Buex, Cuex, huex
   real(8), allocatable, dimension(:) :: Buey, Cuey, huey
   real(8), allocatable, dimension(:) :: Buez, Cuez, huez
   real(8), allocatable, dimension(:) :: Buix, Cuix, huix
   real(8), allocatable, dimension(:) :: Buiy, Cuiy, huiy
   real(8), allocatable, dimension(:) :: Buiz, Cuiz, huiz
   dimension :: tanhX(ncx)
   call const
!-------------------------------------------------------------------------------
!   iOMP = OMP_get_max_threads()
   iOMP = 4
   print *, iOMP
!   stop
   call OMP_SET_NUM_THREADS(iOMP)
   allocate (x(ncx), y(ncy), z(ncz), hx(ncx))
   allocate (uex(nuex), uey(nuey), uez(nuez))
   allocate (uix(nuix), uiy(nuiy), uiz(nuiz))
   allocate (Buex(nuex), Cuex(nuex), huex(nuex))
   allocate (Buey(nuey), Cuey(nuey), huey(nuey))
   allocate (Buez(nuez), Cuez(nuez), huez(nuez))
   allocate (Buix(nuix), Cuix(nuix), huix(nuix))
   allocate (Buiy(nuiy), Cuiy(nuiy), huiy(nuiy))
   allocate (Buiz(nuiz), Cuiz(nuiz), huiz(nuiz))
!===============================================================================
!  set the physical quantities of the background equilibrium state
!
   ak = 0.5d0
   alpha = 0.001d0
!
   C = 200.d0
!   CFL = 0.1d0
!   dt = 0.025d0 !CFL*(1/C)
!   print *, dt 	 
!
   Vx0 = 0.d0
!   Vx1 = 0.d0
!   Vx2 = 0.d0
!   Vy0 = 0.d0
!   Vz0 = 0.d0
!   Vy1 = Vy0
!   Vy2 = Vy0
!   Vz1 = Vz0
!   Vz2 = Vz0
!
   Pe0 = 1.d0
   Pi0 = Pe0
!   Pe1 = Pe0
!   Pi1 = Pe0
!   Pi2 = Pe0
!   P1  = Pe1+pi1
!   Pe2 = P1-Pi2
!
   Rho0 = 1.d0
   Rho1 = Rho0
   Rho2 = Rho0
!
   Vthe0 = dsqrt(Pe0/Rho0)
   Vthi0 = dsqrt(Pi0/Rho0/ami)
!   Vthe1 = dsqrt(Pe1/Rho1)
!   Vthe2 = dsqrt(Pe2/Rho2)
!   Vthi1 = dsqrt(Pi1/Rho1/ami)
!   Vthi2 = dsqrt(Pi2/Rho2/ami)
!
!   beta = 2.d0
!   B1  = dsqrt(2.d0*P1/beta)
!   B2  = B1
!   thetaB1 = 0.d0  !(in degree)
!   thetaB2 = thetaB1
!   Bx1 = B1*dcos(thetaB1/180.d0*pi)
!   By1 = B1*dsin(thetaB1/180.d0*pi)
!   Bz1 = 0.d0
!   Bx2 = B2*dcos(thetaB2/180.d0*pi)
!   By2 = B2*dsin(thetaB2/180.d0*pi)
!   Bz2 = 0.d0
!===============================================================================
!  set the mesh along the x-axis direction
!
   x_start = 0.d0
   x_end   = 4*pi
   delx = (x_end-x_start)/dfloat(ncx-1)
! 
   do i = 1, ncx
     x(i) = x_start+(i-1)*delx
   enddo
!-------------------------------------------------------------------------------
!  set the mesh along the y-axis direction
!
   y_start = -0.d0
   y_end   = +0.d0
   dely = 0.d0 !(y_end-y_start)/dfloat(ncy-1)
! 
   do j = 1, ncy
     y(j) = y_start+(j-1)*dely
   enddo
!-------------------------------------------------------------------------------
!  set the mesh along the z-axis direction
!
   z_start = -0.d0
   z_end   = +0.d0
   delz = 0.d0 !(z_end-z_start)/dfloat(ncz-1)
! 
   do k = 1, ncz
     z(k) = z_start+(k-1)*delz
   enddo
!===============================================================================
!  set the mesh along the uex-axis direction
!
!   uex_max1 = umax(0.1d0,Rho1,Vx0,Vthe0,smalle)
!   uex_max2 = umax(0.1d0,Rho2,Vx2,Vthe2,smalle)
   uex_max  = 6.d0 !max(uex_max1,uex_max2)
!   print *, 'uex_max', uex_max1, uex_max2
!
   uex_start = -uex_max
   uex_end   = +uex_max
   duex = (uex_end-uex_start)/dfloat(nuex-1)
! 
   do iu = 1, nuex
     uex(iu) = uex_start+(iu-1)*duex
!     print *, iu, uex(iu)
   enddo
!   stop
!-------------------------------------------------------------------------------
!  set the mesh along the uey-axis direction
!
!   uey_max1 = umax(0.1d0,Rho1,Vy1,Vthe1,smalle)
!   uey_max2 = umax(0.1d0,Rho2,Vy2,Vthe2,smalle)
!   uey_max  = max(uey_max1,uey_max2)
!   print *, 'uey_max', uey_max1, uey_max2
!
!   uey_start = -uey_max
!   uey_end   = +uey_max
!   duey = (uey_end-uey_start)/dfloat(nuey-1)
! 
!   do ju = 1, nuey
!     uey(ju) = uey_start+(ju-1)*duey
!     print *, ju, uey(ju)
!   enddo
!   stop
!-------------------------------------------------------------------------------
!  set the mesh along the uez-axis direction
!
!   uez_max1 = umax(0.1d0,Rho1,Vz1,Vthe1,smalle)
!   uez_max2 = umax(0.1d0,Rho2,Vz2,Vthe2,smalle)
!   uez_max  = max(uez_max1,uez_max2)
!   print *, 'uez_max', uez_max1, uez_max2
!
!   uez_start = -uez_max
!   uez_end   = +uez_max
!   duez = (uez_end-uez_start)/dfloat(nuez-1)
! 
!   do ku = 1, nuez
!     uez(ku) = uez_start+(ku-1)*duez
!     print *, ku, uez(ku)
!   enddo
!   stop
!===============================================================================
!  set the mesh along the uix-axis direction
!
   uix_max1 = umax(0.001d0,Rho0,Vx0,Vthi0,smalli)
!   uix_max2 = umax(0.001d0,Rho2,Vx2,Vthi2,smalli)
   uix_max  = max(uix_max1,uix_max2)
!   print *, 'uix_max', uix_max1, uix_max2
!
   uix_start = -uix_max
   uix_end   = +uix_max
   duix = (uix_end-uix_start)/dfloat(nuix-1)
! 
   do iu = 1, nuix
     uix(iu) = uix_start+(iu-1)*duix
!     print *, iu, uix(iu)
   enddo
!   stop
!-------------------------------------------------------------------------------
!  set the mesh along the uiy-axis direction
!
!   uiy_max1 = umax(0.001d0,Rho1,Vy1,Vthi1,smalli)
!   uiy_max2 = umax(0.001d0,Rho2,Vy2,Vthi2,smalli)
!   uiy_max  = max(uiy_max1,uiy_max2)
!   print *, 'uiy_max', uiy_max1, uiy_max2
!
!   uiy_start = -uiy_max
!   uiy_end   = +uiy_max
!   duiy = (uiy_end-uiy_start)/dfloat(nuiy-1)
! 
!   do ju = 1, nuiy
!     uiy(ju) = uiy_start+(ju-1)*duiy
!     print *, ju, uiy(ju)
!   enddo
!   stop
!-------------------------------------------------------------------------------
!  set the mesh along the uiz-axis direction
!
!   uiz_max1 = umax(0.001d0,Rho1,Vz1,Vthi1,smalli)
!   uiz_max2 = umax(0.001d0,Rho2,Vz2,Vthi2,smalli)
!   uiz_max  = max(uiz_max1,uiz_max2)
!   print *, 'uiz_max', uiz_max1, uiz_max2
!
!   uiz_start = -uiz_max
!   uiz_end   = +uiz_max
!   duiz = (uiz_end-uiz_start)/dfloat(nuiz-1)
! 
!   do ku = 1, nuiz
!     uiz(ku) = uiz_start+(ku-1)*duiz
!     print *, ku, uiz(ku)
!   enddo
!   stop
!===============================================================================
!   open (3,file='axis.bin',status='unknown',form='binary')
   open (3,file='axis.bin',status='unknown',form='unformatted')
   write (3) ncx, ncy, ncz, nuex, nuey, nuez, nuix, nuiy, nuiz, C
   write (3) x, y, z, uex, uey, uez, uix, uiy, uiz 
   close (3)
   print *, 'save axis'
!===============================================================================
   dx = x(2)-x(1)
!   CFL = 0.5 ! dt = 0.0423d0
!   dt = CFL*dx/uex(nuex)
   dt = 0.005d0
   print *, 'dt = ', dt
   print *, 'dx = ', dx
!   print *, 'CFLx = ', (uex(nuex)*dt)/dx
!   stop
!===============================================================================
!  prepare cubic spline table
!
   call geth(uex,huex,nuex)
!   call geth(uey,huey,nuey)
!   call geth(uez,huez,nuez)
   call FCSYP0(huex,Buex,Cuex,nuex)
!   call FCSYP0(huey,Buey,Cuey,nuey)
!   call FCSYP0(huez,Buez,Cuez,nuez)
!
   call geth(uix,huix,nuix)
!   call geth(uiy,huiy,nuiy)
!   call geth(uiz,huiz,nuiz)
   call FCSYP0(huix,Buix,Cuix,nuix)
!   call FCSYP0(huiy,Buiy,Cuiy,nuiy)
!   call FCSYP0(huiz,Buiz,Cuiz,nuiz)
!===============================================================================
!   Bxsum = Bx1+Bx2
!   Bxdif = Bx2-Bx1   
!   Bysum = By1+By2
!   Bydif = By2-By1
!   Bzsum = Bz1+Bz2
!   Bzdif = Bz2-Bz1
!
!   Rhosum = Rho1+Rho2
!   Rhodif = Rho2-Rho1
!
!   ADe1   = dsqrt(Rho1**3/Pe1)
!   ADe2   = dsqrt(Rho2**3/Pe2)
!   ADesum = ADe1+ADe2
!   ADedif = ADe2-ADe1  
!
!   ADi1   = dsqrt(Rho1**3/Pi1)
!   ADi2   = dsqrt(Rho2**3/Pi2)
!   ADisum = ADi1+ADi2
!   ADidif = ADi2-ADi1
!
!   W0   = 5.d0
!   X0   = 0.d0
!   Y0   = 0.d0
!   Z0   = 0.d0
!===============================================================================
!   do i = 1, ncx     
!     tanhX(i) = dtanh((x(i)-X0)/W0)
!     write (6,*) x(i), tanhX(i)
!   enddo
!-------------------------------------------------------------------------------
   open (1,file='init.bin',status='unknown',form='unformatted')
   write (1) ncx, ncy, ncz, nuex, nuey, nuez, nuix, nuiy, nuiz, ami, C, 0.d0, 0, dt
!===============================================================================
!  get the initial condition of fe
!
   allocate (fe(nofe))
!$OMP parallel do private(ii,ijk,ijku,Rho,Vx,Ade,Pe,Vthe,wkfex,wkfey,wkfez) collapse(3)
   do 11 k = 1, ncz
   do 11 j = 1, ncy
   do 11 i = 1, ncx
!
     ijk  = i+(j-1)*ncx+(k-1)*ncxy
     Rho  = Rho0
     Vx   = Vx0
!     Pe   = Pe0
     Vthe = 1.d0 !dsqrt(Pe/Rho)
!     
     do 12 ku = 1, nuez
     do 12 ju = 1, nuey
     do 12 iu = 1, nuex
       ijku  = iu+(ju-1)*nuex+(ku-1)*nuexy
       ii    = ijku+nuexyz*(ijk-1)
       wkfex = Gaussian(uex(iu),Vx,Vthe)
       fe(ii) = Rho*wkfex*(1.d0+alpha*dcos(ak*x(i)))
       if (iu == 1 .or. iu == nuex) fe(ii) = 0.d0    
12   continue
11 continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
!  calculate the number density for the electron
!
   allocate (wkRhoe(ncxyz), wkJex(ncxyz))
!
   call density (ncx,ncy,ncz,nuex,nuey,nuez,fe, &
                 Buex,Cuex,huex,Buey,Cuey,huey,Buez,Cuez,huez,wkRhoe)
   call currentX(ncx,ncy,ncz,nuex,nuey,nuez,fe,uex, &
                 Buex,Cuex,huex,Buey,Cuey,huey,Buez,Cuez,huez,wkJex)
!
!   do i = 1, ncx
!     print *, i, wkRhoe(i), wkJex(i)
!   enddo
!   print *, ' '
!   stop
!-------------------------------------------------------------------------------
!  use Poisson's equation to correct number density of electron 
!    distirbution function
!
!$OMP parallel do private(ii,ijk,ijku,Rho) collapse(3)
   do 16 k = 1, ncz
   do 16 j = 1, ncy
   do 16 i = 1, ncx
!
     ijk  = i+(j-1)*ncx+(k-1)*ncxy
     Rho  = 1.d0+alpha*dcos(ak*x(i))
!
     do 17 ku = 1, nuez
     do 17 ju = 1, nuey
     do 17 iu = 1, nuex
       ijku  = iu+(ju-1)*nuex+(ku-1)*nuexy
       ii    = ijku+nuexyz*(ijk-1)
       fe(ii) = fe(ii)/wkRhoe(ijk)*Rho     
17   continue
16 continue
!$OMP end parallel do
!
   print *, 'save fe'
   write (1) fe
!
   call density(ncx,ncy,ncz,nuex,nuey,nuez,fe, &
                Buex,Cuex,huex,Buey,Cuey,huey,Buez,Cuez,huez,wkRhoe)
   call currentX(ncx,ncy,ncz,nuex,nuey,nuez,fe,uex, &
                 Buex,Cuex,huex,Buey,Cuey,huey,Buez,Cuez,huez,wkJex)
!
!   do i = 1, ncx
!     print *, i, wkRhoe(i), wkJex(i)
!   enddo
!   print *, ' '
!
   deallocate (fe,wkJex)
!   stop
!===============================================================================
!  get the initial condition of fi
!
   allocate (fi(nofi))
!$OMP parallel do private(ii,ijk,ijku,Rho,Vx,ADi,aPi,Vthi,wkfix,wkfiy,wkfiz) collapse(3)
   do 21 k = 1, ncz
   do 21 j = 1, ncy
   do 21 i = 1, ncx
!
     ijk  = i+(j-1)*ncx+(k-1)*ncxy
     Rho  = Rho0
     Vx   = Vx0
     aPi  = Pi0
     Vthi = dsqrt(aPi/Rho/ami)
!     
     do 22 ku = 1, nuiz
     do 22 ju = 1, nuiy
     do 22 iu = 1, nuix
       ijku  = iu+(ju-1)*nuix+(ku-1)*nuixy
       ii    = ijku+nuixyz*(ijk-1)
       wkfix = Gaussian(uix(iu),Vx ,Vthi)
       fi(ii) = Rho*wkfix
       if (iu == 1 .or. iu == nuix) fi(ii) = 0.d0     
22   continue
21 continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
!  calculate the number density for the ion
!
   allocate (wkRhoi(ncxyz), wkJix(ncxyz))
!
   call density (ncx,ncy,ncz,nuix,nuiy,nuiz,fi, &
                 Buix,Cuix,huix,Buiy,Cuiy,huiy,Buiz,Cuiz,huiz,wkRhoi)
   call currentX(ncx,ncy,ncz,nuix,nuiy,nuiz,fi,uix, &
                 Buix,Cuix,huix,Buiy,Cuiy,huiy,Buiz,Cuiz,huiz,wkJix)
!
!   do i = 1, ncx
!     print *, i, wkRhoi(i), wkJix(i)
!   enddo
!   stop
!-------------------------------------------------------------------------------
!  use Poisson's equation to correct number density of ion 
!    distirbution function
!
!$OMP parallel do private(ii,ijk,ijku,Rho) collapse(3)
   do 26 k = 1, ncz
   do 26 j = 1, ncy
   do 26 i = 1, ncx
!
     ijk  = i+(j-1)*ncx+(k-1)*ncxy
     Rho  = Rho0
!
     do 27 ku = 1, nuiz
     do 27 ju = 1, nuiy
     do 27 iu = 1, nuix
       ijku  = iu+(ju-1)*nuix+(ku-1)*nuixy
       ii    = ijku+nuixyz*(ijk-1)
       fi(ii) = fi(ii)/wkRhoi(ijk)*Rho     
27   continue
26 continue
!$OMP end parallel do
!
   print *, 'save fi'
   write (1) fi

   call density(ncx,ncy,ncz,nuix,nuiy,nuiz,fi, &
                Buix,Cuix,huix,Buiy,Cuiy,huiy,Buiz,Cuiz,huiz,wkRhoi)
   call currentX(ncx,ncy,ncz,nuix,nuiy,nuiz,fi,uix, &
                 Buix,Cuix,huix,Buiy,Cuiy,huiy,Buiz,Cuiz,huiz,wkJix)
!
!   do i = 1, ncx
!     print *, i, wkRhoi(i), wkJix(i)
!   enddo
!
   deallocate (fi,wkRhoi,wkJix)
!   stop
!===============================================================================
!  get the initial condition of Bx
!
   allocate (Bx(ncxyz))
   do 30 k = 1, ncz
   do 30 j = 1, ncy
   do 30 i = 1, ncx
     ijk = i+(j-1)*ncx+(k-1)*ncxy
     Bx(ijk) = 0.0 !0.5d0*(Bxdif*tanhX(i)+Bxsum)
30 continue
   print *, 'save Bx'
   write (1) Bx
   deallocate (Bx)
!-------------------------------------------------------------------------------
!  get the initial condition of By
!
   allocate (By(ncxyz))
   do 31 k = 1, ncz
   do 31 j = 1, ncy
   do 31 i = 1, ncx
     ijk = i+(j-1)*ncx+(k-1)*ncxy
     By(ijk) = 0.0 !0.5d0*(Bydif*tanhX(i)+Bysum)
31 continue
   print *, 'save By'
   write (1) By
   deallocate (By)
!-------------------------------------------------------------------------------
!  get the initial condition of Bz
!
   allocate (Bz(ncxyz))
   do 32 k = 1, ncz
   do 32 j = 1, ncy
   do 32 i = 1, ncx
     ijk = i+(j-1)*ncx+(k-1)*ncxy
     Bz(ijk) = 0.0 !0.5d0*(Bzdif*tanhX(i)+Bzsum)
32 continue
   print *, 'save Bz'
   write (1) Bz
   deallocate (Bz)
!===============================================================================
!  get the initial condition of Ex
!
   allocate (Ex(ncx))
!
   do i = 1, ncx
     Ex(i) = -alpha*dsin(ak*x(i))/ak
   enddo
!
   print *, 'save Ex'
   write (1) Ex
   print *, 'CFL = ', (uex(nuex)/dx+maxval(Ex)/duex)*dt
   deallocate (Ex)
!-------------------------------------------------------------------------------
!  get the initial condition of Ey
!
   allocate (Ey(ncxyz))
   Ey = 0.d0
   print *, 'save Ey'
   write (1) Ey
   deallocate (Ey)
!-------------------------------------------------------------------------------
!  get the initial condition of Ez
!
   allocate (Ez(ncxyz))
   Ez = 0.d0
   print *, 'save Ez'
   write (1) Ez
   deallocate (Ez)
   close (1)
!===============================================================================
!  get the initial condition of etaB
!
   open (2,file='eta.bin',status='unknown',form='unformatted')
   allocate (etaB(ncxyz))
   etaB = 0.00d0
   print *, 'save etaB'
   write (2) etaB
   deallocate (etaB)
!-------------------------------------------------------------------------------
!  get the initial condition of etafe
!
   allocate (etafe(ncxyz))
   etafe = 0.00d0
   print *, 'save etafe'
   write (2) etafe
   deallocate (etafe)
!-------------------------------------------------------------------------------
!  get the initial condition of etafi
!
   allocate (etafi(ncxyz))
   etafi = 0.000d0 
   print *, 'save etafi'
   write (2) etafi
   deallocate (etafi)
   close (2)
!===============================================================================
   deallocate (x, y, z)
   deallocate (uex, uey, uez)
   deallocate (uix, uiy, uiz)
end program initial_Vlasov_ES1D1V
!===============================================================================
