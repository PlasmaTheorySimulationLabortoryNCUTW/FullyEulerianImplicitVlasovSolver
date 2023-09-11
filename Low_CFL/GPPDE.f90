!===============================================================================
module GPPDE
!
!  the kernel code of the General-Purpose PDE solver in Cartesian coordinate
!
     use MPI3D
     use coef
     use omp_lib
   contains
!===============================================================================
!--SUB. TVDRKDE4th
!  Fourth-order total variation diminishing (TVD) Runge-Kutta method
!    for solving system ODEs
!-------------------------------------------------------------------------------
   subroutine TVDRKDE4th(funcfp,f,funcin,funcout,t,dt,Ndim)
!
!  Fourth order Runge-Kutta method for solving 1st-order system 
!    differential equations of the form df/dt = funcfp
!
!  Input Data : dt
!  Input & Output Data  : t
!  Input & Output Array : f(Ndim)
!  Input Array    : funcin(Ndim)
!  Output Array   : funcout(Ndim)
!
!  Working Arrays : wkf1(Ndim), wkf2(Ndim), wkf3(Ndim)
!  Working Arrays : wkfunc1(Ndim), wkfunc2(Ndim)
!
!  External Function : funcfp
!
     implicit double precision (a-h,o-z)
     dimension :: f(1), funcin(1), funcout(1)
     real*8, allocatable, dimension(:) :: wkf1, wkf2, wkf3
     external :: funcfp
     allocate (wkf1(Ndim), wkf2(Ndim), wkf3(Ndim))
!
     aq = 0.25d0
     ah = 0.5d0
     a13 = 1.d0/3.d0
     a16 = 1.d0/6.d0
     a19 = 1.d0/9.d0
     a29 = 2.d0/9.d0
     a23 = 2.d0/3.d0
!
     do i = 1, Ndim
       wkf1(i) = DIFAB(f(i),-ah*dt*funcin(i))
!      wkf1(i) = f(i)+ah*dt*funcin(i)
!
       wkf2(i) = DIFAB(ah*f(i),aq*dt*funcin(i))
!      wkf2(i) = ah*f(i)-aq*dt*funcin(i)
!
       wkf3(i) = a19*DIFAB(f(i),dt*funcin(i))
!      wkf3(i) = a19*(f(i)-dt*funcin(i))
     enddo
!
     twork = t+ah
     call funcfp(twork,wkf1,funcout,Ndim)
     do i = 1, Ndim
       wkf2(i) = sumABC(wkf2(i),ah*wkf1(i),ah*dt*funcout(i))
!      wkf2(i) = wkf2(i)+ah*wkf1(i)+ah*dt*funcout(i)
!
       wkf3(i) = sumABC(wkf3(i),a29*wkf1(i),-a13*dt*funcout(i))
!      wkf3(i) = wkf3(i)+a29*wkf1(i)-a13*dt*funcout(i)
!
       f(i)    = DIFAB(a13*wkf1(i),-a16*dt*funcout(i))
!      f(i)    = a13*wkf1(i)+a16*dt*funcout(i)
     enddo
!
     call funcfp(twork,wkf2,funcout,Ndim)
     do i = 1, Ndim
       wkf3(i) = sumABC(wkf3(i),a23*wkf2(i),dt*funcout(i))
!      wkf3(i) = wkf3(i)+a23*wkf2(i)+dt*funcout(i)
!
       f(i)    = DIFAB(f(i),-a13*wkf2(i))
!      f(i)    = f(i)+a13*wkf2(i)
     enddo
!
     call funcfp(twork,wkf3,funcout,Ndim)
     do i = 1, Ndim
       f(i) = sumABC(f(i),a13*wkf3(i),a16*dt*funcout(i))     
!      f(i) = f(i)+a13*wkf3(i)+a16*dt*funcout(i)
     enddo     
!
     call funcfp(t,f,funcout,Ndim)
     deallocate (wkf1,wkf2,wkf3)
!
   end subroutine TVDRKDE4th
!===============================================================================
!===============================================================================
!--SUB. RKDE4th(Fourth order Runge-Kutta method for solving system ODEs)
!-------------------------------------------------------------------------------
   subroutine RKDE4th(funcfp,f,funcin,funcout,wkf1,wkf2,t,dt,Ndim)
!
!  Fourth order Runge-Kutta method for solving 1st-order system 
!    differential equations of the form df/dt = funcfp
!
!  Input Data : dt
!  Input & Output Data  : t
!  Input & Output Array : f(Ndim)
!  Input Array    : funcin(Ndim)
!  Output Array   : funcout(Ndim)
!  Working Arrays : wkf1(Ndim), wkf2(Ndim)
!
!  External Function : funcfp
!
     implicit double precision (a-h,o-z)
     dimension :: f(1), funcin(1), funcout(1), wkf1(1), wkf2(1)
     external :: funcfp
!
     ah = 0.5d0*dt
     a1v3 = dt/3.d0
     a1v6 = dt/6.d0
     wkf1(1:Ndim) = f(1:Ndim)
!     do i = 1, Ndim
!       wkf1(i) = f(i)
!     enddo
!
!$OMP parallel do private(temp1,temp2)
     do i = 1, Ndim
       temp1 = -ah*funcin(i)
       wkf2(i) = DIFAB(f(i),temp1)
!      wkf2(i) = f(i)+ah*funcin(i)
       temp2 = -a1v6*funcin(i)
       wkf1(i) = DIFAB(wkf1(i),temp2)
!      wkf1(i) = wkf1(i)+a1v6*funcin(i)
     enddo
!$OMP end parallel do
!
     twork = t+ah
     call funcfp(twork,wkf2,funcout,Ndim)
!
!$OMP parallel do private(temp1,temp2)
     do i = 1, Ndim
       temp1 = -ah*funcout(i)
       wkf2(i) = DIFAB(f(i),temp1)
!      wkf2(i) = f(i)+ah*funcout(i)
       temp2 = -a1v3*funcout(i)
       wkf1(i) = DIFAB(wkf1(i),temp2)
!      wkf1(i) = wkf1(i)+a1v3*funcout(i)
     enddo
!$OMP end parallel do
!
     call funcfp(twork,wkf2,funcout,Ndim)
!
!$OMP parallel do private(temp1,temp2)
     do i = 1, Ndim
       temp1 = -dt*funcout(i)
       wkf2(i) = DIFAB(f(i),temp1)
!      wkf2(i) = f(i)+dt*funcout(i)
       temp2 = -a1v3*funcout(i)
       wkf1(i) = DIFAB(wkf1(i),temp2)
!      wkf1(i) = wkf1(i)+a1v3*funcout(i)
     enddo
!$OMP end parallel do
!
     t = t+dt
     call funcfp(t,wkf2,funcout,Ndim)
!
!$OMP parallel do private(temp1)
     do i = 1, Ndim
       temp1 = -a1v6*funcout(i)
       f(i) = DIFAB(wkf1(i),temp1) 
!      f(i) = wkf1(i)+a1v6*funcout(i)
     enddo
!$OMP end parallel do
!
     call funcfp(t,f,funcout,Ndim)
!
   end subroutine RKDE4th
!===============================================================================
!===============================================================================
!--SUB. RKDE2nd(Second order Runge-Kutta method for solving system ODEs)
!-------------------------------------------------------------------------------
   subroutine RKDE2nd(funcfp,f,funcin,funcout,wkf1,wkf2,t,dt,Ndim)
!
!  Second order Runge-Kutta method for solving 1st-order system 
!    differential equations of the form df/dt = funcfp
!
!  Input Data : dt
!  Input & Output Data  : t
!  Input & Output Array : f(Ndim)
!  Input Array    : funcin(Ndim)
!  Output Array   : funcout(Ndim)
!  Working Arrays : wkf1(Ndim), wkf2(Ndim)
!
!  External Function : funcfp
!
     implicit double precision (a-h,o-z)
     dimension :: f(1), funcin(1), funcout(1), wkf1(1), wkf2(1)
     external :: funcfp
!
     ah = 0.5d0*dt
!
     do i = 1, Ndim
       temp = -ah*funcin(i)
       wkf2(i) = DIFAB(f(i),temp)
!      wkf2(i) = f(i)+ah*funcin(i)
     enddo
!
     twork = t+ah
     call funcfp(twork,wkf2,funcout,Ndim)
     do i = 1, Ndim
       temp = -dt*funcout(i)
       f(i) = DIFAB(f(i),temp)
!      f(i) = f(i)+dt*funcout(i)
     enddo
!
     t = t+dt
     call funcfp(t,f,funcin,Ndim)
!
   end subroutine RKDE2nd
!===============================================================================
!===============================================================================
!  sub. PCDE4th0
!-------------------------------------------------------------------------------
   subroutine PCDE4th0(funcfp,f,func,funcm1,funcm2,funcm3,wkf1,wkf2,t,dt,Ndim)
!
!  Using 4th order Predict-Corrector Method to solve first order
!    differential equation of the form df/dt = funcfp
!    Preparing func, funcm1, funcm2, funcm3 from 4th-order Runge-Kutta method.
!    At the same time, advance t for three time steps (t=0 --> 3*dt).
!
!  Input Data : dt
!  Input & Output Data  : t
!  Input & Output Array : f(Ndim)
!  Output Arrays  : func(Ndim), funcm1(Ndim), funcm2(Ndim), funcm3(Ndim)
!  Working Arrays : wkf1(Ndim), wkf2(Ndim)
!
!  External Function : funcfp
!
     implicit double precision (a-h,o-z)
     dimension :: f(1), func(1), funcm1(1), funcm2(1), funcm3(1)
     dimension :: wkf1(1), wkf2(1)
     external :: funcfp
!
!  f and funcm3 at t = 0
!
     call funcfp(t,f,funcm3,Ndim)
!
!  f and funcm2 at t = dt
!
     call RKDE4th(funcfp,f,funcm3,funcm2,wkf1,wkf2,t,dt,Ndim)
!
!  f and funcm1 at t = 2*dt
!
     call RKDE4th(funcfp,f,funcm2,funcm1,wkf1,wkf2,t,dt,Ndim)
!
!  f and func at t = 3*dt
!
     call RKDE4th(funcfp,f,funcm1,func,wkf1,wkf2,t,dt,Ndim)
!
   end subroutine PCDE4th0
!===============================================================================
!===============================================================================
!--sub. PCDE4th1
!-------------------------------------------------------------------------------
   subroutine PCDE4th1(funcfp,f,func,funcm1,funcm2,funcm3,wkfsum,wkf,t,dt,Ndim)
!
!  Using 4th order Predict-Corrector Method to solve first order
!    differential equation of the form df/dt = funcfp
!
!  Input Data : dt
!  Input & Output Data   : t
!  Input & Output Arrays : f(Ndim), func(Ndim), funcm1(Ndim),
!                          funcm2(Ndim), funcm3(Ndim)
!  Working Arrays : wkfsum(Ndim), wkf(Ndim)
!
!  External Function : funcfp
!
     implicit double precision (a-h,o-z)
     dimension :: f(1), func(1), funcm1(1), funcm2(1), funcm3(1)
     dimension :: wkfsum(1), wkf(1)
     logical :: index
     external :: funcfp
!
     a1v24 = dt/24.d0
     a3v8 = 0.375d0*dt
     a55 = 55.d0
     a59 = 59.d0
     a37 = 37.d0
     a19 = 19.d0
     a09 = 9.d0
     a05 = 5.d0
     t = t+dt
!
!$OMP parallel do private(t55,t59,t37,t09,t19,t05,sum,sum1,sum2)
     do i = 1, Ndim
!-------------------------------------------------------------------------------
!  4th order Adams' Open Formula
!
       t55 = a55*func(i)
       t59 = a59*funcm1(i)
       t37 = a37*funcm2(i)
       t09 = a09*funcm3(i)
       sum1 = DIFAB(t55,t59)
       sum2 = DIFAB(t09,t37)
       sum  = DIFAB(sum1,sum2)
       sum  = -a1v24*sum
       wkf(i) = DIFAB(f(i),sum)
!      wkf(i) = f(i)+a1v24*(a55*func(i)-a59*funcm1(i) &
!              +a37*funcm2(i)-a09*funcm3(i))
!-------------------------------------------------------------------------------
!  4th order Adams' Closed Formula (prepare)
!
       t19 = a19*func(i)
       t05 = a05*funcm1(i)
       sum1 = DIFAB(t19,t05)
       sum  = DIFAB(sum1,-funcm2(i))
       sum  = -a1v24*sum
       wkfsum(i) = DIFAB(f(i),sum)
!      wkfsum(i) = f(i)+a1v24*(a19*func(i)-a05*funcm1(i)+funcm2(i))
     enddo
!$OMP end parallel do
!-------------------------------------------------------------------------------
     ii = 0
  16 continue
     ii = ii+1
     if (MyID == 0) write (*,100) 'iteration for checking =', ii
 100 format (a24,i3)
!-------------------------------------------------------------------------------
!  Main loop of iteration
!
     kk = iter  ! number of iteration per cycle
     kkm1 = kk-1
     do k = 1, kkm1
       call funcfp(t,wkf,funcm3,Ndim)
!$OMP parallel do private(t38)
       do i = 1, Ndim
         t38 = a3v8*funcm3(i)
         wkf(i) = DIFAB(t38,-wkfsum(i))
!        wkf(i) = a3v8*funcm3(i)+wkfsum(i)
       enddo
!$OMP end parallel do
     enddo
!-------------------------------------------------------------------------------
!  Preparing for converged checking
!
     f(1:Ndim) = wkf(1:Ndim)
!     do i = 1, Ndim
!       f(i) = wkf(i)
!     enddo
!
     call funcfp(t,wkf,funcm3,Ndim)
!
!$OMP parallel do private(t38)
     do i = 1, Ndim
       t38 = a3v8*funcm3(i)
       wkf(i) = DIFAB(t38,-wkfsum(i))
!       wkf(i) = a3v8*funcm3(i)+wkfsum(i)
     enddo
!$OMP end parallel do
!
     if (ccheck == 'N') goto 17 ! Don't check the convergence condition.
!
!  Converged checking
!
     call check_converge(wkf,f,index,Ndim)
!
     if (index) then
       goto 17
     else
!
       if (ii >= 10)	then
         !print *, 'Program stop at SUB. PCDE4th1 !!!'
         !if (MyID == 0) print *, 'Program do not converge !!!'
         call MPI_FINALIZE(IERR)
         stop
       endif
!
       goto 16
     endif
!
  17 continue
!
!  update f, func, funcm1, funcm2 and funcm3
!
     funcm3(1:Ndim) = funcm2(1:Ndim)
     funcm2(1:Ndim) = funcm1(1:Ndim)
     funcm1(1:Ndim) = func(1:Ndim)
     f(1:Ndim)      = wkf(1:Ndim)
!
!     do i = 1, Ndim
!       funcm3(i) = funcm2(i)
!     enddo
!
!     do i = 1, Ndim
!       funcm2(i) = funcm1(i)
!     enddo
!
!     do i = 1, Ndim
!       funcm1(i) = func(i)
!     enddo
!
!     do i = 1, Ndim
!       f(i) = wkf(i)
!     enddo
!
     call funcfp(t,f,func,Ndim)
!
   end subroutine PCDE4th1
!===============================================================================
!===============================================================================
!--SUB. Check_converge
!-------------------------------------------------------------------------------
   subroutine Check_converge(f1,f2,index,Ndim)
!
!  Convergence checking for the 4th order Adams' Closed Formula
!
!  Input Arrays   : f1(Ndim), f2(Ndim)
!  Output Logical : index
!
     implicit double precision (a-h,o-z)
     dimension :: f1(1), f2(1)
     dimension :: fmax(1), fmax_global(1)
     dimension :: ii(1), ii_global(1)
     logical :: index
!
     index = .false.
!
     fmax1 = maxval(f1)
     fmax2 = maxval(f2)
     fmax(1) = max(fmax1,fmax2,1.d0)
!
     call MPI_ALLREDUCE(fmax(1),fmax_global(1),1, &
                        MPI_REAL8,MPI_MAX,MPI_WORLD,IERR)

     abserr = fmax_global(1)*relerr
!
     isum = 0
!
!$OMP parallel do private(err) reduction(+:ii)
     do k = 1, Ndim
       err = dabs(f1(k)-f2(k))
       if (err > abserr) isum = isum+1
     enddo
!$OMP end parallel do
!
     ii(1) = isum
!
     call MPI_ALLREDUCE(ii(1),ii_global(1),1, &
                        MPI_INTEGER,MPI_MAX,MPI_WORLD,IERR)
!
     if (ii_global(1) == 0) index = .true.
!
   end subroutine Check_converge
!===============================================================================
!===============================================================================
!  sub. PCDE2nd0
!-------------------------------------------------------------------------------
   subroutine PCDE2nd0(funcfp,f,func,funcm1,wkf1,wkf2,t,dt,Ndim)
!
!  Using 2nd-order Predict-Corrector Method to solve first order
!    differential equation of the form df/dt = funcfp
!    Preparing func, funcm1 from 4th-order Runge-Kutta method.
!    At the same time, advance t for three time steps (t=0 --> dt).
!
!  Input Data : dt
!  Input & Output Data  : t
!  Input & Output Array : f(Ndim)
!  Output Arrays  : func(Ndim), funcm1(Ndim)
!  Working Arrays : wkf1(Ndim), wkf2(Ndim)
!
!  External Function : funcfp
!
     implicit double precision (a-h,o-z)
     dimension :: f(1), func(1), funcm1(1)
     dimension :: wkf1(1), wkf2(1)
     external :: funcfp
!
!  f and funcm3 at t = 0
!
     call funcfp(t,f,funcm1,Ndim)
!
!  f and func at t = dt
!
     call RKDE4th(funcfp,f,funcm1,func,wkf1,wkf2,t,dt,Ndim)
!
   end subroutine PCDE2nd0
!===============================================================================
!===============================================================================
!--sub. PCDE2nd1
!-------------------------------------------------------------------------------
   subroutine PCDE2nd1(funcfp,f,func,funcm1,wkfsum,wkf,t,dt,Ndim)
!
!  Using 2nd order Predict-Corrector Method to solve first order
!    differential equation of the form df/dt = funcfp
!
!  Input Data : dt
!  Input & Output Data   : t
!  Input & Output Arrays : f(Ndim), func(Ndim), funcm1(Ndim)
!  Working Arrays : wkfsum(Ndim), wkf(Ndim)
!
!  External Function : funcfp
!
     implicit double precision (a-h,o-z)
     dimension :: f(1), func(1), funcm1(1)
     dimension :: wkfsum(1), wkf(1)
     logical :: index
     external :: funcfp
!
     ah   = 0.5d0*dt
     a3 = 3.d0
     t = t+dt
!
     do i = 1, Ndim
!-------------------------------------------------------------------------------
!  2nd order Adams' Open Formula
!
       temp = a3*func(i)
       sum = -ah*DIFAB(temp,funcm1(i)) 
       wkf(i) = DIFAB(f(i),sum)
!      wkf(i) = f(i)+ah*(a3*func(i)-funcm1(i))
!-------------------------------------------------------------------------------
!  2nd order Adams' Closed Formula (prepare)
!
       temp = -ah*func(i)
       wkfsum(i) = DIFAB(f(i),temp)
!      wkfsum(i) = f(i)+ah*func(i)
     enddo
!-------------------------------------------------------------------------------
     ii = 0
  16 continue
     ii = ii+1
!     if (MyID == 0) write (*,100) 'iteration for checking =', ii
 100 format (a24,i3)
!-------------------------------------------------------------------------------
!  Main loop of iteration
!
     kk = iter  ! number of iteration per cycle
     kkm1 = kk-1
     do k = 1, kkm1
       call funcfp(t,wkf,funcm1,Ndim)
       do i = 1, Ndim
         temp = ah*funcm1(i)
         wkf(i) = DIFAB(temp,-wkfsum(i))
!        wkf(i) = ah*funcm1(i)+wkfsum(i)
       enddo
     enddo
!-------------------------------------------------------------------------------
!  Preparing for converged checking
!
     do i = 1, Ndim
       f(i) = wkf(i)
     enddo
!
     call funcfp(t,wkf,funcm1,Ndim)
!
     do i = 1, Ndim
       temp = ah*funcm1(i)
       wkf(i) = DIFAB(temp,-wkfsum(i))
!      wkf(i) = ah*funcm1(i)+wkfsum(i)
     enddo
!
     if (ccheck == 'N') goto 17 ! Don't check the convergence condition.
!
!  Converged checking
!
     call check_converge(wkf,f,index,Ndim)
!
     if (index) then
       goto 17
     else
!
       if (ii >= 10)	then
         !print *, 'Program stop at SUB. PCDE2nd1 !!!'
         !if (MyID == 0) print *, 'Program do not converge !!!'
         call MPI_FINALIZE(IERR)
         stop
       endif
!
       goto 16
     endif
!
  17 continue
!
!  update f, func, funcm1
!
     do i = 1, Ndim
       funcm1(i) = func(i)
     enddo
!
     do i = 1, Ndim
       f(i) = wkf(i)
     enddo
!
     call funcfp(t,f,func,Ndim)
!
   end subroutine PCDE2nd1
!===============================================================================
end module GPPDE
!===============================================================================