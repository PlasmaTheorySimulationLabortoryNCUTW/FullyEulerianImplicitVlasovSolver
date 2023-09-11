!===============================================================================
module vector3D
!
!  the kernel code of the General-Purpose PDE solver in Cartesian coordinate
!
     use MPI3D
     use CenDif3D_MPI
     use coef
     use WENO3D_MPI
   contains
!===============================================================================
!  sub. dot3d
   subroutine dot3d(Ax,Ay,Az,Bx,By,Bz,C,Nx,Ny,Nz)
!
!  C = A dot B
!
!  Input Arrays : Ax(Nx,Ny,Nz), Ay(Nx,Ny,Nz), Az(Nx,Ny,Nz),
!                 Bx(Nx,Ny,Nz), By(Nx,Ny,Nz), Bz(Nx,Ny,Nz)
!  Output Array : C(Nx,Ny,Nz)
!
!  for 2D Arrays : Nz = 1
!  for 1D Arrays : Ny = Nz = 1
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: Ax(1), Ay(1), Az(1)
     dimension :: Bx(1), By(1), Bz(1)
     dimension :: C(1)
!-------------------------------------------------------------------------------
     Nxy = Nx*Ny
     do 10 k = 1, Nz
     do 10 j = 1, Ny
     do 10 i = 1, Nx
       ijk = i+(j-1)*Nx+(k-1)*Nxy
       sumx = Ax(ijk)*Bx(ijk)
       sumy = Ay(ijk)*By(ijk)
       sumz = Az(ijk)*Bz(ijk)
       sumxy = DIFAB(sumx,-sumy)
       C(ijk) = DIFAB(sumxy,-sumz)
!      C(ijk) = sumx+sumy+sumz
  10 continue
!-------------------------------------------------------------------------------
   end subroutine dot3d
!===============================================================================
!===============================================================================
!  sub. corss3d
   subroutine cross3d(Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Nx,Ny,Nz)
!
!  C = A cross B
!
!  Input Arrays  : Ax(Nx,Ny,Nz), Ay(Nx,Ny,Nz), Az(Nx,Ny,Nz),
!                  Bx(Nx,Ny,Nz), By(Nx,Ny,Nz), Bz(Nx,Ny,Nz)
!  Output Arrays : Cx(Nx,Ny,Nz), Cy(Nx,Ny,Nz), Cz(Nx,Ny,Nz)
!
!  for 2D Arrays : Nz = 1
!  for 1D Arrays : Ny = Nz = 1
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: Ax(1), Ay(1), Az(1)
     dimension :: Bx(1), By(1), Bz(1)
     dimension :: Cx(1), Cy(1), Cz(1)
!
     Nxy = Nx*Ny
     do 10 k = 1, Nz
     do 10 j = 1, Ny
     do 10 i = 1, Nx
       ijk = i+(j-1)*Nx+(k-1)*Nxy
       sum1 = Ay(ijk)*Bz(ijk)
       sum2 = Az(ijk)*By(ijk) 
       Cx(ijk) = DIFAB(sum1,sum2)
!
       sum1 = Az(ijk)*Bx(ijk)
       sum2 = Ax(ijk)*Bz(ijk)
       Cy(ijk) = DIFAB(sum1,sum2)
!
       sum1 = Ax(ijk)*By(ijk)
       sum2 = Ay(ijk)*Bx(ijk)
       Cz(ijk) = DIFAB(sum1,sum2)
!      Cx(ijk) = Ay(ijk)*Bz(ijk)-Az(ijk)*By(ijk) 
!      Cy(ijk) = Az(ijk)*Bx(ijk)-Ax(ijk)*Bz(ijk)
!      Cz(ijk) = Ax(ijk)*By(ijk)-Ay(ijk)*Bx(ijk)
  10 continue
!-------------------------------------------------------------------------------      
   end subroutine cross3d
!===============================================================================
!===============================================================================
!  sub. grad3d_1st
   subroutine grad3d_1st(Nx,Ny,Nz,A,Bx,By,Bz,hx,WAx,WBx,WCx, & 
                         hy,WAy,WBy,WCy,hz,WAz,WBz,WCz,ixp,iyp,izp)
!-------------------------------------------------------------------------------
!  B = gradient A   use 1st-order central difference
!
!  Input Arrays   : A(Nx,Ny,Nz),
!                   hx(Nx), WAx(Nx), WBx(Nx), WCx(Nx) 
!                   hy(Ny), WAy(Ny), WBy(Ny), WCy(Ny) 
!                   hz(Nz), WAz(Nz), WBz(Nz), WCz(Nz)
!
!  Output Arrays  : Bx(Nx,Ny,Nz), By(Nx,Ny,Nz), Bz(Nx,Ny,Nz)
!
!  ixp, iyp, izp \= 1 : Free Boundary Condition
!  ixp, iyp, izp =  1 : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension A(1), Bx(1), By(1), Bz(1)
     dimension hx(1), WAx(1), WBx(1), WCx(1)
     dimension hy(1), WAy(1), WBy(1), WCy(1)
     dimension hz(1), WAz(1), WBz(1), WCz(1)
!-------------------------------------------------------------------------------
! 
!  Bx = dA/dx
!
     if (ixp .eq. 1) then
!      Periodic Boundary Condition along x-direction
       call PCDx1st3D_MPI(A,WAx,WBx,WCx,Bx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)
     else
!      Free Boundary Condition along x-direction
       call FCDx1st3D_MPI(A,WAx,WBx,WCx,Bx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)
     endif
!-------------------------------------------------------------------------------
! 
!  By = dA/dy
!
     if (Ny == 1 .and. Nz == 1) then
     	 By(1) = 0.d0
     	 Bz(1) = 0.d0
     	 return
     endif
!
     if (iyp .eq. 1) then 
!      Periodic Boundary Condition along y-direction
       call PCDy1st3D_MPI(A,WAy,WBy,WCy,By,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)
     else 
!      Non-Periodic Boundary Condition along y-direction
       call FCDy1st3D_MPI(A,WAy,WBy,WCy,By,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)
     endif
!-------------------------------------------------------------------------------
! 
!  Bz = dA/dz
!
     if (Nz == 1) then
     	 Bz(1) = 0.d0
     	 return
     endif
!
     if (izp .eq. 1) then 
!      Periodic Boundary Condition along z-direction
       call PCDz1st3D_MPI(A,WAz,WBz,WCz,Bz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)
     else
!      Non-Periodic Boundary Condition along z-direction
       call FCDz1st3D_MPI(A,WAz,WBz,WCz,Bz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)
     endif
!
   end subroutine grad3d_1st
!===============================================================================
!===============================================================================
!  sub. grad3d_3rd
   subroutine grad3d_3rd(Nx,Ny,Nz,A,Bx,By,Bz,    &
                         hx,WAx,WBx,WCx,WDx,WEx, & 
                         hy,WAy,WBy,WCy,WDy,WEy, &
                         hz,WAz,WBz,WCz,WDz,WEz, &
                         ixp,iyp,izp)
!-------------------------------------------------------------------------------
!  B = gradient A   use 3rd-order central difference
!
!  Input Arrays   : A(Nx,Ny,Nz),
!                   hx(Nx), WAx(Nx), WBx(Nx), WCx(Nx), WDx(Nx), WEx(Nx) 
!                   hy(Ny), WAy(Ny), WBy(Ny), WCy(Ny), WDy(Ny), WEy(Ny) 
!                   hz(Nz), WAz(Nz), WBz(Nz), WCz(Nz), WDz(Nz), WEz(Nz)
!
!  Output Arrays  : Bx(Nx,Ny,Nz), By(Nx,Ny,Nz), Bz(Nx,Ny,Nz)
!
!  ixp, iyp, izp \= 1 : Free Boundary Condition
!  ixp, iyp, izp =  1 : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension A(1), Bx(1), By(1), Bz(1)
     dimension hx(1), WAx(1), WBx(1), WCx(1), WDx(1), WEx(1)
     dimension hy(1), WAy(1), WBy(1), WCy(1), WDy(1), WEy(1)
     dimension hz(1), WAz(1), WBz(1), WCz(1), WDz(1), WEz(1)
!-------------------------------------------------------------------------------
! 
!  Bx = dA/dx
!
     if (ixp .eq. 1) then
!      Periodic Boundary Condition along x-direction
       call PCDx3rd3D_MPI(A,WAx,WBx,WCx,WDx,WEx,Bx,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)
     else
!      Free Boundary Condition along x-direction
       call FCDx3rd3D_MPI(A,WAx,WBx,WCx,WDx,WEx,Bx,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)
     endif
!-------------------------------------------------------------------------------
! 
!  By = dA/dy
!
     if (Ny == 1 .and. Nz == 1) then
     	 By(1) = 0.d0
     	 Bz(1) = 0.d0
     	 return
     endif
!
     if (iyp .eq. 1) then 
!      Periodic Boundary Condition along y-direction
       call PCDy3rd3D_MPI(A,WAy,WBy,WCy,WDy,WEy,By,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)
     else 
!      Non-Periodic Boundary Condition along y-direction
       call FCDy3rd3D_MPI(A,WAy,WBy,WCy,WDy,WEy,By,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)
     endif
!-------------------------------------------------------------------------------
! 
!  Bz = dA/dz
!
!
     if (Nz == 1) then
     	 Bz(1) = 0.d0
     	 return
     endif
!
     if (izp .eq. 1) then 
!      Periodic Boundary Condition along z-direction
       call PCDz3rd3D_MPI(A,WAz,WBz,WCz,WDz,WEz,Bz,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
     else
!      Non-Periodic Boundary Condition along z-direction
       call FCDz3rd3D_MPI(A,WAz,WBz,WCz,WDz,WEz,Bz,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
     endif
!
   end subroutine grad3d_3rd
!===============================================================================
!===============================================================================
!  sub. div3d_1st
   subroutine div3d_1st(Nx,Ny,Nz,Ax,Ay,Az,B, &
                        hx,WAx,WBx,WCx,hy,WAy,WBy,WCy,hz,WAz,WBz,WCz, &
                        ixp,iyp,izp)
!-------------------------------------------------------------------------------
!  B = divergence A  use 1st-order central difference
!
!  Input Arrays   : Ax(Nx,Ny,Nz), Ay(Nx,Ny,Nz), Az(Nx,Ny,Nz)
!                   hx(Nx), WAx(Nx), WBx(Nx), WCx(Nx), 
!                   hy(Ny), WAy(Ny), WBy(Ny), WCy(Ny), 
!                   hz(Nz), WAz(Nz), WBz(Nz), WCz(Nz)
!
!  Output Arrays  : B(Nx,Ny,Nz)
!
!  ixp, iyp, izp \= 1 : Free Boundary Condition
!  ixp, iyp, izp =  1 : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: Ax(1), Ay(1), Az(1), B(1)
     dimension :: hx(1), WAx(1), WBx(1), WCx(1)
     dimension :: hy(1), WAy(1), WBy(1), WCy(1)
     dimension :: hz(1), WAz(1), WBz(1), WCz(1)
     real*8, allocatable, dimension(:) :: Bsum
     Nxyz = Nx*Ny*Nz
     allocate (Bsum(Nxyz))
!-------------------------------------------------------------------------------
!  B = Bx+By+Bz = dAx/dx+dAy/dy+dAz/dz
!  Bx = dAx/dx
!
     if (ixp .eq. 1) then
!      Periodic Boundary Condition along x-direction
       call PCDx1st3D_MPI(Ax,WAx,WBx,WCx,B,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)
     else
!      Free Boundary Condition along x-direction
       call FCDx1st3D_MPI(Ax,WAx,WBx,WCx,B,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)
     endif
!-------------------------------------------------------------------------------
!  B = Bx+By+Bz = dAx/dx+dAy/dy+dAz/dz
!  Bx+By = dAx/dx+dAy/dy
!
     if (Ny == 1 .and. Nz == 1) return
!
     if (iyp .eq. 1) then 
!      Periodic Boundary Condition along y-direction
       call PCDy1st3D_MPI(Ay,WAy,WBy,WCy,Bsum,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)
     else 
!      Non-Periodic Boundary Condition along y-direction
       call FCDy1st3D_MPI(Ay,WAy,WBy,WCy,Bsum,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)
     endif
!
     do ijk = 1, Nxyz
       B(ijk) = DIFAB(B(ijk),-Bsum(ijk))
!      B(ijk) = B(ijk)+Bsum(ijk)
     enddo
!-------------------------------------------------------------------------------
!  B = Bx+By+Bz = dAx/dx+dAy/dy+dAz/dz
!
     if (Nz == 1) return
!
     if (izp .eq. 1) then 
!      Periodic Boundary Condition along z-direction
       call PCDz1st3D_MPI(Az,WAz,WBz,WCz,Bsum,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
     else
!      Non-Periodic Boundary Condition along z-direction
       call FCDz1st3D_MPI(Az,WAz,WBz,WCz,Bsum,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
     endif
!
     do ijk = 1, Nxyz
       B(ijk) = DIFAB(B(ijk),-Bsum(ijk))
!      B(ijk) = B(ijk)+Bsum(ijk)
     enddo
!
     deallocate (Bsum)
   end subroutine div3d_1st
!===============================================================================
!===============================================================================
!  sub. div3d_3rd
   subroutine div3d_3rd(Nx,Ny,Nz,Ax,Ay,Az,B,    &
                        hx,WAx,WBx,WCx,WDx,WEx, &
                        hy,WAy,WBy,WCy,WDy,WEy, &    
                        hz,WAz,WBz,WCz,WDz,WEz, &
                        ixp,iyp,izp)
!-------------------------------------------------------------------------------
!  B = divergence A  use 3rd-order central difference
!
!  Input Arrays   : Ax(Nx,Ny,Nz), Ay(Nx,Ny,Nz), Az(Nx,Ny,Nz)
!                   hx(Nx), WAx(Nx), WBx(Nx), WCx(Nx), WDx(Nx), WEx(Nx), 
!                   hy(Ny), WAy(Ny), WBy(Ny), WCy(Ny), WDy(Ny), WEy(Ny), 
!                   hz(Nz), WAz(Nz), WBz(Nz), WCz(Nz), WDz(Nz), WEz(Nz)
!
!  Output Arrays  : B(Nx,Ny,Nz)
!
!  ixp, iyp, izp \= 1 : Free Boundary Condition
!  ixp, iyp, izp =  1 : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: Ax(1), Ay(1), Az(1), B(1)
     dimension :: hx(1), WAx(1), WBx(1), WCx(1), WDx(1), WEx(1)
     dimension :: hy(1), WAy(1), WBy(1), WCy(1), WDy(1), WEy(1)
     dimension :: hz(1), WAz(1), WBz(1), WCz(1), WDz(1), WEz(1)
     real*8, allocatable, dimension(:) :: Bsum
     Nxyz = Nx*Ny*Nz
     allocate (Bsum(Nxyz))
!-------------------------------------------------------------------------------
!  B = Bx+By+Bz = dAx/dx+dAy/dy+dAz/dz
!  Bx = dAx/dx
!
     if (ixp .eq. 1) then
!      Periodic Boundary Condition along x-direction
       call PCDx3rd3D_MPI(Ax,WAx,WBx,WCx,WDx,WEx,B,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)
     else
!      Free Boundary Condition along x-direction
       call FCDx3rd3D_MPI(Ax,WAx,WBx,WCx,WDx,WEx,B,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)
     endif
!-------------------------------------------------------------------------------
!  B = Bx+By+Bz = dAx/dx+dAy/dy+dAz/dz
!  Bx+By = dAx/dx+dAy/dy
!
     if (Ny == 1 .and. Nz == 1) return
!
     if (iyp .eq. 1) then 
!      Periodic Boundary Condition along y-direction
       call PCDy3rd3D_MPI(Ay,WAy,WBy,WCy,WDy,WEy,Bsum,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)
     else 
!      Non-Periodic Boundary Condition along y-direction
       call FCDy3rd3D_MPI(Ay,WAy,WBy,WCy,WDy,WEy,Bsum,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)
     endif
!
     do ijk = 1, Nxyz
       B(ijk) = DIFAB(B(ijk),-Bsum(ijk))
!      B(ijk) = B(ijk)+Bsum(ijk)
     enddo
!-------------------------------------------------------------------------------
!  B = Bx+By+Bz = dAx/dx+dAy/dy+dAz/dz
!
!
     if (Nz == 1) return
!
     if (izp .eq. 1) then 
!      Periodic Boundary Condition along z-direction
       call PCDz3rd3D_MPI(Az,WAz,WBz,WCz,WDz,WEz,Bsum,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
     else
!      Non-Periodic Boundary Condition along z-direction
       call FCDz3rd3D_MPI(Az,WAz,WBz,WCz,WDz,WEz,Bsum,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
     endif
!
     do ijk = 1, Nxyz
       B(ijk) = DIFAB(B(ijk),-Bsum(ijk))
!      B(ijk) = B(ijk)+Bsum(ijk)
     enddo
!
     deallocate (Bsum)
   end subroutine div3d_3rd
!===============================================================================
!===============================================================================
!  sub. curl3d_1st
   subroutine curl3d_1st(Nx,Ny,Nz,Ax,Ay,Az,Bx,By,Bz, &
                         hx,WAx,WBx,WCx,hy,WAy,WBy,WCy,hz,WAz,WBz,WCz,ixp,iyp,izp)
!-------------------------------------------------------------------------------
!  B = curl A  using 1st-order finite difference
!
!  Input Arrays   : Ax(Nx,Ny,Nz), Ay(Nx,Ny,Nz), Az(Nx,Ny,Nz)
!                   hx(Nx), WAx(Nx), WBx(Nx), WCx(Nx) 
!                   hy(Ny), WAy(Ny), WBy(Ny), WCy(Ny) 
!                   hz(Nz), WAz(Nz), WBz(Nz), WCz(Nz)
!  Output Arrays  : Bx(Nx,Ny,Nz), By(Nx,Ny,Nz), Bz(Nx,Ny,Nz)
!
!  ixp, iyp, izp \= 1 : Non-Periodic Boundary Condition
!  ixp, iyp, izp = 1  : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: Ax(1), Ay(1), Az(1), Bx(1), By(1), Bz(1)
     dimension :: hx(1), WAx(1), WBx(1), WCx(1)
     dimension :: hy(1), WAy(1), WBy(1), WCy(1)
     dimension :: hz(1), WAz(1), WBz(1), WCz(1)
     real*8, allocatable, dimension(:) :: Bx_sum, By_sum, Bz_sum
!
     Nxyz = Nx*Ny*Nz
     allocate (Bx_sum(Nxyz), By_sum(Nxyz), Bz_sum(Nxyz))
!-------------------------------------------------------------------------------
!  Bx = dAz/dy-dAy/dz
!  By = dAx/dz-dAz/dx
!  Bz = dAy/dx-dAx/dy
!-------------------------------------------------------------------------------
!  dAy/dx and dAz/dx
!
     if (ixp .eq. 1) then 
!  Periodic Boundary Condition along x-direction
!  dAy/dx
       call PCDx1st3D_MPI(Ay,WAx,WBx,WCx,Bz_sum,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)  
!  dAz/dx
       call PCDx1st3D_MPI(Az,WAx,WBx,WCx,By_sum,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)
     else
!  Free Boundary Condition along x-direction
!  dAy/dx
       call FCDx1st3D_MPI(Ay,WAx,WBx,WCx,Bz_sum,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)
!  dAz/dx
       call FCDx1st3D_MPI(Az,WAx,WBx,WCx,By_sum,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)
     endif
!-------------------------------------------------------------------------------
!  dAx/dy and dAz/dy
!
     if (Ny == 1 .and. Nz == 1) then
     	 do ijk = 1, Nxyz
         Bx(ijk)     = 0.d0
         Bx_sum(ijk) = 0.d0
         By(ijk)     = 0.d0
         Bz(ijk)     = 0.d0
       enddo
!
       goto 999
     endif
!
     if (iyp .eq. 1) then 
!  Periodic Boundary Condition along y-direction
!  dAx/dy
       call PCDy1st3D_MPI(Ax,WAy,WBy,WCy,Bz,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)      
!  dAz/dy
       call PCDy1st3D_MPI(Az,WAy,WBy,WCy,Bx_sum,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)
	   else 
!  Free Boundary Condition along y-direction
!  dAx/dy
       call FCDy1st3D_MPI(Ax,WAy,WBy,WCy,Bz,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)
!  dAz/dy
       call FCDy1st3D_MPI(Az,WAy,WBy,WCy,Bx_sum,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)
     endif
!-------------------------------------------------------------------------------
!  dAx/dz and dAy/dz
!
     if (Nz == 1) then
     	 do ijk = 1, Nxyz
         Bx(ijk) = 0.d0
       enddo
!
       do ijk = 1, Nxyz
         By(ijk) = 0.d0
       enddo
!
       goto 999
     endif
!
     if (izp .eq. 1) then 
!  Periodic Boundary Condition along z-direction
!  dAx/dz
       call PCDz1st3D_MPI(Ax,WAz,WBz,WCz,By,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
!  dAy/dz
       call PCDz1st3D_MPI(Ay,WAz,WBz,WCz,Bx,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
     else 
!  Free Boundary Condition along z-direction
!  dAx/dz
       call FCDz1st3D_MPI(Ax,WAz,WBz,WCz,By,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
!  dAy/dz
       call FCDz1st3D_MPI(Ay,WAz,WBz,WCz,Bx,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
     endif
!
 999 continue
!-------------------------------------------------------------------------------
!  Bx = dAz/dy-dAy/dz
!
     do ijk = 1, Nxyz
       Bx(ijk) = DIFAB(Bx_sum(ijk),Bx(ijk))
!      Bx(ijk) = Bx_sum(ijk)-Bx(ijk)
     enddo
!-------------------------------------------------------------------------------
!  By = dAx/dz-dAz/dx
!
     do ijk = 1, Nxyz
       By(ijk) = DIFAB(By(ijk),By_sum(ijk))
!      By(ijk) = By(ijk)-By_sum(ijk)
     enddo
!-------------------------------------------------------------------------------
!  Bz = dAy/dx-dAx/dy
!
     do ijk = 1, Nxyz
       Bz(ijk) = DIFAB(Bz_sum(ijk),Bz(ijk))
!      Bz(ijk) = Bz_sum(ijk)-Bz(ijk)
     enddo
!-------------------------------------------------------------------------------
     deallocate (Bx_sum, By_sum, Bz_sum)
   end subroutine curl3d_1st
!===============================================================================
!===============================================================================
!  sub. curl3d_3rd
   subroutine curl3d_3rd(Nx,Ny,Nz,Ax,Ay,Az,Bx,By,Bz, &
                         hx,WAx,WBx,WCx,WDx,WEx, &
                         hy,WAy,WBy,WCy,WDy,WEy, &
                         hz,WAz,WBz,WCz,WDz,WEz, &
                         ixp,iyp,izp)
!-------------------------------------------------------------------------------
!  B = curl A  using 3rd-order finite difference
!
!  Input Arrays   : Ax(Nx,Ny,Nz), Ay(Nx,Ny,Nz), Az(Nx,Ny,Nz)
!                   hx(Nx), WAx(Nx), WBx(Nx), WCx(Nx), WDx(Nx), WEx(Nx), 
!                   hy(Ny), WAy(Ny), WBy(Ny), WCy(Ny), WDy(Ny), WEy(Ny), 
!                   hz(Nz), WAz(Nz), WBz(Nz), WCz(Nz), WDz(Nz), WEz(Nz)
!  Output Arrays  : Bx(Nx,Ny,Nz), By(Nx,Ny,Nz), Bz(Nx,Ny,Nz)
!
!  ixp, iyp, izp \= 1 : Non-Periodic Boundary Condition
!  ixp, iyp, izp = 1  : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: Ax(1), Ay(1), Az(1), Bx(1), By(1), Bz(1)
     dimension :: hx(1), WAx(1), WBx(1), WCx(1), WDx(1), WEx(1)
     dimension :: hy(1), WAy(1), WBy(1), WCy(1), WDy(1), WEy(1)
     dimension :: hz(1), WAz(1), WBz(1), WCz(1), WDz(1), WEz(1)
     real*8, allocatable, dimension(:) :: Bx_sum, By_sum, Bz_sum
!
     Nxyz = Nx*Ny*Nz
     allocate (Bx_sum(Nxyz), By_sum(Nxyz), Bz_sum(Nxyz))
!-------------------------------------------------------------------------------
!  Bx = dAz/dy-dAy/dz
!  By = dAx/dz-dAz/dx
!  Bz = dAy/dx-dAx/dy
!-------------------------------------------------------------------------------
!  dAy/dx and dAz/dx
!
     if (ixp .eq. 1) then 
!  Periodic Boundary Condition along x-direction
!  dAy/dx
       call PCDx3rd3D_MPI(Ay,WAx,WBx,WCx,WDx,WEx,Bz_sum,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)  
!  dAz/dx
       call PCDx3rd3D_MPI(Az,WAx,WBx,WCx,WDx,WEx,By_sum,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)
     else
!  Free Boundary Condition along x-direction
!  dAy/dx
       call FCDx3rd3D_MPI(Ay,WAx,WBx,WCx,WDx,WEx,Bz_sum,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)
!  dAz/dx
       call FCDx3rd3D_MPI(Az,WAx,WBx,WCx,WDx,WEx,By_sum,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)
     endif
!-------------------------------------------------------------------------------
!  dAx/dy and dAz/dy
!
     if (Ny == 1 .and. Nz == 1) then
     	 do ijk = 1, Nxyz
         Bx(ijk)     = 0.d0
         Bx_sum(ijk) = 0.d0
         By(ijk)     = 0.d0
         Bz(ijk)     = 0.d0
       enddo
!
       goto 999
     endif
!
     if (iyp .eq. 1) then 
!  Periodic Boundary Condition along y-direction
!  dAx/dy
       call PCDy3rd3D_MPI(Ax,WAy,WBy,WCy,WDy,WEy,Bz,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)      
!  dAz/dy
       call PCDy3rd3D_MPI(Az,WAy,WBy,WCy,WDy,WEy,Bx_sum,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)
	   else 
!  Free Boundary Condition along y-direction
!  dAx/dy
       call FCDy3rd3D_MPI(Ax,WAy,WBy,WCy,WDy,WEy,Bz,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)
!  dAz/dy
       call FCDy3rd3D_MPI(Az,WAy,WBy,WCy,WDy,WEy,Bx_sum,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)
     endif
!-------------------------------------------------------------------------------
!  dAx/dz and dAy/dz
!
     if (Nz == 1) then
     	 do ijk = 1, Nxyz
         Bx(ijk) = 0.d0
       enddo
!
       do ijk = 1, Nxyz
         By(ijk) = 0.d0
       enddo
!
       goto 999
     endif
!
     if (izp .eq. 1) then 
!  Periodic Boundary Condition along z-direction
!  dAx/dz
       call PCDz3rd3D_MPI(Ax,WAz,WBz,WCz,WDz,WEz,By,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
!  dAy/dz
       call PCDz3rd3D_MPI(Ay,WAz,WBz,WCz,WDz,WEz,Bx,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
     else 
!  Free Boundary Condition along z-direction
!  dAx/dz
       call FCDz3rd3D_MPI(Ax,WAz,WBz,WCz,WDz,WEz,By,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
!  dAy/dz
       call FCDz3rd3D_MPI(Ay,WAz,WBz,WCz,WDz,WEz,Bx,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
     endif
!
 999 continue
!-------------------------------------------------------------------------------
!  Bx = dAz/dy-dAy/dz
!
     do ijk = 1, Nxyz
       Bx(ijk) = DIFAB(Bx_sum(ijk),Bx(ijk))
!      Bx(ijk) = Bx_sum(ijk)-Bx(ijk)
     enddo
!-------------------------------------------------------------------------------
!  By = dAx/dz-dAz/dx
!
     do ijk = 1, Nxyz
       By(ijk) = DIFAB(By(ijk),By_sum(ijk))
!      By(ijk) = By(ijk)-By_sum(ijk)
     enddo
!-------------------------------------------------------------------------------
!  Bz = dAy/dx-dAx/dy
!
     do ijk = 1, Nxyz
       Bz(ijk) = DIFAB(Bz_sum(ijk),Bz(ijk))
!      Bz(ijk) = Bz_sum(ijk)-Bz(ijk)
     enddo
!-------------------------------------------------------------------------------
     deallocate (Bx_sum, By_sum, Bz_sum)
   end subroutine curl3d_3rd
!===============================================================================
!===============================================================================
!  sub. divgrad3d_1st
   subroutine divgrad3d_1st(Nx,Ny,Nz,A,B, &
                            hx,WAx,WBx,WCx,hy,WAy,WBy,WCy,hz,WAz,WBz,WCz, &
                            ixp,iyp,izp)
!-------------------------------------------------------------------------------
!  B = div ( grad A )  using 1st-order finite difference
!
!  Input Arrays   : A(Nx,Ny,Nz)
!                   hx(Nx), WAx(Nx), WBx(Nx), WCx(Nx), 
!                   hy(Ny), WAy(Ny), WBy(Ny), WCy(Ny), 
!                   hz(Nz), WAz(Nz), WBz(Nz), WCz(Nz) 
!  Output Arrays  : B(Nx,Ny,Nz)
!
!  ixp, iyp, izp \= 1 : Free Boundary Condition
!  ixp, iyp, izp = 1  : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: A(1), B(1)
     dimension :: hx(1), WAx(1), WBx(1), WCx(1)
     dimension :: hy(1), WAy(1), WBy(1), WCy(1)
     dimension :: hz(1), WAz(1), WBz(1), WCz(1)
     real*8, allocatable, dimension(:) :: Bsum, Bp
     Nxyz = Nx*Ny*Nz
     allocate (Bsum(Nxyz), Bp(Nxyz))
!-------------------------------------------------------------------------------
!  B = ddA/dxx+ddA/dyy+ddA/dzz
!  ddA/dxx
!
     if (ixp .eq. 1) then 
!  Periodic Boundary Condition along x-direction
       call PCDx1st3D_MPI(A,WAx,WBx,WCx,Bp,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)
       call PCDx1st3D_MPI(Bp,WAx,WBx,WCx,B,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)                          
     else 
!  Free Boundary Condition along x-direction
       call FCDx1st3D_MPI(A,WAx,WBx,WCx,Bp,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)
       call FCDx1st3D_MPI(Bp,WAx,WBx,WCx,B,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)                            
     endif
!-------------------------------------------------------------------------------
!  ddA/dyy
!
     if (Ny == 1 .and. Nz == 1) return
!
     if (iyp .eq. 1) then 
!  Periodic Boundary Condition along y-direction
       call PCDy1st3D_MPI(A,WAy,WBy,WCy,Bp,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)
       call PCDy1st3D_MPI(Bp,WAy,WBy,WCy,Bsum,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)                          
     else 
!  Free Boundary Condition along y-direction
       call FCDy1st3D_MPI(A,WAy,WBy,WCy,Bp,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)
       call FCDy1st3D_MPI(Bp,WAy,WBy,WCy,Bsum,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)                           
     endif
!-------------------------------------------------------------------------------
!  ddA/dxx+ddA/dyy
!
     do ijk = 1, Nxyz
       B(ijk) = DIFAB(B(ijk),-Bsum(ijk))
!      B(ijk) = B(ijk)+Bsum(ijk)
     enddo
!-------------------------------------------------------------------------------
!  ddA/dzz
!
     if (Nz == 1) return
!
     if (izp .eq. 1) then 
!  Periodic Boundary Condition along z-direction
       call PCDz1st3D_MPI(A,WAz,WBz,WCz,Bp,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
       call PCDz1st3D_MPI(Bp,WAz,WBz,WCz,Bsum,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)                          
     else 
!  Free Boundary Condition along z-direction
       call FCDz1st3D_MPI(A,WAz,WBz,WCz,Bp,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
       call FCDz1st3D_MPI(Bp,WAz,WBz,WCz,Bsum,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
     endif
!-------------------------------------------------------------------------------
!  B = ddA/dxx+ddA/dyy+ddA/dzz
!
     do ijk = 1, Nxyz
       B(ijk) = DIFAB(B(ijk),-Bsum(ijk))
!      B(ijk) = B(ijk)+Bsum(ijk)
     enddo
!-------------------------------------------------------------------------------
     deallocate (Bsum, Bp)
   end subroutine divgrad3d_1st
!===============================================================================
!===============================================================================
!  sub. divgrad3d_3rd
   subroutine divgrad3d_3rd(Nx,Ny,Nz,A,B, &
                            hx,WAx,WBx,WCx,WDx,WEx, &
                            hy,WAy,WBy,WCy,WDy,WEy, &
                            hz,WAz,WBz,WCz,WDz,WEz, &
                            ixp,iyp,izp)
!-------------------------------------------------------------------------------
!  B = div ( grad A )  using 3rd-order finite difference
!
!  Input Arrays   : A(Nx,Ny,Nz)
!                   hx(Nx), WAx(Nx), WBx(Nx), WCx(Nx), WDx(Nx), WEx(Nx), 
!                   hy(Ny), WAy(Ny), WBy(Ny), WCy(Ny), WDy(Ny), WEy(Ny), 
!                   hz(Nz), WAz(Nz), WBz(Nz), WCz(Nz), WDz(Nz), WEz(Nz) 
!  Output Arrays  : B(Nx,Ny,Nz)
!
!  ixp, iyp, izp \= 1 : Free Boundary Condition
!  ixp, iyp, izp = 1  : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: A(1), B(1)
     dimension :: hx(1), WAx(1), WBx(1), WCx(1), WDx(1), WEx(1)
     dimension :: hy(1), WAy(1), WBy(1), WCy(1), WDy(1), WEy(1)
     dimension :: hz(1), WAz(1), WBz(1), WCz(1), WDz(1), WEz(1)
     real*8, allocatable, dimension(:) :: Bsum, Bp
     Nxyz = Nx*Ny*Nz
     allocate (Bsum(Nxyz), Bp(Nxyz))
!-------------------------------------------------------------------------------
!  B = ddA/dxx+ddA/dyy+ddA/dzz
!  ddA/dxx
!
     if (ixp .eq. 1) then 
!  Periodic Boundary Condition along x-direction
       call PCDx3rd3D_MPI(A,WAx,WBx,WCx,WDx,WEx,Bp,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)
       call PCDx3rd3D_MPI(Bp,WAx,WBx,WCx,WDx,WEx,B,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)                          
     else 
!  Free Boundary Condition along x-direction
       call FCDx3rd3D_MPI(A,WAx,WBx,WCx,WDx,WEx,Bp,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)
       call FCDx3rd3D_MPI(Bp,WAx,WBx,WCx,WDx,WEx,B,Nx,Ny,Nz, &
                          ipm1,MyID,ipp1,ipstart,ipend)                            
     endif
!-------------------------------------------------------------------------------
!  ddA/dyy
!
     if (Ny == 1 .and. Nz == 1) return
!
     if (iyp .eq. 1) then 
!  Periodic Boundary Condition along y-direction
       call PCDy3rd3D_MPI(A,WAy,WBy,WCy,WDy,WEy,Bp,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)
       call PCDy3rd3D_MPI(Bp,WAy,WBy,WCy,WDy,WEy,Bsum,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)                          
     else 
!  Free Boundary Condition along y-direction
       call FCDy3rd3D_MPI(A,WAy,WBy,WCy,WDy,WEy,Bp,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)
       call FCDy3rd3D_MPI(Bp,WAy,WBy,WCy,WDy,WEy,Bsum,Nx,Ny,Nz, &
                          iqm1,MyID,iqp1,iqstart,iqend)                           
     endif
!-------------------------------------------------------------------------------
!  ddA/dxx+ddA/dyy
!
     do ijk = 1, Nxyz
       B(ijk) = DIFAB(B(ijk),-Bsum(ijk))
!      B(ijk) = B(ijk)+Bsum(ijk)
     enddo
!-------------------------------------------------------------------------------
!  ddA/dzz
!
     if (Nz == 1) return
!
     if (izp .eq. 1) then 
!  Periodic Boundary Condition along z-direction
       call PCDz3rd3D_MPI(A,WAz,WBz,WCz,WDz,WEz,Bp,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
       call PCDz3rd3D_MPI(Bp,WAz,WBz,WCz,WDz,WEz,Bsum,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)                          
     else 
!  Free Boundary Condition along z-direction
       call FCDz3rd3D_MPI(A,WAz,WBz,WCz,WDz,WEz,Bp,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
       call FCDz3rd3D_MPI(Bp,WAz,WBz,WCz,WDz,WEz,Bsum,Nx,Ny,Nz, &
                          irm1,MyID,irp1,irstart,irend)
     endif
!-------------------------------------------------------------------------------
!  B = ddA/dxx+ddA/dyy+ddA/dzz
!
     do ijk = 1, Nxyz
       B(ijk) = DIFAB(B(ijk),-Bsum(ijk))
!      B(ijk) = B(ijk)+Bsum(ijk)
     enddo
!-------------------------------------------------------------------------------
     deallocate (Bsum, Bp)
   end subroutine divgrad3d_3rd
!===============================================================================
!===============================================================================
!  sub. INTP3dm_5th
   subroutine INTP3dm_5th(Nx,Ny,Nz,Ax,Ay,Az,Bx,By,Bz,ixp,iyp,izp)
!-------------------------------------------------------------------------------
!  5th-order interpolation method
!    B(i) = (3*A(i-3)-25*A(i-2)+150*A(i-1)+150*A(i)-25*A(i+1)+3*A(i))/256
!
!  Input Arrays   : Ax(Nx,Ny,Nz), Ay(Nx,Ny,Nz), Az(Nx,Ny,Nz)
!
!  Output Arrays  : Bx(Nx,Ny,Nz), By(Nx,Ny,Nz), Bz(Nx,Ny,Nz)
!
!  ixp, iyp, izp \= 1 : Free Boundary Condition
!  ixp, iyp, izp =  1 : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: Ax(1), Ay(1), Az(1)
     dimension :: Bx(1), By(1), Bz(1)
!     Nxyz = Nx*Ny*Nz
!-------------------------------------------------------------------------------
!  interpolation along x-direction
!
     if (ixp .eq. 1) then
!      Periodic Boundary Condition along x-direction
       call PINTPxm5th3D_MPI(Ax,Bx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)
     else
!      Free Boundary Condition along x-direction
       call FINTPxm5th3D_MPI(Ax,Bx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)
     endif
!-------------------------------------------------------------------------------
! 
!  By = dA/dy
!
     if (Ny == 1 .and. Nz == 1) then
     	 By(1) = 0.d0
     	 Bz(1) = 0.d0
     	 return
     endif
!
     if (iyp .eq. 1) then 
!      Periodic Boundary Condition along y-direction
       call PINTPym5th3D_MPI(Ay,By,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)
     else 
!      Non-Periodic Boundary Condition along y-direction
       call FINTPym5th3D_MPI(Ay,By,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)
     endif
!-------------------------------------------------------------------------------
! 
!  Bz = dA/dz
!
!
     if (Nz == 1) then
     	 Bz(1) = 0.d0
     	 return
     endif
!
     if (izp .eq. 1) then 
!      Periodic Boundary Condition along z-direction
       call PINTPzm5th3D_MPI(Az,Bz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)
     else
!      Non-Periodic Boundary Condition along z-direction
       call FINTPzm5th3D_MPI(Az,Bz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)
     endif
!
   end subroutine INTP3dm_5th
!===============================================================================
!===============================================================================
!  sub. INTP3dp_5th
   subroutine INTP3dp_5th(Nx,Ny,Nz,Ax,Ay,Az,Bx,By,Bz,ixp,iyp,izp)
!-------------------------------------------------------------------------------
!  5th-order interpolation method
!    B(i) = (3*A(i-3)-25*A(i-2)+150*A(i-1)+150*A(i)-25*A(i+1)+3*A(i))/256
!
!  Input Arrays   : Ax(Nx,Ny,Nz), Ay(Nx,Ny,Nz), Az(Nx,Ny,Nz)
!
!  Output Arrays  : Bx(Nx,Ny,Nz), By(Nx,Ny,Nz), Bz(Nx,Ny,Nz)
!
!  ixp, iyp, izp \= 1 : Free Boundary Condition
!  ixp, iyp, izp =  1 : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: Ax(1), Ay(1), Az(1)
     dimension :: Bx(1), By(1), Bz(1)
!     Nxyz = Nx*Ny*Nz
!-------------------------------------------------------------------------------
!  interpolation along x-direction
!
     if (ixp .eq. 1) then
!      Periodic Boundary Condition along x-direction
       call PINTPxp5th3D_MPI(Ax,Bx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)
     else
!      Free Boundary Condition along x-direction
       call FINTPxp5th3D_MPI(Ax,Bx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)
     endif
!-------------------------------------------------------------------------------
! 
!  By = dA/dy
!
     if (Ny == 1 .and. Nz == 1) then
     	 By(1) = 0.d0
     	 Bz(1) = 0.d0
     	 return
     endif
!
     if (iyp .eq. 1) then 
!      Periodic Boundary Condition along y-direction
       call PINTPyp5th3D_MPI(Ay,By,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)
     else 
!      Non-Periodic Boundary Condition along y-direction
       call FINTPyp5th3D_MPI(Ay,By,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)
     endif
!-------------------------------------------------------------------------------
! 
!  Bz = dA/dz
!
!
     if (Nz == 1) then
     	 Bz(1) = 0.d0
     	 return
     endif
!
     if (izp .eq. 1) then 
!      Periodic Boundary Condition along z-direction
       call PINTPzp5th3D_MPI(Az,Bz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)
     else
!      Non-Periodic Boundary Condition along z-direction
       call FINTPzp5th3D_MPI(Az,Bz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)
     endif
!
   end subroutine INTP3dp_5th
!===============================================================================
!===============================================================================
!  sub. grad3d_5th
   subroutine grad3d_5th(Nx,Ny,Nz,A,Bx,By,Bz,dx,dy,dz,ixp,iyp,izp)
!-------------------------------------------------------------------------------
!  B = gradient A   use 5th-order central difference
!
!  Input Arrays   : A(Nx,Ny,Nz)
!  Input datas    : dx, dy, dz
!
!  Output Arrays  : Bx(Nx,Ny,Nz), By(Nx,Ny,Nz), Bz(Nx,Ny,Nz)
!
!  ixp, iyp, izp \= 1 : Free Boundary Condition
!  ixp, iyp, izp =  1 : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: A(1), Bx(1), By(1), Bz(1)
     Nxyz = Nx*Ny*Nz
!-------------------------------------------------------------------------------
! 
!  Bx = dA/dx
!
     if (ixp .eq. 1) then
!      Periodic Boundary Condition along x-direction
       call PCDx5th3D_MPI(A,Bx,dx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)
!       call WENOx5th6D_MPI(A,Bx,dx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)
     else
!      Free Boundary Condition along x-direction
       call FCDx5th3D_MPI(A,Bx,dx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)
     endif
!-------------------------------------------------------------------------------
! 
!  By = dA/dy
!
     if (Ny == 1 .and. Nz == 1) then
     	 By(1) = 0.d0
     	 Bz(1) = 0.d0
     	 return
     endif
!
     if (iyp .eq. 1) then 
!      Periodic Boundary Condition along y-direction
       call PCDy5th3D_MPI(A,By,dy,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)
     else 
!      Non-Periodic Boundary Condition along y-direction
       call FCDy5th3D_MPI(A,By,dy,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)
     endif
!-------------------------------------------------------------------------------
! 
!  Bz = dA/dz
!
!
     if (Nz == 1) then
     	 Bz(1) = 0.d0
     	 return
     endif
!
     if (izp .eq. 1) then 
!      Periodic Boundary Condition along z-direction
       call PCDz5th3D_MPI(A,Bz,dz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)
     else
!      Non-Periodic Boundary Condition along z-direction
       call FCDz5th3D_MPI(A,Bz,dz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)
     endif
!
   end subroutine grad3d_5th
!===============================================================================
!===============================================================================
!  sub. gradu_5th
   subroutine gradu_5th(Nx,Ny,Nz,A,Bx,By,Bz,dx,dy,dz)
!-------------------------------------------------------------------------------
!  B = gradient A   use 5th-order central difference no MPI version
!                   for free boundary condition
!
!  Input Arrays   : A(Nx,Ny,Nz)
!  Input datas    : dx, dy, dz
!
!  Output Arrays  : Bx(Nx,Ny,Nz), By(Nx,Ny,Nz), Bz(Nx,Ny,Nz)
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: A(1), Bx(1), By(1), Bz(1)
     Nxyz = Nx*Ny*Nz
!-------------------------------------------------------------------------------
! 
!  Bx = dA/dx
!
     call FCDx5th3D(A,Bx,dx,Nx,Ny,Nz)
!-------------------------------------------------------------------------------
! 
!  By = dA/dy
!
     if (Ny == 1 .and. Nz == 1) then
     	 By(1) = 0.d0
     	 Bz(1) = 0.d0
     	 return
     endif
!
     call FCDy5th3D(A,By,dy,Nx,Ny,Nz)
!-------------------------------------------------------------------------------
! 
!  Bz = dA/dz
!
!
     if (Nz == 1) then
     	 Bz(1) = 0.d0
     	 return
     endif
!
     call FCDz5th3D(A,Bz,dz,Nx,Ny,Nz)
!
   end subroutine gradu_5th
!===============================================================================
!===============================================================================
!  sub. gradu_5th
   subroutine gradu_10th(Nx,Ny,Nz,A,Bx,By,Bz,dx,dy,dz)
!-------------------------------------------------------------------------------
!  B = gradient A   use 5th-order central difference no MPI version
!                   for free boundary condition
!
!  Input Arrays   : A(Nx,Ny,Nz)
!  Input datas    : dx, dy, dz
!
!  Output Arrays  : Bx(Nx,Ny,Nz), By(Nx,Ny,Nz), Bz(Nx,Ny,Nz)
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: A(1), Bx(1), By(1), Bz(1)
     Nxyz = Nx*Ny*Nz
!-------------------------------------------------------------------------------
! 
!  Bx = dA/dx
!
     call FCDx10th3D(A,Bx,dx,Nx,Ny,Nz)
!-------------------------------------------------------------------------------
! 
!  By = dA/dy
!
     if (Ny == 1 .and. Nz == 1) then
     	 By(1) = 0.d0
     	 Bz(1) = 0.d0
     	 return
     endif
!
     call FCDy10th3D(A,By,dy,Nx,Ny,Nz)
!-------------------------------------------------------------------------------
! 
!  Bz = dA/dz
!
!
     if (Nz == 1) then
     	 Bz(1) = 0.d0
     	 return
     endif
!
     call FCDz10th3D(A,Bz,dz,Nx,Ny,Nz)
!
   end subroutine gradu_10th
!===============================================================================
!===============================================================================
!--SUB. dfdvx in 3D velocity space
   subroutine gradux_CWENO(f,fp,du,nux,nuy,nuz)
!
!  calculate the first derivative of the distribution function
!    along x-axis in 3D velocity space
!
!  input array   : f(nux,nuy,nuz)
!  input array   : Bx(nx), Cx(nx), hx(nx)
!  working array : Rx(nx), FPx(nx), FP0x(nx)
!  output array  : fp(nux,nuy,nuz)
!
     implicit double precision (A-H,O-Z)

     dimension :: f(1), fp(1)
     nuxy = nux*nuy
!-------------------------------------------------------------------------------
     call CWENOux6th(f,fp,du,nux,nuy,nuz)
!

!-------------------------------------------------------------------------------  
   end subroutine gradux_CWENO
!===============================================================================
!===============================================================================
!--SUB. dfdvx in 3D velocity space
   subroutine gradux_WENO5(f,fp,du,nux,nuy,nuz,E)
!
!  calculate the first derivative of the distribution function
!    along x-axis in 3D velocity space
!
!  input array   : f(nux,nuy,nuz)
!  input array   : Bx(nx), Cx(nx), hx(nx)
!  working array : Rx(nx), FPx(nx), FP0x(nx)
!  output array  : fp(nux,nuy,nuz)
!
     implicit double precision (A-H,O-Z)

     dimension :: f(1), fp(1)
     nuxy = nux*nuy
!-------------------------------------------------------------------------------
     call WENOux5th(f,fp,du,nux,nuy,nuz,E)
!

!-------------------------------------------------------------------------------  
   end subroutine gradux_WENO5
!===============================================================================
!===============================================================================
!--SUB. dfdvx in 3D velocity space
   subroutine FP_ux_CWENO(f,fp,du,nux,nuy,nuz,vx)
!
!  calculate the first derivative of the distribution function
!    along x-axis in 3D velocity space
!
!  input array   : f(nux,nuy,nuz)
!  input array   : Bx(nx), Cx(nx), hx(nx)
!  working array : Rx(nx), FPx(nx), FP0x(nx)
!  output array  : fp(nux,nuy,nuz)
!
     implicit double precision (A-H,O-Z)

     dimension :: f(1), fp(1),vx(1)
     real(8), allocatable, dimension(:) :: wkfx, wkfy, wkfz
     nuxyz = nux*nuy*nuz
     allocate (wkfx(nuxyz), wkfy(nuxyz), wkfz(nuxyz))
     wkfx = 0.d0
     wkfy = 0.d0
     wkfz = 0.d0
!-------------------------------------------------------------------------------
     call CWENOux6th(f,wkfx,du,nux,nuy,nuz)
     do iu =1,nux
        wkfx(iu) = vx(iu)*f(iu)+wkfx(iu)
     enddo
     call CWENOux6th(wkfx,fp,du,nux,nuy,nuz)

!------------------------------------------------------------------------------- 
   deallocate (wkfx, wkfy, wkfz) 
   end subroutine FP_ux_CWENO
!===============================================================================
!===============================================================================
!  sub. div3d_5th
   subroutine div3d_5th(Nx,Ny,Nz,Ax,Ay,Az,B,dx,dy,dz,ixp,iyp,izp)
!-------------------------------------------------------------------------------
!  B = divergence A  use 5th-order central difference
!
!  Input Arrays   : Ax(Nx,Ny,Nz), Ay(Nx,Ny,Nz), Az(Nx,Ny,Nz)
!
!  Output Arrays  : B(Nx,Ny,Nz)
!
!  ixp, iyp, izp \= 1 : Free Boundary Condition
!  ixp, iyp, izp =  1 : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: Ax(1), Ay(1), Az(1), B(1)
     real*8, allocatable, dimension(:) :: Bsum
     Nxyz = Nx*Ny*Nz
     allocate (Bsum(Nxyz))
!-------------------------------------------------------------------------------
!  B = Bx+By+Bz = dAx/dx+dAy/dy+dAz/dz
!  Bx = dAx/dx
!
     if (ixp .eq. 1) then
!      Periodic Boundary Condition along x-direction
       call PCDx5th3D_MPI(Ax,B,dx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)
     else
!      Free Boundary Condition along x-direction
       call FCDx5th3D_MPI(Ax,B,dx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)
     endif
!-------------------------------------------------------------------------------
!  B = Bx+By+Bz = dAx/dx+dAy/dy+dAz/dz
!  Bx+By = dAx/dx+dAy/dy
!
     if (Ny == 1 .and. Nz == 1) return
!
     if (iyp .eq. 1) then 
!      Periodic Boundary Condition along y-direction
       call PCDy5th3D_MPI(Ay,Bsum,dy,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)
     else 
!      Non-Periodic Boundary Condition along y-direction
       call FCDy5th3D_MPI(Ay,Bsum,dy,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)
     endif
!
!$OMP parallel do
     do ijk = 1, Nxyz
       B(ijk) = DIFAB(B(ijk),-Bsum(ijk))
!      B(ijk) = B(ijk)+Bsum(ijk)
     enddo
!$OMP end parallel do
!-------------------------------------------------------------------------------
!  B = Bx+By+Bz = dAx/dx+dAy/dy+dAz/dz
!
!
     if (Nz == 1) return
!
     if (izp .eq. 1) then 
!      Periodic Boundary Condition along z-direction
       call PCDz5th3D_MPI(Az,Bsum,dz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)
     else
!      Non-Periodic Boundary Condition along z-direction
       call FCDz5th3D_MPI(Az,Bsum,dz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)
     endif
!
!$OMP parallel do
     do ijk = 1, Nxyz
       B(ijk) = DIFAB(B(ijk),-Bsum(ijk))
!      B(ijk) = B(ijk)+Bsum(ijk)
     enddo
!$OMP end parallel do
!
     deallocate (Bsum)
   end subroutine div3d_5th
!===============================================================================
!===============================================================================
!  sub. curl3d_5th
   subroutine curl3d_5th(Nx,Ny,Nz,Ax,Ay,Az,Bx,By,Bz,dx,dy,dz,ixp,iyp,izp)
!-------------------------------------------------------------------------------
!  B = curl A using 5th-order finite central difference
!
!  Input Arrays  : Ax(Nx,Ny,Nz), Ay(Nx,Ny,Nz), Az(Nx,Ny,Nz)
!  Input datas   : dx, dy, dz
!
!  Output Arrays : Bx(Nx,Ny,Nz), By(Nx,Ny,Nz), Bz(Nx,Ny,Nz)
!
!  ixp, iyp, izp \= 1 : Non-Periodic Boundary Condition
!  ixp, iyp, izp = 1  : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: Ax(1), Ay(1), Az(1), Bx(1), By(1), Bz(1)
     real*8, allocatable, dimension(:) :: Bx_sum, By_sum, Bz_sum
!
     Nxyz = Nx*Ny*Nz
     allocate (Bx_sum(Nxyz), By_sum(Nxyz), Bz_sum(Nxyz))
!-------------------------------------------------------------------------------
!  Bx = dAz/dy-dAy/dz
!  By = dAx/dz-dAz/dx
!  Bz = dAy/dx-dAx/dy
!-------------------------------------------------------------------------------
!  dAy/dx and dAz/dx
!
     if (ixp .eq. 1) then 
!  Periodic Boundary Condition along x-direction
!  dAy/dx
       call PCDx5th3D_MPI(Ay,Bz_sum,dx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)  
!  dAz/dx
       call PCDx5th3D_MPI(Az,By_sum,dx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)
     else
!  Free Boundary Condition along x-direction
!  dAy/dx
       call FCDx5th3D_MPI(Ay,Bz_sum,dx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)
!  dAz/dx
       call FCDx5th3D_MPI(Az,By_sum,dx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)
     endif
!-------------------------------------------------------------------------------
!  dAx/dy and dAz/dy
!
     if (Ny == 1 .and. Nz == 1) then
!       Bx(1:Nxyz)     = 0.d0
!       Bx_sum(1:Nxyz) = 0.d0
!       By(1:Nxyz)     = 0.d0
!       Bz(1:Nxyz)     = 0.d0
     	 do ijk = 1, Nxyz
         Bx(ijk)     = 0.d0
         Bx_sum(ijk) = 0.d0
         By(ijk)     = 0.d0
         Bz(ijk)     = 0.d0
       enddo
!
       goto 999
     endif
!
     if (iyp .eq. 1) then 
!  Periodic Boundary Condition along y-direction
!  dAx/dy
       call PCDy5th3D_MPI(Ax,Bz,dy,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)      
!  dAz/dy
       call PCDy5th3D_MPI(Az,Bx_sum,dy,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)
	   else 
!  Free Boundary Condition along y-direction
!  dAx/dy
       call FCDy5th3D_MPI(Ax,Bz,dy,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)
!  dAz/dy
       call FCDy5th3D_MPI(Az,Bx_sum,dy,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)
     endif
!-------------------------------------------------------------------------------
!  dAx/dz and dAy/dz
!
     if (Nz == 1) then
!       Bx(1:Nxyz) = 0.d0
!       By(1:Nxyz) = 0.d0
     	 do ijk = 1, Nxyz
         Bx(ijk) = 0.d0
       enddo
!
       do ijk = 1, Nxyz
         By(ijk) = 0.d0
       enddo
!
       goto 999
     endif
!
     if (izp .eq. 1) then 
!  Periodic Boundary Condition along z-direction
!  dAx/dz
       call PCDz5th3D_MPI(Ax,By,dz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)
!  dAy/dz
       call PCDz5th3D_MPI(Ay,Bx,dz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)
     else 
!  Free Boundary Condition along z-direction
!  dAx/dz
       call FCDz5th3D_MPI(Ax,By,dz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)
!  dAy/dz
       call FCDz5th3D_MPI(Ay,Bx,dz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)
     endif
!
 999 continue
!-------------------------------------------------------------------------------
!  Bx = dAz/dy-dAy/dz
!  By = dAx/dz-dAz/dx
!  Bz = dAy/dx-dAx/dy
!
!$OMP parallel do
     do ijk = 1, Nxyz
       Bx(ijk) = DIFAB(Bx_sum(ijk),Bx(ijk))
!      Bx(ijk) = Bx_sum(ijk)-Bx(ijk)
       By(ijk) = DIFAB(By(ijk),By_sum(ijk))
!      By(ijk) = By(ijk)-By_sum(ijk)
       Bz(ijk) = DIFAB(Bz_sum(ijk),Bz(ijk))
!      Bz(ijk) = Bz_sum(ijk)-Bz(ijk)
     enddo
!$OMP end parallel do
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
     deallocate (Bx_sum, By_sum, Bz_sum)
   end subroutine curl3d_5th
!===============================================================================
!===============================================================================
!  sub. divgrad3d_5th
   subroutine divgrad3d_5th(Nx,Ny,Nz,A,B,dx,dy,dz,ixp,iyp,izp)
!-------------------------------------------------------------------------------
!  B = div ( grad A )  using 5th-order central finite difference
!
!  Input Arrays  : A(Nx,Ny,Nz)
!  Input datas   : dx, dy, dz
! 
!  Output Arrays  : B(Nx,Ny,Nz)
!
!  ixp, iyp, izp \= 1 : Free Boundary Condition
!  ixp, iyp, izp = 1  : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: A(1), B(1)
     real*8, allocatable, dimension(:) :: Bsum
     Nxyz = Nx*Ny*Nz
     allocate (Bsum(Nxyz))
!-------------------------------------------------------------------------------
!  B = ddA/dxx+ddA/dyy+ddA/dzz
!  ddA/dxx
!
     if (ixp .eq. 1) then 
!  Periodic Boundary Condition along x-direction
       call PCDxx5th3D_MPI(A,B,dx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)                          
     else 
!  Free Boundary Condition along x-direction
       call FCDxx5th3D_MPI(A,B,dx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)
     endif
!-------------------------------------------------------------------------------
!  ddA/dyy
!
     if (Ny == 1 .and. Nz == 1) return
!
     if (iyp .eq. 1) then 
!  Periodic Boundary Condition along y-direction
       call PCDyy5th3D_MPI(A,Bsum,dy,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)                          
     else 
!  Free Boundary Condition along y-direction
       call FCDyy5th3D_MPI(A,Bsum,dy,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)                           
     endif
!-------------------------------------------------------------------------------
!  ddA/dxx+ddA/dyy
!
!$OMP parallel do
     do ijk = 1, Nxyz
       B(ijk) = DIFAB(B(ijk),-Bsum(ijk))
!      B(ijk) = B(ijk)+Bsum(ijk)
     enddo
!$OMP end parallel do
!-------------------------------------------------------------------------------
!  ddA/dzz
!
     if (Nz == 1) return
!
     if (izp .eq. 1) then 
!  Periodic Boundary Condition along z-direction
       call PCDzz5th3D_MPI(A,Bsum,dz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)                          
     else 
!  Free Boundary Condition along z-direction
       call FCDzz5th3D_MPI(A,Bsum,dz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)
     endif
!-------------------------------------------------------------------------------
!  B = ddA/dxx+ddA/dyy+ddA/dzz
!
!$OMP parallel do
     do ijk = 1, Nxyz
       B(ijk) = DIFAB(B(ijk),-Bsum(ijk))
!      B(ijk) = B(ijk)+Bsum(ijk)
     enddo
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (Bsum)
   end subroutine divgrad3d_5th
!===============================================================================
!===============================================================================
!  sub. divgradx3d_5th
   subroutine divgradx3d_5th(Nx,Ny,Nz,A,B,dx,ixp)
!-------------------------------------------------------------------------------
!  B = d2A/dx2 )  using 5th-order central finite difference
!
!  Input Arrays  : A(Nx,Ny,Nz)
!  Input datas   : dx
! 
!  Output Arrays  : B(Nx,Ny,Nz)
!
!  ixp, iyp, izp \= 1 : Free Boundary Condition
!  ixp, iyp, izp = 1  : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: A(1), B(1)
     Nxyz = Nx*Ny*Nz
!-------------------------------------------------------------------------------
!  B = ddA/dxx
!
     if (ixp .eq. 1) then 
!  Periodic Boundary Condition along x-direction
       call PCDxx5th3D_MPI(A,B,dx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend)                          
     else 
!  Free Boundary Condition along x-direction
       call FCDxx5th3D_MPI(A,B,dx,Nx,Ny,Nz,ipm1,MyID,ipp1,ipstart,ipend) 
     endif
!-------------------------------------------------------------------------------
   end subroutine divgradx3d_5th
!===============================================================================
!===============================================================================
!  sub. divgrady3d_5th
   subroutine divgrady3d_5th(Nx,Ny,Nz,A,B,dy,iyp)
!-------------------------------------------------------------------------------
!  B = d2A/dy2 )  using 5th-order central finite difference
!
!  Input Arrays  : A(Nx,Ny,Nz)
!  Input datas   : dy
! 
!  Output Arrays  : B(Nx,Ny,Nz)
!
!  ixp, iyp, izp \= 1 : Free Boundary Condition
!  ixp, iyp, izp = 1  : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: A(1), B(1)
     Nxyz = Nx*Ny*Nz
!-------------------------------------------------------------------------------
!  B = ddA/dyy
!
     B(1:Nxyz) = 0.d0     
     if (Ny == 1) return
!
     if (iyp .eq. 1) then 
!  Periodic Boundary Condition along y-direction
       call PCDyy5th3D_MPI(A,B,dy,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)                         
     else 
!  Free Boundary Condition along y-direction
       call FCDyy5th3D_MPI(A,B,dy,Nx,Ny,Nz,iqm1,MyID,iqp1,iqstart,iqend)                           
     endif
!-------------------------------------------------------------------------------
   end subroutine divgrady3d_5th
!===============================================================================
!===============================================================================
!  sub. divgradz3d_5th
   subroutine divgradz3d_5th(Nx,Ny,Nz,A,B,dz,izp)
!-------------------------------------------------------------------------------
!  B = d2A/dz2  using 5th-order central finite difference
!
!  Input Arrays  : A(Nx,Ny,Nz)
!  Input datas   : dz
! 
!  Output Arrays  : B(Nx,Ny,Nz)
!
!  ixp, iyp, izp \= 1 : Free Boundary Condition
!  ixp, iyp, izp = 1  : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: A(1), B(1)
     Nxyz = Nx*Ny*Nz
!-------------------------------------------------------------------------------
!  B = ddA/dzz
!
     B(1:Nxyz) = 0.d0
     if (Nz == 1) return
!
     if (izp .eq. 1) then 
!  Periodic Boundary Condition along z-direction
       call PCDzz5th3D_MPI(A,B,dz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)                         
     else 
!  Free Boundary Condition along z-direction
       call FCDzz5th3D_MPI(A,B,dz,Nx,Ny,Nz,irm1,MyID,irp1,irstart,irend)
     endif
!-------------------------------------------------------------------------------
   end subroutine divgradz3d_5th
!===============================================================================
!===============================================================================
!  sub. grad6d_5th
   subroutine grad6d_5th(Nx,Ny,Nz,Nux,Nuy,Nuz,A,Bx,By,Bz,Bxx,Byy,Bzz,dx,dy,dz,ixp,iyp,izp)
!-------------------------------------------------------------------------------
!  B = gradient A      using 5th-order central difference
!  B = div ( grad A )  using 5th-order central finite difference
!
!  Input Arrays   : A(Nx,Ny,Nz,Nux,Nuy,Nuz)
!  Input datas    : dx, dy, dz
!
!  Output Arrays  : Bx(Nx,Ny,Nz,Nux,Nuy,Nuz)
!                   By(Nx,Ny,Nz,Nux,Nuy,Nuz)
!                   Bz(Nx,Ny,Nz,Nux,Nuy,Nuz)
!                   Bxx(Nx,Ny,Nz,Nux,Nuy,Nuz)
!                   Byy(Nx,Ny,Nz,Nux,Nuy,Nuz)
!                   Bzz(Nx,Ny,Nz,Nux,Nuy,Nuz)
!
!  ixp, iyp, izp \= 1 : Free Boundary Condition
!  ixp, iyp, izp =  1 : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: A(1), Bx(1), By(1), Bz(1), Bxx(1), Byy(1), Bzz(1)
!     Nxyz = Nx*Ny*Nz
!-------------------------------------------------------------------------------
! 
!  Bx = dA/dx & Bxx = d2A/dx2
!
     if (ixp .eq. 1) then
!      Periodic Boundary Condition along x-direction
       call PCDx5th6D_MPI(A,Bx,Bxx,dx,Nx,Ny,Nz,Nux,Nuy,Nuz,ipm1,MyID,ipp1,ipstart,ipend)
     else
!      Free Boundary Condition along x-direction
       call FCDx5th6D_MPI(A,Bx,Bxx,dx,Nx,Ny,Nz,Nux,Nuy,Nuz,ipm1,MyID,ipp1,ipstart,ipend)
     endif
!-------------------------------------------------------------------------------
! 
!  By = dA/dy & Byy = d2A/dy2
!
     if (Ny == 1 .and. Nz == 1) then
     	 By(1)  = 0.d0
         Byy(1) = 0.d0
     	 Bz(1)  = 0.d0
         Bzz(1) = 0.d0
     	 return
     endif
!
     if (iyp .eq. 1) then 
!      Periodic Boundary Condition along y-direction
       call PCDy5th6D_MPI(A,By,Byy,dy,Nx,Ny,Nz,Nux,Nuy,Nuz,iqm1,MyID,iqp1,iqstart,iqend)
     else 
!      Non-Periodic Boundary Condition along y-direction
       call FCDy5th6D_MPI(A,By,Byy,dy,Nx,Ny,Nz,Nux,Nuy,Nuz,iqm1,MyID,iqp1,iqstart,iqend)
     endif
!-------------------------------------------------------------------------------
! 
!  Bz = dA/dz & Bzz = d2A/dz2
!
     if (Nz == 1) then
     	 Bz(1) = 0.d0
         Bzz(1) = 0.d0
     	 return
     endif
!
     if (izp .eq. 1) then 
!      Periodic Boundary Condition along z-direction
       call PCDz5th6D_MPI(A,Bz,Bzz,dz,Nx,Ny,Nz,Nux,Nuy,Nuz,irm1,MyID,irp1,irstart,irend)
     else
!      Non-Periodic Boundary Condition along z-direction
       call FCDz5th6D_MPI(A,Bz,Bzz,dz,Nx,Ny,Nz,Nux,Nuy,Nuz,irm1,MyID,irp1,irstart,irend)
     endif
!
   end subroutine grad6d_5th
!===============================================================================
!===============================================================================
!===============================================================================
!  sub. grad6d_WENO
   subroutine grad6d_WENO(Nx,Ny,Nz,Nux,Nuy,Nuz,A,Bx,By,Bz,Bxx,Byy,Bzz,dx,dy,dz,ixp,iyp,izp)
!-------------------------------------------------------------------------------
!  B = gradient A      using 5th-order central difference
!  B = div ( grad A )  using 5th-order central finite difference
!
!  Input Arrays   : A(Nx,Ny,Nz,Nux,Nuy,Nuz)
!  Input datas    : dx, dy, dz
!
!  Output Arrays  : Bx(Nx,Ny,Nz,Nux,Nuy,Nuz)
!                   By(Nx,Ny,Nz,Nux,Nuy,Nuz)
!                   Bz(Nx,Ny,Nz,Nux,Nuy,Nuz)
!                   Bxx(Nx,Ny,Nz,Nux,Nuy,Nuz)
!                   Byy(Nx,Ny,Nz,Nux,Nuy,Nuz)
!                   Bzz(Nx,Ny,Nz,Nux,Nuy,Nuz)
!
!  ixp, iyp, izp \= 1 : Free Boundary Condition
!  ixp, iyp, izp =  1 : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: A(1), Bx(1), By(1), Bz(1), Bxx(1), Byy(1), Bzz(1)
!     Nxyz = Nx*Ny*Nz
!-------------------------------------------------------------------------------
! 
!  Bx = dA/dx & Bxx = d2A/dx2
!
     if (ixp .eq. 1) then
!      Periodic Boundary Condition along x-direction
       call WENOx5th6D_MPI(A,Bx,dx,Nx,Ny,Nz,Nux,Nuy,Nuz,ipm1,MyID,ipp1,ipstart,ipend)
     else
!      Free Boundary Condition along x-direction
       call FCDx5th6D_MPI(A,Bx,Bxx,dx,Nx,Ny,Nz,Nux,Nuy,Nuz,ipm1,MyID,ipp1,ipstart,ipend)
     endif
!-------------------------------------------------------------------------------
! 
!  By = dA/dy & Byy = d2A/dy2
!
     if (Ny == 1 .and. Nz == 1) then
       By(1)  = 0.d0
         Byy(1) = 0.d0
       Bz(1)  = 0.d0
         Bzz(1) = 0.d0
       return
     endif
!
     if (iyp .eq. 1) then 
!      Periodic Boundary Condition along y-direction
       call PCDy5th6D_MPI(A,By,Byy,dy,Nx,Ny,Nz,Nux,Nuy,Nuz,iqm1,MyID,iqp1,iqstart,iqend)
     else 
!      Non-Periodic Boundary Condition along y-direction
       call FCDy5th6D_MPI(A,By,Byy,dy,Nx,Ny,Nz,Nux,Nuy,Nuz,iqm1,MyID,iqp1,iqstart,iqend)
     endif
!-------------------------------------------------------------------------------
! 
!  Bz = dA/dz & Bzz = d2A/dz2
!
     if (Nz == 1) then
       Bz(1) = 0.d0
         Bzz(1) = 0.d0
       return
     endif
!
     if (izp .eq. 1) then 
!      Periodic Boundary Condition along z-direction
       call PCDz5th6D_MPI(A,Bz,Bzz,dz,Nx,Ny,Nz,Nux,Nuy,Nuz,irm1,MyID,irp1,irstart,irend)
     else
!      Non-Periodic Boundary Condition along z-direction
       call FCDz5th6D_MPI(A,Bz,Bzz,dz,Nx,Ny,Nz,Nux,Nuy,Nuz,irm1,MyID,irp1,irstart,irend)
     endif
!
   end subroutine grad6d_WENO
!===============================================================================
!===============================================================================
!  sub. grad6d_WENO
   subroutine grad6d_CWENO(Nx,Ny,Nz,Nux,Nuy,Nuz,A,Bx,By,Bz,Bxx,Byy,Bzz,dx,dy,dz,ixp,iyp,izp)
!-------------------------------------------------------------------------------
!  B = gradient A      using 5th-order central difference
!  B = div ( grad A )  using 5th-order central finite difference
!
!  Input Arrays   : A(Nx,Ny,Nz,Nux,Nuy,Nuz)
!  Input datas    : dx, dy, dz
!
!  Output Arrays  : Bx(Nx,Ny,Nz,Nux,Nuy,Nuz)
!                   By(Nx,Ny,Nz,Nux,Nuy,Nuz)
!                   Bz(Nx,Ny,Nz,Nux,Nuy,Nuz)
!                   Bxx(Nx,Ny,Nz,Nux,Nuy,Nuz)
!                   Byy(Nx,Ny,Nz,Nux,Nuy,Nuz)
!                   Bzz(Nx,Ny,Nz,Nux,Nuy,Nuz)
!
!  ixp, iyp, izp \= 1 : Free Boundary Condition
!  ixp, iyp, izp =  1 : Periodic Boundary Condition
!-------------------------------------------------------------------------------
     implicit double precision (a-h,o-z)
     dimension :: A(1), Bx(1), By(1), Bz(1), Bxx(1), Byy(1), Bzz(1)
!     Nxyz = Nx*Ny*Nz
!-------------------------------------------------------------------------------
! 
!  Bx = dA/dx & Bxx = d2A/dx2
!
     if (ixp .eq. 1) then
!      Periodic Boundary Condition along x-direction
       call CWENOx6th6D_MPI(A,Bx,dx,Nx,Ny,Nz,Nux,Nuy,Nuz,ipm1,MyID,ipp1,ipstart,ipend)
     else
!      Free Boundary Condition along x-direction
       call FCDx5th6D_MPI(A,Bx,Bxx,dx,Nx,Ny,Nz,Nux,Nuy,Nuz,ipm1,MyID,ipp1,ipstart,ipend)
     endif
!-------------------------------------------------------------------------------
! 
!  By = dA/dy & Byy = d2A/dy2
!
     if (Ny == 1 .and. Nz == 1) then
       By(1)  = 0.d0
         Byy(1) = 0.d0
       Bz(1)  = 0.d0
         Bzz(1) = 0.d0
       return
     endif
!
     if (iyp .eq. 1) then 
!      Periodic Boundary Condition along y-direction
       call PCDy5th6D_MPI(A,By,Byy,dy,Nx,Ny,Nz,Nux,Nuy,Nuz,iqm1,MyID,iqp1,iqstart,iqend)
     else 
!      Non-Periodic Boundary Condition along y-direction
       call FCDy5th6D_MPI(A,By,Byy,dy,Nx,Ny,Nz,Nux,Nuy,Nuz,iqm1,MyID,iqp1,iqstart,iqend)
     endif
!-------------------------------------------------------------------------------
! 
!  Bz = dA/dz & Bzz = d2A/dz2
!
     if (Nz == 1) then
       Bz(1) = 0.d0
         Bzz(1) = 0.d0
       return
     endif
!
     if (izp .eq. 1) then 
!      Periodic Boundary Condition along z-direction
       call PCDz5th6D_MPI(A,Bz,Bzz,dz,Nx,Ny,Nz,Nux,Nuy,Nuz,irm1,MyID,irp1,irstart,irend)
     else
!      Non-Periodic Boundary Condition along z-direction
       call FCDz5th6D_MPI(A,Bz,Bzz,dz,Nx,Ny,Nz,Nux,Nuy,Nuz,irm1,MyID,irp1,irstart,irend)
     endif
!
   end subroutine grad6d_CWENO
!===============================================================================
!===============================================================================
end module vector3D
