module CenDif3D_MPI
    use MPI3D
    use coef
    use omp_lib
  contains
!===============================================================================
!  calculate the grid size along the x-direction
!
  subroutine geth_MPI(x,h,N,iDm1,iD,iDp1,iDstart,iDend)
    implicit double precision (A-H,O-Z)
    dimension :: x(1), h(1), xp(1)
    dimension :: hmax_local(1), hmin_local(1), hmax(1), hmin(1)
    character*10 :: dates, times, zone
    integer, dimension(:) :: values(8)
!-------------------------------------------------------------------------------
!  send and receive the boundary
!
    call MPI_SENDRECV(x(1) ,1,MPI_REAL8,iDm1,100, &
                      xp(1),1,MPI_REAL8,iDp1,100, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
    do i = 1, N-1
      h(i) = DIFAB(x(i+1),x(i))
!     h(i) = x(i+1)-x(i)
    enddo
!
    call MPI_BARRIER(MPI_WORLD,IERR)
!-------------------------------------------------------------------------------
!  send the boundary from ip = 1 and receive the boundary at ip = P
!    for periodic case
!
    if (iDend == iDstart) then
      h(N) = h(1)
      goto 123
    endif
!
    if (iD == iDend) then
      call MPI_RECV(h(N),1,MPI_REAL8,iDstart,110,MPI_WORLD,ISTATUS,IERR)
    else 
      h(N) = DIFAB(xp(1),x(N))
!     h(N) = xp(1)-x(N)
      if (iD==iDstart) call MPI_SEND(h(1),1,MPI_REAL8,iDend,110,MPI_WORLD,IERR)
    endif
!-------------------------------------------------------------------------------
!  find the max and min vlaue in the array h
!
123 continue
    htemp = 0.d0
    do i = 1, N
      if (h(i) >= htemp) htemp = h(i)
    enddo
    hmax_local(1) = htemp
!
    do i = 1, N
      if (h(i) <= htemp) htemp = h(i)
    enddo
    hmin_local(1) = htemp
!
    call MPI_BARRIER(MPI_WORLD,IERR) 
!
    call MPI_ALLREDUCE(hmax_local(1),hmax(1),1,MPI_REAL8,MPI_MAX,MPI_WORLD,IERR)
    call MPI_ALLREDUCE(hmin_local(1),hmin(1),1,MPI_REAL8,MPI_MIN,MPI_WORLD,IERR)
!-------------------------------------------------------------------------------
!  check the hmin
!
    if (hmin(1)/hmax(1) <= relerr) then
      if (MyID == 0) then
        call Date_and_Time(dates,times,zone,values)
        open (1,file='error.txt',status='unknown')
        write (1,11) values(1), '-', values(2), '-', values(3), ' ', &
                     values(5), ':', values(6)
        write (1,*) 'hmax =', hmax(1)
        write (1,*) 'hmin =', hmin(1)
        write (1,*) 'hmin is too small !!'
        write (1,*) 'STOP program at geth_MPI !!'
        close (1)
      endif      
!      
      call MPI3D_FINALIZE
      stop
    endif
!
 11 format (1x,i4,4(a1,i2.2))
  end subroutine geth_MPI
!===============================================================================
!===============================================================================
!  This is the first step of first-order central difference
!    to prepare the table A, B, C.
!
  subroutine PCD1st0_MPI(h,A,B,C,N,iDm1,iD,iDp1,iDstart,iDend)
    implicit double precision (A-H,O-Z)
    dimension :: h(1), A(1), B(1), C(1)
    dimension :: hm(1)
!-------------------------------------------------------------------------------
!  send and receive the boundary
!
    call MPI_SENDRECV(h(N) ,1,MPI_REAL8,iDp1,100, &
                      hm(1),1,MPI_REAL8,iDm1,100, &
                      MPI_WORLD,ISTATUS,IERR)
!
    if (iDend == iDstart) then
      hm(1) = h(N-1)
      goto 123
    endif
!
    if (iD == iDstart) call MPI_RECV(hm(1),1,MPI_REAL8,iDend,110, &
                                     MPI_WORLD,ISTATUS,IERR)
    if (iD == iDend)   call MPI_SEND(h(N-1),1,MPI_REAL8,iDstart,110, &
                                     MPI_WORLD,IERR)
123 continue
!-------------------------------------------------------------------------------
    do i = 2, N
      h1 = 1.d0/h(i)
      h2 = 1.d0/h(i-1)
!     hp = DIFAB(h(i),-h(i-1))
!     h3 = 1.d0/hp
      h3 = 1.d0/(h(i)+h(i-1))
      A(i) = h(i-1)*h1*h3
      B(i) = h2-h1
      C(i) = -h(i)*h2*h3
    enddo
!-------------------------------------------------------------------------------
!  for i = 1
    h1 = 1.d0/h(1)
    h2 = 1.d0/hm(1)
!   hp = DIFAB(h(1),-hm(1))
!   h3 = 1.d0/hp
    h3 = 1.d0/(h(1)+hm(1))
    A(1) = hm(1)*h1*h3
    B(1) = h2-h1
    C(1) = -h(1)*h2*h3
!
    call MPI_BARRIER(MPI_WORLD,IERR)
!-------------------------------------------------------------------------------        
  end subroutine PCD1st0_MPI
!===============================================================================
!===============================================================================
!  This is the second step of first-order central difference
!    to calculate the first derivative fp along x-direction
!    for periodic boundary condition.
!
  subroutine PCDx1st3D_MPI(f,A,B,C,fp,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
    implicit double precision (A-H,O-Z)
    dimension :: f(1), fp(1)
    dimension :: A(1), B(1), C(1)
    real*8, allocatable, dimension(:) :: fm1, fp1, fsend
    Nxy = Nx*Ny
    Nyz = Ny*Nz
    allocate (fm1(Nyz), fp1(Nyz), fsend(Nyz))
!-------------------------------------------------------------------------------
!  if P = 1, it doesn't transfer data in MPI
!
    if (iDend == iDstart) then
      call getBCx3D(f(1),fp1(1),2,Nx,Ny,Nz)
      call getBCx3D(f(1),fm1(1),Nx-1,Nx,Ny,Nz)
      goto 123
    endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
    call getBCx3D(f(1),fsend(1),1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nyz,MPI_REAL8,iDm1,100, &
                      fp1(1)  ,Nyz,MPI_REAL8,iDp1,100, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
    call getBCx3D(f(1),fsend(1),Nx,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nyz,MPI_REAL8,iDp1,110, &
                      fm1(1)  ,Nyz,MPI_REAL8,iDm1,110, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
    if (iD == iDstart) then
      call getBCx3D(f(1),fsend(1),2,Nx,Ny,Nz)
      call MPI_SEND(fsend(1),Nyz,MPI_REAL8,iDend,120,MPI_WORLD,IERR)
      call MPI_RECV(fm1(1),Nyz,MPI_REAL8,iDend,130,MPI_WORLD,ISTATUS,IERR)      
    endif
!
    if (iD == iDend) then
      call getBCx3D(f(1),fsend(1),Nx-1,Nx,Ny,Nz)      
      call MPI_RECV(fp1(1),Nyz,MPI_REAL8,iDstart,120,MPI_WORLD,ISTATUS,IERR)
      call MPI_SEND(fsend(1),Nyz,MPI_REAL8,iDstart,130,MPI_WORLD,IERR)      
    endif
!    
123 continue                                          
!-------------------------------------------------------------------------------
    ii = 0
    do 30 k = 1, Nz
    do 30 j = 1, Ny
      ii = ii+1
      jk = (j-1)*Nx+(k-1)*Nxy
!      
      do i = 2, Nx-1
        ijk   = i  +jk
        ijkm1 = i-1+jk
        ijkp1 = i+1+jk
        temp1 = A(i)*f(ijkp1)
        temp2 = B(i)*f(ijk)
        temp3 = C(i)*f(ijkm1)
        sum   = DIFAB(temp1,-temp2)
        fp(ijk) = DIFAB(sum,-temp3)
!       fp(ijk) = A(i)*f(ijkp1)+B(i)*f(ijk)+C(i)*f(ijkm1)
      enddo
!
      ijk   = 1+jk
      ijkp1 = 2+jk
      temp1 = A(1)*f(ijkp1)
      temp2 = B(1)*f(ijk)
      temp3 = C(1)*fm1(ii)
      sum = DIFAB(temp1,-temp2)
      fp(ijk) = DIFAB(sum,-temp3)
!     fp(ijk) = A(1)*f(ijkp1)+B(1)*f(ijk)+C(1)*fm1(ii)
!
      ijk   = Nx  +jk
      ijkm1 = Nx-1+jk
      temp1 = A(Nx)*fp1(ii)
      temp2 = B(Nx)*f(ijk)
      temp3 = C(Nx)*f(ijkm1)
      sum = DIFAB(temp1,-temp2)
      fp(ijk) = DIFAB(sum,-temp3)
!     fp(ijk) = A(Nx)*fp1(ii)+B(Nx)*f(ijk)+C(Nx)*f(ijkm1)
 30 continue
!-------------------------------------------------------------------------------
    deallocate (fm1, fp1, fsend)  
  end subroutine PCDx1st3D_MPI
!===============================================================================
!  return the Boundary values along the x-direction
!
  subroutine getBCx3D(f,fBC,ix,Nx,Ny,Nz)
    implicit double precision (A-H,O-Z)
    dimension :: f(1), fBC(1)
    Nxy = Nx*Ny
    ii = 0
    do 10 k = 1, Nz
    do 10 j = 1, Ny
      ii = ii+1
      ijk = ix+(j-1)*Nx+(k-1)*Nxy
      fBC(ii) = f(ijk)        
 10 continue
  end subroutine getBCx3D
!===============================================================================
!===============================================================================
!  This is the second step of first-order central difference
!    to calculate the first derivative fp along y-direction
!    for periodic boundary condition.
!
  subroutine PCDy1st3D_MPI(f,A,B,C,fp,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
    implicit double precision (A-H,O-Z)
    dimension :: f(1), fp(1)
    dimension :: A(1), B(1), C(1)
    real*8, allocatable, dimension(:) :: fm1, fp1, fsend
    Nxy = Nx*Ny
    Nxz = Nx*Nz
    allocate (fm1(Nxz), fp1(Nxz), fsend(Nxz))
!-------------------------------------------------------------------------------
!  if Q = 1, it doesn't transfer data in MPI
!
    if (iDend == iDstart) then
      call getBCy3D(f(1),fp1(1),2,Nx,Ny,Nz)
      call getBCy3D(f(1),fm1(1),Ny-1,Nx,Ny,Nz)
      goto 123
    endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
    call getBCy3D(f(1),fsend(1),1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxz,MPI_REAL8,iDm1,100, &
                      fp1(1)  ,Nxz,MPI_REAL8,iDp1,100, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
    call getBCy3D(f(1),fsend(1),Ny,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxz,MPI_REAL8,iDp1,110, &
                      fm1(1)  ,Nxz,MPI_REAL8,iDm1,110, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
    if (iD == iDstart) then
      call getBCy3D(f(1),fsend(1),2,Nx,Ny,Nz)
      call MPI_SEND(fsend(1),Nxz,MPI_REAL8,iDend,120,MPI_WORLD,IERR)
      call MPI_RECV(fm1(1),Nxz,MPI_REAL8,iDend,130,MPI_WORLD,ISTATUS,IERR)      
    endif
!
    if (iD == iDend) then
      call getBCy3D(f(1),fsend(1),Ny-1,Nx,Ny,Nz)      
      call MPI_RECV(fp1(1),Nxz,MPI_REAL8,iDstart,120,MPI_WORLD,ISTATUS,IERR)
      call MPI_SEND(fsend(1),Nxz,MPI_REAL8,iDstart,130,MPI_WORLD,IERR)      
    endif
!    
123 continue                                          
!-------------------------------------------------------------------------------
!    goto 124
    ii = 0
    do 30 k = 1, Nz
    do 30 i = 1, Nx
      ii = ii+1
      ik = i+(k-1)*Nxy
!      
      do j = 2, Ny-1
        ijk   = ( j   -1)*Nx+ik
        ijkm1 = ((j-1)-1)*Nx+ik
        ijkp1 = ((j+1)-1)*Nx+ik
        temp1 = A(j)*f(ijkp1)
        temp2 = B(j)*f(ijk)
        temp3 = C(j)*f(ijkm1)
        sum   = DIFAB(temp1,-temp2)
        fp(ijk) = DIFAB(sum,-temp3)
!       fp(ijk) = A(j)*f(ijkp1)+B(j)*f(ijk)+C(j)*f(ijkm1)
      enddo
!
      ijk   = (1-1)*Nx+ik
      ijkp1 = (2-1)*Nx+ik
      temp1 = A(1)*f(ijkp1)
      temp2 = B(1)*f(ijk)
      temp3 = C(1)*fm1(ii)
      sum   = DIFAB(temp1,-temp2)
      fp(ijk) = DIFAB(sum,-temp3)
!     fp(ijk) = A(1)*f(ijkp1)+B(1)*f(ijk)+C(1)*fm1(ii)
!
      ijk   = ( Ny   -1)*Nx+ik
      ijkm1 = ((Ny-1)-1)*Nx+ik
      temp1 = A(Ny)*fp1(ii)
      temp2 = B(Ny)*f(ijk)
      temp3 = C(Ny)*f(ijkm1)
      sum   = DIFAB(temp1,-temp2)
      fp(ijk) = DIFAB(sum,-temp3)
!     fp(ijk) = A(Ny)*fp1(ii)+B(Ny)*f(ijk)+C(Ny)*f(ijkm1)
 30 continue
!-------------------------------------------------------------------------------
124 continue
    deallocate (fm1, fp1, fsend)  
  end subroutine PCDy1st3D_MPI
!===============================================================================
!  return the Boundary values along the y-direction
!
  subroutine getBCy3D(f,fBC,jy,Nx,Ny,Nz)
    implicit double precision (A-H,O-Z)
    dimension :: f(1), fBC(1)
    Nxy = Nx*Ny
    ii = 0
    do 10 k = 1, Nz
    do 10 i = 1, Nx
      ii = ii+1
      ijk = i+(jy-1)*Nx+(k-1)*Nxy
      fBC(ii) = f(ijk)        
 10 continue
  end subroutine getBCy3D
!===============================================================================
!===============================================================================
!  This is the second step of first-order central difference
!    to calculate the first derivative fp along z-direction
!    for periodic boundary condition.
!
  subroutine PCDz1st3D_MPI(f,A,B,C,fp,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
    implicit double precision (A-H,O-Z)
    dimension :: f(1), fp(1)
    dimension :: A(1), B(1), C(1)
    real*8, allocatable, dimension(:) :: fm1, fp1, fsend
    Nxy = Nx*Ny
    allocate (fm1(Nxy), fp1(Nxy), fsend(Nxy))
!-------------------------------------------------------------------------------
!  if P = 1, it doesn't transfer data in MPI
!
    if (iDend == iDstart) then
      call getBCz3D(f(1),fp1(1),2,Nx,Ny,Nz)
      call getBCz3D(f(1),fm1(1),Nz-1,Nx,Ny,Nz)
      goto 123
    endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
    call getBCz3D(f(1),fsend(1),1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxy,MPI_REAL8,iDm1,100, &
                      fp1(1)  ,Nxy,MPI_REAL8,iDp1,100, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
    call getBCz3D(f(1),fsend(1),Nz,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxy,MPI_REAL8,iDp1,110, &
                      fm1(1)  ,Nxy,MPI_REAL8,iDm1,110, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
    if (iD == iDstart) then
      call getBCz3D(f(1),fsend(1),2,Nx,Ny,Nz)
      call MPI_SEND(fsend(1),Nxy,MPI_REAL8,iDend,120,MPI_WORLD,IERR)
      call MPI_RECV(fm1(1),Nxy,MPI_REAL8,iDend,130,MPI_WORLD,ISTATUS,IERR)      
    endif
!
    if (iD == iDend) then
      call getBCz3D(f(1),fsend(1),Nz-1,Nx,Ny,Nz)      
      call MPI_RECV(fp1(1),Nxy,MPI_REAL8,iDstart,120,MPI_WORLD,ISTATUS,IERR)
      call MPI_SEND(fsend(1),Nxy,MPI_REAL8,iDstart,130,MPI_WORLD,IERR)      
    endif
!    
123 continue                                          
!-------------------------------------------------------------------------------
    ii = 0
    do 30 j = 1, Ny
    do 30 i = 1, Nx
      ii = ii+1
      ij = i+(j-1)*Nx
!      
      do k = 2, Nz-1
        ijk   = ij+( k   -1)*Nxy
        ijkm1 = ij+((k-1)-1)*Nxy
        ijkp1 = ij+((k+1)-1)*Nxy
        temp1 = A(k)*f(ijkp1)
        temp2 = B(k)*f(ijk)
        temp3 = C(k)*f(ijkm1)
        sum = DIFAB(temp1,-temp2)
        fp(ijk) = DIFAB(sum,-temp3)
!       fp(ijk) = A(k)*f(ijkp1)+B(k)*f(ijk)+C(k)*f(ijkm1)
      enddo
!
      ijk   = ij+(1-1)*Nxy
      ijkp1 = ij+(2-1)*Nxy
      temp1 = A(1)*f(ijkp1)
      temp2 = B(1)*f(ijk)
      temp3 = C(1)*fm1(ii)
      sum = DIFAB(temp1,-temp2)
      fp(ijk) = DIFAB(sum,-temp3)
!     fp(ijk) = A(1)*f(ijkp1)+B(1)*f(ijk)+C(1)*fm1(ii)
!
      ijk   = ij+( Nz   -1)*Nxy
      ijkm1 = ij+((Nz-1)-1)*Nxy
      temp1 = A(Nz)*fp1(ii)
      temp2 = B(Nz)*f(ijk)
      temp3 = C(Nz)*f(ijkm1)
      sum = DIFAB(temp1,-temp2)
      fp(ijk) = DIFAB(sum,-temp3)
!     fp(ijk) = A(Nz)*fp1(ii)+B(Nz)*f(ijk)+C(Nz)*f(ijkm1)
 30 continue
!-------------------------------------------------------------------------------
    deallocate (fm1, fp1, fsend)  
  end subroutine PCDz1st3D_MPI
!===============================================================================
!  return the Boundary values along the z-direction
!
  subroutine getBCz3D(f,fBC,kz,Nx,Ny,Nz)
    implicit double precision (A-H,O-Z)
    dimension :: f(1), fBC(1)
    Nxy = Nx*Ny
    ii = 0
    do 10 j = 1, Ny
    do 10 i = 1, Nx
      ii = ii+1
      ijk = i+(j-1)*Nx+(kz-1)*Nxy
      fBC(ii) = f(ijk)        
 10 continue
  end subroutine getBCz3D
!===============================================================================
!===============================================================================
!  This is the second step of first-order central difference
!    to calculate the first derivative fp along x-direction
!    for free boundary condition.
!
  subroutine FCDx1st3D_MPI(f,A,B,C,fp,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
    implicit double precision (A-H,O-Z)
    dimension :: f(1), fp(1)
    dimension :: A(1), B(1), C(1)
    real*8, allocatable, dimension(:) :: fm1, fp1, fsend
    Nxy = Nx*Ny
    Nyz = Ny*Nz
    allocate (fm1(Nyz), fp1(Nyz), fsend(Nyz))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
    call getBCx3D(f(1),fsend(1),1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nyz,MPI_REAL8,iDm1,100, &
                      fp1(1)  ,Nyz,MPI_REAL8,iDp1,100, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
    call getBCx3D(f(1),fsend(1),Nx,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nyz,MPI_REAL8,iDp1,110, &
                      fm1(1)  ,Nyz,MPI_REAL8,iDm1,110, &
                      MPI_WORLD,ISTATUS,IERR)                                       
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
    if (iD == iDstart) then
      call getBCx3D(f(1),fm1(1),1,Nx,Ny,Nz)      
    endif
!
    if (iD == iDend) then           
      call getBCx3D(f(1),fp1(1),Nx,Nx,Ny,Nz)      
    endif  
!-------------------------------------------------------------------------------
    ii = 0
    do 30 k = 1, Nz
    do 30 j = 1, Ny
      ii = ii+1
      jk = (j-1)*Nx+(k-1)*Nxy
!      
      do i = 2, Nx-1
        ijk   = i  +jk
        ijkm1 = i-1+jk
        ijkp1 = i+1+jk
        temp1 = A(i)*f(ijkp1)
        temp2 = B(i)*f(ijk)
        temp3 = C(i)*f(ijkm1)
        sum   = DIFAB(temp1,-temp2)
        fp(ijk) = DIFAB(sum,-temp3)
!       fp(ijk) = A(i)*f(ijkp1)+B(i)*f(ijk)+C(i)*f(ijkm1)
      enddo
!
      ijk   = 1+jk
      ijkp1 = 2+jk
!      if (iD == iDstart) then
!        fp(ijk) = 0.d0
!      else
        temp1 = A(1)*f(ijkp1)
        temp2 = B(1)*f(ijk)
        temp3 = C(1)*fm1(ii)
        sum   = DIFAB(temp1,-temp2)
        fp(ijk) = DIFAB(sum,-temp3)
!       fp(ijk) = A(1)*f(ijkp1)+B(1)*f(ijk)+C(1)*fm1(ii)
!      endif
!
      ijk   = Nx  +jk
      ijkm1 = Nx-1+jk
!      if (iD == iDend) then
!        fp(ijk) = 0.d0
!      else
        temp1 = A(Nx)*fp1(ii)
        temp2 = B(Nx)*f(ijk)
        temp3 = C(Nx)*f(ijkm1)
        sum   = DIFAB(temp1,-temp2)
        fp(ijk) = DIFAB(sum,-temp3)
!       fp(ijk) = A(Nx)*fp1(ii)+B(Nx)*f(ijk)+C(Nx)*f(ijkm1)
!      endif
 30 continue
!-------------------------------------------------------------------------------
    deallocate (fm1, fp1, fsend)  
  end subroutine FCDx1st3D_MPI
!===============================================================================!===============================================================================
!  This is the second step of first-order central difference
!    to calculate the first derivative fp along y-direction
!    for free boundary condition.
!
  subroutine FCDy1st3D_MPI(f,A,B,C,fp,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
    implicit double precision (A-H,O-Z)
    dimension :: f(1), fp(1)
    dimension :: A(1), B(1), C(1)
    real*8, allocatable, dimension(:) :: fm1, fp1, fsend
    Nxy = Nx*Ny
    Nxz = Nx*Nz
    allocate (fm1(Nxz), fp1(Nxz), fsend(Nxz))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
    call getBCy3D(f(1),fsend(1),1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxz,MPI_REAL8,iDm1,100, &
                      fp1(1)  ,Nxz,MPI_REAL8,iDp1,100, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
    call getBCy3D(f(1),fsend(1),Ny,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxz,MPI_REAL8,iDp1,110, &
                      fm1(1)  ,Nxz,MPI_REAL8,iDm1,110, &
                      MPI_WORLD,ISTATUS,IERR)                                          
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
    if (iD == iDstart) then
      call getBCy3D(f(1),fm1(1),1,Nx,Ny,Nz)       
    endif
!
    if (iD == iDend) then           
      call getBCy3D(f(1),fp1(1),Ny,Nx,Ny,Nz)      
    endif  
!-------------------------------------------------------------------------------
!    goto 124
    ii = 0
    do 30 k = 1, Nz
    do 30 i = 1, Nx
      ii = ii+1
      ik = i+(k-1)*Nxy
!      
      do j = 2, Ny-1
        ijk   = ( j   -1)*Nx+ik
        ijkm1 = ((j-1)-1)*Nx+ik
        ijkp1 = ((j+1)-1)*Nx+ik
        temp1 = A(j)*f(ijkp1)
        temp2 = B(j)*f(ijk)
        temp3 = C(j)*f(ijkm1)
        sum   = DIFAB(temp1,-temp2)
        fp(ijk) = DIFAB(sum,-temp3)
!       fp(ijk) = A(j)*f(ijkp1)+B(j)*f(ijk)+C(j)*f(ijkm1)
      enddo
!
      ijk   = (1-1)*Nx+ik
      ijkp1 = (2-1)*Nx+ik
!      if (iD == iDstart) then
!        fp(ijk) = 0.d0
!      else
        temp1 = A(1)*f(ijkp1)
        temp2 = B(1)*f(ijk)
        temp3 = C(1)*fm1(ii)
        sum   = DIFAB(temp1,-temp2)
        fp(ijk) = DIFAB(sum,-temp3)
!       fp(ijk) = A(1)*f(ijkp1)+B(1)*f(ijk)+C(1)*fm1(ii)
!      endif
!
      ijk   = ( Ny   -1)*Nx+ik
      ijkm1 = ((Ny-1)-1)*Nx+ik
!      if (iD == iDend) then
!        fp(ijk) = 0.d0
!      else
        temp1 = A(Ny)*fp1(ii)
        temp2 = B(Ny)*f(ijk)
        temp3 = C(Ny)*f(ijkm1)
        sum   = DIFAB(temp1,-temp2)
        fp(ijk) = DIFAB(sum,-temp3)
!       fp(ijk) = A(Ny)*fp1(ii)+B(Ny)*f(ijk)+C(Ny)*f(ijkm1)
!      endif
 30 continue
!-------------------------------------------------------------------------------
124 continue
    deallocate (fm1, fp1, fsend)  
  end subroutine FCDy1st3D_MPI
!===============================================================================
!===============================================================================
!  This is the second step of first-order central difference
!    to calculate the first derivative fp along z-direction
!    for free boundary condition.
!
  subroutine FCDz1st3D_MPI(f,A,B,C,fp,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
    implicit double precision (A-H,O-Z)
    dimension :: f(1), fp(1)
    dimension :: A(1), B(1), C(1)
    real*8, allocatable, dimension(:) :: fm1, fp1, fsend
    Nxy = Nx*Ny
    allocate (fm1(Nxy), fp1(Nxy), fsend(Nxy))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
    call getBCz3D(f(1),fsend(1),1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxy,MPI_REAL8,iDm1,100, &
                      fp1(1)  ,Nxy,MPI_REAL8,iDp1,100, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
    call getBCz3D(f(1),fsend(1),Nz,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxy,MPI_REAL8,iDp1,110, &
                      fm1(1)  ,Nxy,MPI_REAL8,iDm1,110, &
                      MPI_WORLD,ISTATUS,IERR)                                          
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
    if (iD == iDstart) then
      call getBCz3D(f(1),fm1(1),1,Nx,Ny,Nz)      
    endif
!
    if (iD == iDend) then           
      call getBCz3D(f(1),fp1(1),Nz,Nx,Ny,Nz)      
    endif  
!-------------------------------------------------------------------------------
    ii = 0
    do 30 j = 1, Ny
    do 30 i = 1, Nx
      ii = ii+1
      ij = i+(j-1)*Nx
!      
      do k = 2, Nz-1
        ijk   = ij+( k   -1)*Nxy
        ijkm1 = ij+((k-1)-1)*Nxy
        ijkp1 = ij+((k+1)-1)*Nxy
        temp1 = A(k)*f(ijkp1)
        temp2 = B(k)*f(ijk)
        temp3 = C(k)*f(ijkm1)
        sum   = DIFAB(temp1,-temp2)
        fp(ijk) = DIFAB(sum,-temp3)
!       fp(ijk) = A(k)*f(ijkp1)+B(k)*f(ijk)+C(k)*f(ijkm1)
      enddo
!
      ijk   = ij+(1-1)*Nxy
      ijkp1 = ij+(2-1)*Nxy
!      if (iD == iDstart) then
!        fp(ijk) = 0.d0
!      else
        temp1 = A(1)*f(ijkp1)
        temp2 = B(1)*f(ijk)
        temp3 = C(1)*fm1(ii)
        sum   = DIFAB(temp1,-temp2)
        fp(ijk) = DIFAB(sum,-temp3)
!       fp(ijk) = A(1)*f(ijkp1)+B(1)*f(ijk)+C(1)*fm1(ii)
!      endif
!
      ijk   = ij+( Nz   -1)*Nxy
      ijkm1 = ij+((Nz-1)-1)*Nxy
!      if (iD == iDend) then
!        fp(ijk) = 0.d0
!      else
        temp1 = A(Nz)*fp1(ii)
        temp2 = B(Nz)*f(ijk)
        temp3 = C(Nz)*f(ijkm1)
        sum   = DIFAB(temp1,-temp2)
        fp(ijk) = DIFAB(sum,-temp3)
!       fp(ijk) = A(Nz)*fp1(ii)+B(Nz)*f(ijk)+C(Nz)*f(ijkm1)
!      endif
 30 continue
!-------------------------------------------------------------------------------
    deallocate (fm1, fp1, fsend)  
  end subroutine FCDz1st3D_MPI
!===============================================================================
!===============================================================================
!  This is the first step of 3rd-order central difference in x-direction
!    to prepare the table A, B, C, D, E
!
  subroutine PCD3rd0_MPI(h,A,B,C,D,E,N,iDm1,iD,iDp1,iDstart,iDend)
    implicit double precision (A-H,O-Z)
    dimension :: h(1), A(1), B(1), C(1), D(1), E(1)  
    real*8, allocatable, dimension(:) :: hm, hp
    allocate (hm(2), hp(1))
    i2 = 2
    i3 = 3
!-------------------------------------------------------------------------------
!  if P = 1, it doesn't transfer data in MPI
!
    if (iDend == iDstart) then
      hp(1) = h(i2)
      hm(1) = h(N-2)
      hm(2) = h(N-1)      
      goto 123
    endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
    call MPI_SENDRECV(h(1) ,1,MPI_REAL8,iDm1,100, &
                      hp(1),1,MPI_REAL8,iDp1,100, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
    call MPI_SENDRECV(h(N-1),2,MPI_REAL8,iDp1,110, &
                      hm(1) ,2,MPI_REAL8,iDm1,110, &
                      MPI_WORLD,ISTATUS,IERR)                                          
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
    if (iD == iDstart) then
      call MPI_SEND(h(i2),1,MPI_REAL8,iDend,120,MPI_WORLD,IERR)
      call MPI_RECV(hm(1),2,MPI_REAL8,iDend,130,MPI_WORLD,ISTATUS,IERR)      
    endif
!
    if (iD == iDend) then     
      call MPI_RECV(hp(1),1,MPI_REAL8,iDstart,120,MPI_WORLD,ISTATUS,IERR)
      call MPI_SEND(h(N-2),2,MPI_REAL8,iDstart,130,MPI_WORLD,IERR)      
    endif
!    
123 continue
!-------------------------------------------------------------------------------
    do i = 3, N-1
      ht1 = h(i+1)+h(i)
      ht2 = h(i)+h(i-1)
      ht3 = h(i-1)+h(i-2)
      ht4 = h(i+1)+h(i)+h(i-1)
      ht5 = h(i)+h(i-1)+h(i-2)
      ht6 = h(i+1)+h(i)+h(i-1)+h(i-2)
      A(i) = -h(i)*h(i-1)*ht3/(h(i+1)*ht1*ht4*ht6)
      B(i) =  h(i-1)*ht1*ht3/(h(i)*h(i+1)*ht2*ht5)
!
      temp1 = 1.d0/ht1
      temp2 = 1.d0/h(i)
      temp3 = 1.d0/h(i-1)
      temp4 = 1.d0/ht3
      sum   = DIFAB(temp3,temp1)
      sum   = DIFAB(sum,temp2)
      C(i)  = DIFAB(sum,-temp4)
!     C(i)  = -1.d0/ht1-1.d0/h(i)+1.d0/h(i-1)+1.d0/ht3
!
      D(i) = -h(i)*ht1*ht3/(h(i-1)*h(i-2)*ht2*ht4)
      E(i) =  h(i)*h(i-1)*ht1/(h(i-2)*ht3*ht5*ht6)
    enddo
!-------------------------------------------------------------------------------
!   for i = 1
    ht1 = h(i2)+h(1)
    ht2 = h(1)+hm(2)
    ht3 = hm(2)+hm(1)
    ht4 = h(i2)+h(1)+hm(2)
    ht5 = h(1)+hm(2)+hm(1)
    ht6 = h(i2)+h(1)+hm(2)+hm(1)
    A(1) = -h(1)*hm(2)*ht3/(h(i2)*ht1*ht4*ht6)
    B(1) =  hm(2)*ht1*ht3/(h(1)*h(i2)*ht2*ht5)
!
    temp1 = 1.d0/ht1
    temp2 = 1.d0/h(1)
    temp3 = 1.d0/hm(2)
    temp4 = 1.d0/ht3
    sum   = DIFAB(temp3,temp1)
    sum   = DIFAB(sum,temp2)
    C(1)  = DIFAB(sum,-temp4)
!   C(1)  = -1.d0/ht1-1.d0/h(1)+1.d0/hm(2)+1.d0/ht3
!
    D(1) = -h(1)*ht1*ht3/(hm(2)*hm(1)*ht2*ht4)
    E(1) =  h(1)*hm(2)*ht1/(hm(1)*ht3*ht5*ht6)
!-------------------------------------------------------------------------------
!   for i = 2
    ht1 = h(i3)+h(i2)
    ht2 = h(i2)+h(1)      
    ht3 = h(1)+hm(2)
    ht4 = h(i3)+h(i2)+h(1)      
    ht5 = h(i2)+h(1)+hm(2)      
    ht6 = h(i3)+h(i2)+h(1)+hm(2)
    A(i2) = -h(i2)*h(1)*ht3/(h(i3)*ht1*ht4*ht6)
    B(i2) =  h(1)*ht1*ht3/(h(i2)*h(i3)*ht2*ht5)
!
    temp1 = 1.d0/ht1
    temp2 = 1.d0/h(i2)
    temp3 = 1.d0/h(1)
    temp4 = 1.d0/ht3
    sum   = DIFAB(temp3,temp1)
    sum   = DIFAB(sum,temp2)
    C(i2) = DIFAB(sum,-temp4)
!   C(i2) = -1.d0/ht1-1.d0/h(i2)+1.d0/h(1)+1.d0/ht3
!
    D(i2) = -h(i2)*ht1*ht3/(h(1)*hm(2)*ht2*ht4)      
    E(i2) =  h(i2)*h(1)*ht1/(hm(2)*ht3*ht5*ht6)
!-------------------------------------------------------------------------------
!   for i = N
    ht1 = hp(1)+h(N)
    ht2 = h(N)+h(N-1)
    ht3 = h(N-1)+h(N-2)
    ht4 = hp(1)+h(N)+h(N-1)
    ht5 = h(N)+h(N-1)+h(N-2)
    ht6 = hp(1)+h(N)+h(N-1)+h(N-2)
    A(N) = -h(N)*h(N-1)*ht3/(hp(1)*ht1*ht4*ht6)
    B(N) =  h(N-1)*ht1*ht3/(h(N)*hp(1)*ht2*ht5)
!
    temp1 = 1.d0/ht1
    temp2 = 1.d0/h(N)
    temp3 = 1.d0/h(N-1)
    temp4 = 1.d0/ht3
    sum   = DIFAB(temp3,temp1)
    sum   = DIFAB(sum,temp2)
    C(N)  = DIFAB(sum,-temp4)
!   C(N)  = -1.d0/ht1-1.d0/h(N)+1.d0/h(N-1)+1.d0/ht3
!
    D(N) = -h(N)*ht1*ht3/(h(N-1)*h(N-2)*ht2*ht4)
    E(N) =  h(N)*h(N-1)*ht1/(h(N-2)*ht3*ht5*ht6)
!-------------------------------------------------------------------------------
    deallocate (hm, hp)
  end subroutine PCD3rd0_MPI  
!===============================================================================
!===============================================================================
!  This is the second step of third-order central difference
!    to calculate the first derivative fp along x-direction
!    for periodic boundary condition.
!
  subroutine PCDx3rd3D_MPI(f,A,B,C,D,E,fp,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
    implicit double precision (A-H,O-Z)
    dimension :: f(1), fp(1)
    dimension :: A(1), B(1), C(1), D(1), E(1)
    real*8, allocatable, dimension(:) :: fm2, fm1, fp1, fp2, fsend
    Nxy = Nx*Ny
    Nyz = Ny*Nz
    allocate (fm2(Nyz), fm1(Nyz), fp1(Nyz), fp2(Nyz), fsend(Nyz))
!-------------------------------------------------------------------------------
!  if P = 1, it doesn't transfer data in MPI
!
    if (iDend == iDstart) then
      call getBCx3D(f(1),fp2(1),3,Nx,Ny,Nz)
      call getBCx3D(f(1),fp1(1),2,Nx,Ny,Nz)
      call getBCx3D(f(1),fm1(1),Nx-1,Nx,Ny,Nz)
      call getBCx3D(f(1),fm2(1),Nx-2,Nx,Ny,Nz)
      goto 123
    endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
    call getBCx3D(f(1),fsend(1),1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nyz,MPI_REAL8,iDm1,100, &
                      fp1(1)  ,Nyz,MPI_REAL8,iDp1,100, &
                      MPI_WORLD,ISTATUS,IERR)
!                      
    call getBCx3D(f(1),fsend(1),2,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nyz,MPI_REAL8,iDm1,105, &
                      fp2(1)  ,Nyz,MPI_REAL8,iDp1,105, &
                      MPI_WORLD,ISTATUS,IERR)                      
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
    call getBCx3D(f(1),fsend(1),Nx,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nyz,MPI_REAL8,iDp1,110, &
                      fm1(1)  ,Nyz,MPI_REAL8,iDm1,110, &
                      MPI_WORLD,ISTATUS,IERR)
!                      
    call getBCx3D(f(1),fsend(1),Nx-1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nyz,MPI_REAL8,iDp1,120, &
                      fm2(1)  ,Nyz,MPI_REAL8,iDm1,120, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
    if (iD == iDstart) then
      call getBCx3D(f(1),fsend(1),2,Nx,Ny,Nz)
      call MPI_SEND(fsend(1),Nyz,MPI_REAL8,iDend,130,MPI_WORLD,IERR)
!      
      call getBCx3D(f(1),fsend(1),3,Nx,Ny,Nz)
      call MPI_SEND(fsend(1),Nyz,MPI_REAL8,iDend,135,MPI_WORLD,IERR)
!      
      call MPI_RECV(fm1(1),Nyz,MPI_REAL8,iDend,140,MPI_WORLD,ISTATUS,IERR)
      call MPI_RECV(fm2(1),Nyz,MPI_REAL8,iDend,150,MPI_WORLD,ISTATUS,IERR)      
    endif
!
    if (iD == iDend) then      
      call MPI_RECV(fp1(1),Nyz,MPI_REAL8,iDstart,130,MPI_WORLD,ISTATUS,IERR)
      call MPI_RECV(fp2(1),Nyz,MPI_REAL8,iDstart,135,MPI_WORLD,ISTATUS,IERR)
!      
      call getBCx3D(f(1),fsend(1),Nx-1,Nx,Ny,Nz)
      call MPI_SEND(fsend(1),Nyz,MPI_REAL8,iDstart,140,MPI_WORLD,IERR)
!      
      call getBCx3D(f(1),fsend(1),Nx-2,Nx,Ny,Nz)
      call MPI_SEND(fsend(1),Nyz,MPI_REAL8,iDstart,150,MPI_WORLD,IERR)      
    endif
!    
123 continue                                          
!-------------------------------------------------------------------------------
    ii = 0
    do 30 k = 1, Nz
    do 30 j = 1, Ny
      ii = ii+1
      jk = (j-1)*Nx+(k-1)*Nxy
!      
      do i = 3, Nx-2
        ijkp2 = (i+2)+jk
        ijkp1 = (i+1)+jk
        ijk   = (i  )+jk
        ijkm1 = (i-1)+jk
        ijkm2 = (i-2)+jk       
!
        temp1 = A(i)*f(ijkp2)
        temp2 = B(i)*f(ijkp1)
        temp3 = C(i)*f(ijk)
        temp4 = D(i)*f(ijkm1)
        temp5 = E(i)*f(ijkm2)
        sum = DIFAB(temp1,-temp2)
        sum = DIFAB(sum,-temp3)
        sum = DIFAB(sum,-temp4)
        fp(ijk) = DIFAB(sum,-temp5)
!       fp(ijk) = A(i)*f(ijkp2)+B(i)*f(ijkp1)+C(i)*f(ijk) &
!                +D(i)*f(ijkm1)+E(i)*f(ijkm2)
      enddo      
!-------------------------------------------------------------------------------
!     for i = 1
      ijkp2 = 3+jk
      ijkp1 = 2+jk
      ijk   = 1+jk
!
      temp1 = A(1)*f(ijkp2)
      temp2 = B(1)*f(ijkp1)
      temp3 = C(1)*f(ijk)
      temp4 = D(1)*fm1(ii)
      temp5 = E(1)*fm2(ii)
      sum = DIFAB(temp1,-temp2)
      sum = DIFAB(sum,-temp3)
      sum = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)      
!     fp(ijk) = A(1)*f(ijkp2)+B(1)*f(ijkp1)+C(1)*f(ijk) &
!              +D(1)*fm1(ii)+E(1)*fm2(ii)
!-------------------------------------------------------------------------------
!     for i = 2
      i2 = 2
      ijkp2 = 4+jk
      ijkp1 = 3+jk
      ijk   = 2+jk
      ijkm1 = 1+jk
!      
      temp1 = A(i2)*f(ijkp2)
      temp2 = B(i2)*f(ijkp1)
      temp3 = C(i2)*f(ijk)
      temp4 = D(i2)*f(ijkm1)
      temp5 = E(i2)*fm1(ii)
      sum = DIFAB(temp1,-temp2)
      sum = DIFAB(sum,-temp3)
      sum = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)     
!     fp(ijk) = A(i2)*f(ijkp2)+B(i2)*f(ijkp1)+C(i2)*f(ijk) &
!              +D(i2)*f(ijkm1)+E(i2)*fm1(ii)
!-------------------------------------------------------------------------------
!     for i = Nx-1
      ijkp1 = (Nx  )+jk
      ijk   = (Nx-1)+jk
      ijkm1 = (Nx-2)+jk
      ijkm2 = (Nx-3)+jk
!
      temp1 = A(Nx-1)*fp1(ii)
      temp2 = B(Nx-1)*f(ijkp1)
      temp3 = C(Nx-1)*f(ijk)
      temp4 = D(Nx-1)*f(ijkm1)
      temp5 = E(Nx-1)*f(ijkm2)
      sum = DIFAB(temp1,-temp2)
      sum = DIFAB(sum,-temp3)
      sum = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)                          
!     fp(ijk) = A(Nx-1)*fp1(ii)+B(Nx-1)*f(ijkp1)+C(Nx-1)*f(ijk) &
!              +D(Nx-1)*f(ijkm1)+E(Nx-1)*f(ijkm2)
!-------------------------------------------------------------------------------
!     for i = Nx
      ijk   = (Nx  )+jk
      ijkm1 = (Nx-1)+jk
      ijkm2 = (Nx-2)+jk              
!
      temp1 = A(Nx)*fp2(ii)
      temp2 = B(Nx)*fp1(ii)
      temp3 = C(Nx)*f(ijk)
      temp4 = D(Nx)*f(ijkm1)
      temp5 = E(Nx)*f(ijkm2)
      sum = DIFAB(temp1,-temp2)
      sum = DIFAB(sum,-temp3)
      sum = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(Nx)*fp2(ii)+B(Nx)*fp1(ii)+C(Nx)*f(ijk) &
!              +D(Nx)*f(ijkm1)+E(Nx)*f(ijkm2) 
 30 continue
!-------------------------------------------------------------------------------
    deallocate (fm2, fm1, fp1, fp2, fsend)  
  end subroutine PCDx3rd3D_MPI
!===============================================================================
!===============================================================================
!  This is the second step of third-order central difference
!    to calculate the first derivative fp along y-direction
!    for periodic boundary condition.
!
  subroutine PCDy3rd3D_MPI(f,A,B,C,D,E,fp,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
    implicit double precision (A-H,O-Z)
    dimension :: f(1), fp(1)
    dimension :: A(1), B(1), C(1), D(1), E(1)
    real*8, allocatable, dimension(:) :: fm2, fm1, fp1, fp2, fsend
    Nxy = Nx*Ny
    Nxz = Nx*Nz
    allocate (fm2(Nxz), fm1(Nxz), fp1(Nxz), fp2(Nxz), fsend(Nxz))
!-------------------------------------------------------------------------------
!  if P = 1, it doesn't transfer data in MPI
!
    if (iDend == iDstart) then
      call getBCy3D(f(1),fp2(1),3,Nx,Ny,Nz)
      call getBCy3D(f(1),fp1(1),2,Nx,Ny,Nz)
      call getBCy3D(f(1),fm1(1),Ny-1,Nx,Ny,Nz)
      call getBCy3D(f(1),fm2(1),Ny-2,Nx,Ny,Nz)
      goto 123
    endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
    call getBCy3D(f(1),fsend(1),1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxz,MPI_REAL8,iDm1,100, &
                      fp1(1)  ,Nxz,MPI_REAL8,iDp1,100, &
                      MPI_WORLD,ISTATUS,IERR)
!                      
    call getBCy3D(f(1),fsend(1),2,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxz,MPI_REAL8,iDm1,105, &
                      fp2(1)  ,Nxz,MPI_REAL8,iDp1,105, &
                      MPI_WORLD,ISTATUS,IERR)                      
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
    call getBCy3D(f(1),fsend(1),Ny,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxz,MPI_REAL8,iDp1,110, &
                      fm1(1)  ,Nxz,MPI_REAL8,iDm1,110, &
                      MPI_WORLD,ISTATUS,IERR)
!                      
    call getBCy3D(f(1),fsend(1),Ny-1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxz,MPI_REAL8,iDp1,120, &
                      fm2(1)  ,Nxz,MPI_REAL8,iDm1,120, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
    if (iD == iDstart) then
      call getBCy3D(f(1),fsend(1),2,Nx,Ny,Nz)
      call MPI_SEND(fsend(1),Nxz,MPI_REAL8,iDend,130,MPI_WORLD,IERR)
!      
      call getBCy3D(f(1),fsend(1),3,Nx,Ny,Nz)
      call MPI_SEND(fsend(1),Nxz,MPI_REAL8,iDend,135,MPI_WORLD,IERR)
!      
      call MPI_RECV(fm1(1),Nxz,MPI_REAL8,iDend,140,MPI_WORLD,ISTATUS,IERR)
      call MPI_RECV(fm2(1),Nxz,MPI_REAL8,iDend,150,MPI_WORLD,ISTATUS,IERR)      
    endif
!
    if (iD == iDend) then      
      call MPI_RECV(fp1(1),Nxz,MPI_REAL8,iDstart,130,MPI_WORLD,ISTATUS,IERR)
      call MPI_RECV(fp2(1),Nxz,MPI_REAL8,iDstart,135,MPI_WORLD,ISTATUS,IERR)
!      
      call getBCy3D(f(1),fsend(1),Ny-1,Nx,Ny,Nz)
      call MPI_SEND(fsend(1),Nxz,MPI_REAL8,iDstart,140,MPI_WORLD,IERR)
!      
      call getBCy3D(f(1),fsend(1),Ny-2,Nx,Ny,Nz)
      call MPI_SEND(fsend(1),Nxz,MPI_REAL8,iDstart,150,MPI_WORLD,IERR)      
    endif
!    
123 continue                                          
!-------------------------------------------------------------------------------
    ii = 0
    do 30 k = 1, Nz
    do 30 i = 1, Nx
      ii = ii+1
      ik = i+(k-1)*Nxy
!      
      do j = 3, Ny-2
        ijkp2 = ((j+2)-1)*Nx+ik
        ijkp1 = ((j+1)-1)*Nx+ik
        ijk   = ((j  )-1)*Nx+ik
        ijkm1 = ((j-1)-1)*Nx+ik
        ijkm2 = ((j-2)-1)*Nx+ik       
!
        temp1 = A(j)*f(ijkp2)
        temp2 = B(j)*f(ijkp1)
        temp3 = C(j)*f(ijk)
        temp4 = D(j)*f(ijkm1)
        temp5 = E(j)*f(ijkm2)
        sum = DIFAB(temp1,-temp2)
        sum = DIFAB(sum,-temp3)
        sum = DIFAB(sum,-temp4)
        fp(ijk) = DIFAB(sum,-temp5)
!       fp(ijk) = A(j)*f(ijkp2)+B(j)*f(ijkp1)+C(j)*f(ijk) &
!                +D(j)*f(ijkm1)+E(j)*f(ijkm2)
      enddo      
!-------------------------------------------------------------------------------
!     for j = 1
      ijkp2 = (3-1)*Nx+ik
      ijkp1 = (2-1)*Nx+ik
      ijk   = (1-1)*Nx+ik
!
      temp1 = A(1)*f(ijkp2)
      temp2 = B(1)*f(ijkp1)
      temp3 = C(1)*f(ijk)
      temp4 = D(1)*fm1(ii)
      temp5 = E(1)*fm2(ii)
      sum = DIFAB(temp1,-temp2)
      sum = DIFAB(sum,-temp3)
      sum = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)        
!     fp(ijk) = A(1)*f(ijkp2)+B(1)*f(ijkp1)+C(1)*f(ijk) &
!              +D(1)*fm1(ii)+E(1)*fm2(ii)
!-------------------------------------------------------------------------------
!     for j = 2
      j2 = 2
      ijkp2 = (4-1)*Nx+ik
      ijkp1 = (3-1)*Nx+ik
      ijk   = (2-1)*Nx+ik
      ijkm1 = (1-1)*Nx+ik    
!
      temp1 = A(j2)*f(ijkp2)
      temp2 = B(j2)*f(ijkp1)
      temp3 = C(j2)*f(ijk)
      temp4 = D(j2)*f(ijkm1)
      temp5 = E(j2)*fm1(ii)
      sum = DIFAB(temp1,-temp2)
      sum = DIFAB(sum,-temp3)
      sum = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5) 
!     fp(ijk) = A(j2)*f(ijkp2)+B(j2)*f(ijkp1)+C(j2)*f(ijk) &
!              +D(j2)*f(ijkm1)+E(j2)*fm1(ii)
!-------------------------------------------------------------------------------
!     for j = Ny-1
      ijkp1 = ( Ny   -1)*Nx+ik
      ijk   = ((Ny-1)-1)*Nx+ik
      ijkm1 = ((Ny-2)-1)*Nx+ik
      ijkm2 = ((Ny-3)-1)*Nx+ik             
!
      temp1 = A(Ny-1)*fp1(ii)
      temp2 = B(Ny-1)*f(ijkp1)
      temp3 = C(Ny-1)*f(ijk)
      temp4 = D(Ny-1)*f(ijkm1)
      temp5 = E(Ny-1)*f(ijkm2)
      sum = DIFAB(temp1,-temp2)
      sum = DIFAB(sum,-temp3)
      sum = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(Ny-1)*fp1(ii)+B(Ny-1)*f(ijkp1)+C(Ny-1)*f(ijk) &
!              +D(Ny-1)*f(ijkm1)+E(Ny-1)*f(ijkm2)
!-------------------------------------------------------------------------------
!     for j = Ny
      ijk   = ((Ny  )-1)*Nx+ik
      ijkm1 = ((Ny-1)-1)*Nx+ik
      ijkm2 = ((Ny-2)-1)*Nx+ik             
!
      temp1 = A(Ny)*fp2(ii)
      temp2 = B(Ny)*fp1(ii)
      temp3 = C(Ny)*f(ijk)
      temp4 = D(Ny)*f(ijkm1)
      temp5 = E(Ny)*f(ijkm2)
      sum = DIFAB(temp1,-temp2)
      sum = DIFAB(sum,-temp3)
      sum = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(Ny)*fp2(ii)+B(Ny)*fp1(ii)+C(Ny)*f(ijk) &
!              +D(Ny)*f(ijkm1)+E(Ny)*f(ijkm2) 
 30 continue
!-------------------------------------------------------------------------------
    deallocate (fm2, fm1, fp1, fp2, fsend)  
  end subroutine PCDy3rd3D_MPI
!===============================================================================
!===============================================================================
!  This is the second step of third-order central difference
!    to calculate the first derivative fp along y-direction
!    for periodic boundary condition.
!
  subroutine PCDz3rd3D_MPI(f,A,B,C,D,E,fp,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
    implicit double precision (A-H,O-Z)
    dimension :: f(1), fp(1)
    dimension :: A(1), B(1), C(1), D(1), E(1)
    real*8, allocatable, dimension(:) :: fm2, fm1, fp1, fp2, fsend
    Nxy = Nx*Ny
    allocate (fm2(Nxy), fm1(Nxy), fp1(Nxy), fp2(Nxy), fsend(Nxy))
!-------------------------------------------------------------------------------
!  if P = 1, it doesn't transfer data in MPI
!
    if (iDend == iDstart) then
      call getBCz3D(f(1),fp2(1),3,Nx,Ny,Nz)
      call getBCz3D(f(1),fp1(1),2,Nx,Ny,Nz)
      call getBCz3D(f(1),fm1(1),Nz-1,Nx,Ny,Nz)
      call getBCz3D(f(1),fm2(1),Nz-2,Nx,Ny,Nz)
      goto 123
    endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
    call getBCz3D(f(1),fsend(1),1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxy,MPI_REAL8,iDm1,100, &
                      fp1(1)  ,Nxy,MPI_REAL8,iDp1,100, &
                      MPI_WORLD,ISTATUS,IERR)
!                      
    call getBCz3D(f(1),fsend(1),2,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxy,MPI_REAL8,iDm1,105, &
                      fp2(1)  ,Nxy,MPI_REAL8,iDp1,105, &
                      MPI_WORLD,ISTATUS,IERR)                      
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
    call getBCz3D(f(1),fsend(1),Nz,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxy,MPI_REAL8,iDp1,110, &
                      fm1(1)  ,Nxy,MPI_REAL8,iDm1,110, &
                      MPI_WORLD,ISTATUS,IERR)
!                      
    call getBCz3D(f(1),fsend(1),Nz-1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxy,MPI_REAL8,iDp1,120, &
                      fm2(1)  ,Nxy,MPI_REAL8,iDm1,120, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
    if (iD == iDstart) then
      call getBCz3D(f(1),fsend(1),2,Nx,Ny,Nz)
      call MPI_SEND(fsend(1),Nxy,MPI_REAL8,iDend,130,MPI_WORLD,IERR)
!      
      call getBCz3D(f(1),fsend(1),3,Nx,Ny,Nz)
      call MPI_SEND(fsend(1),Nxy,MPI_REAL8,iDend,135,MPI_WORLD,IERR)
!      
      call MPI_RECV(fm1(1),Nxy,MPI_REAL8,iDend,140,MPI_WORLD,ISTATUS,IERR)
      call MPI_RECV(fm2(1),Nxy,MPI_REAL8,iDend,150,MPI_WORLD,ISTATUS,IERR)      
    endif
!
    if (iD == iDend) then      
      call MPI_RECV(fp1(1),Nxy,MPI_REAL8,iDstart,130,MPI_WORLD,ISTATUS,IERR)
      call MPI_RECV(fp2(1),Nxy,MPI_REAL8,iDstart,135,MPI_WORLD,ISTATUS,IERR)
!      
      call getBCz3D(f(1),fsend(1),Nz-1,Nx,Ny,Nz)
      call MPI_SEND(fsend(1),Nxy,MPI_REAL8,iDstart,140,MPI_WORLD,IERR)
!      
      call getBCz3D(f(1),fsend(1),Nz-2,Nx,Ny,Nz)
      call MPI_SEND(fsend(1),Nxy,MPI_REAL8,iDstart,150,MPI_WORLD,IERR)      
    endif
!    
123 continue                                          
!-------------------------------------------------------------------------------
    ii = 0
    do 30 j = 1, Ny
    do 30 i = 1, Nx
      ii = ii+1
      !ijk = i+(j-1)*Nx+(k-1)*Nxy
      ij = i+(j-1)*Nx
!      
      do k = 3, Nz-2
        ijkp2 = ((k+2)-1)*Nxy+ij
        ijkp1 = ((k+1)-1)*Nxy+ij
        ijk   = ((k  )-1)*Nxy+ij
        ijkm1 = ((k-1)-1)*Nxy+ij
        ijkm2 = ((k-2)-1)*Nxy+ij       
!
        temp1 = A(k)*f(ijkp2)
        temp2 = B(k)*f(ijkp1)
        temp3 = C(k)*f(ijk)
        temp4 = D(k)*f(ijkm1)
        temp5 = E(k)*f(ijkm2)
        sum = DIFAB(temp1,-temp2)
        sum = DIFAB(sum,-temp3)
        sum = DIFAB(sum,-temp4)
        fp(ijk) = DIFAB(sum,-temp5)
!       fp(ijk) = A(k)*f(ijkp2)+B(k)*f(ijkp1)+C(k)*f(ijk) &
!                +D(k)*f(ijkm1)+E(k)*f(ijkm2)
      enddo      
!-------------------------------------------------------------------------------
!     for k = 1
      ijkp2 = (3-1)*Nxy+ij
      ijkp1 = (2-1)*Nxy+ij
      ijk   = (1-1)*Nxy+ij  
!
      temp1 = A(1)*f(ijkp2)
      temp2 = B(1)*f(ijkp1)
      temp3 = C(1)*f(ijk)
      temp4 = D(1)*fm1(ii)
      temp5 = E(1)*fm2(ii)
      sum = DIFAB(temp1,-temp2)
      sum = DIFAB(sum,-temp3)
      sum = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(1)*f(ijkp2)+B(1)*f(ijkp1)+C(1)*f(ijk) &
!              +D(1)*fm1(ii)+E(1)*fm2(ii)
!-------------------------------------------------------------------------------
!     for k = 2
      k2 = 2
      ijkp2 = (4-1)*Nxy+ij
      ijkp1 = (3-1)*Nxy+ij
      ijk   = (2-1)*Nxy+ij
      ijkm1 = (1-1)*Nxy+ij    
!
      temp1 = A(k2)*f(ijkp2)
      temp2 = B(k2)*f(ijkp1)
      temp3 = C(k2)*f(ijk)
      temp4 = D(k2)*f(ijkm1)
      temp5 = E(k2)*fm1(ii)
      sum = DIFAB(temp1,-temp2)
      sum = DIFAB(sum,-temp3)
      sum = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(k2)*f(ijkp2)+B(k2)*f(ijkp1)+C(k2)*f(ijk) &
!              +D(k2)*f(ijkm1)+E(k2)*fm1(ii)
!-------------------------------------------------------------------------------
!     for k = Nz-1
      ijkp1 = ((Nz  )-1)*Nxy+ij
      ijk   = ((Nz-1)-1)*Nxy+ij
      ijkm1 = ((Nz-2)-1)*Nxy+ij
      ijkm2 = ((Nz-3)-1)*Nxy+ij             
!
      temp1 = A(Nz-1)*fp1(ii)
      temp2 = B(Nz-1)*f(ijkp1)
      temp3 = C(Nz-1)*f(ijk)
      temp4 = D(Nz-1)*f(ijkm1)
      temp5 = E(Nz-1)*f(ijkm2)
      sum = DIFAB(temp1,-temp2)
      sum = DIFAB(sum,-temp3)
      sum = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(Nz-1)*fp1(ii)+B(Nz-1)*f(ijkp1)+C(Nz-1)*f(ijk) &
!              +D(Nz-1)*f(ijkm1)+E(Nz-1)*f(ijkm2)
!-------------------------------------------------------------------------------
!     for k = Nz
      ijk   = ((Nz  )-1)*Nxy+ij
      ijkm1 = ((Nz-1)-1)*Nxy+ij
      ijkm2 = ((Nz-2)-1)*Nxy+ij             
!
      temp1 = A(Nz)*fp2(ii)
      temp2 = B(Nz)*fp1(ii)
      temp3 = C(Nz)*f(ijk)
      temp4 = D(Nz)*f(ijkm1)
      temp5 = E(Nz)*f(ijkm2) 
      sum = DIFAB(temp1,-temp2)
      sum = DIFAB(sum,-temp3)
      sum = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(Nz)*fp2(ii)+B(Nz)*fp1(ii)+C(Nz)*f(ijk) &
!              +D(Nz)*f(ijkm1)+E(Nz)*f(ijkm2) 
 30 continue
!-------------------------------------------------------------------------------
    deallocate (fm2, fm1, fp1, fp2, fsend)  
  end subroutine PCDz3rd3D_MPI
!===============================================================================
!===============================================================================
!  This is the first step of 3rd-order central difference
!    to prepare the table A, B, C, D, E for free boundary condition
!
  subroutine FCD3rd0_MPI(h,A,B,C,D,E,N,iDm1,iD,iDp1,iDstart,iDend)
    implicit double precision (A-H,O-Z)
    dimension :: h(1), A(1), B(1), C(1), D(1), E(1)  
    real*8, allocatable, dimension(:) :: hm, hp
    allocate (hm(2), hp(1))
    i2 = 2
    i3 = 3
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
    call MPI_SENDRECV(h(1) ,1,MPI_REAL8,iDm1,100, &
                      hp(1),1,MPI_REAL8,iDp1,100, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
    call MPI_SENDRECV(h(N-1),2,MPI_REAL8,iDp1,110, &
                      hm(1) ,2,MPI_REAL8,iDm1,110, &
                      MPI_WORLD,ISTATUS,IERR)                                          
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
    if (iD == iDstart) then
      hm(1) = h(1)
      hm(2) = h(1)
    endif
!
    if (iD == iDend) then
      h(N)  = h(N)
      hp(1) = h(N)
    endif
!-------------------------------------------------------------------------------
    do i = 3, N-2
      ht1 = h(i+1)+h(i)
      ht2 = h(i)+h(i-1)
      ht3 = h(i-1)+h(i-2)
      ht4 = h(i+1)+h(i)+h(i-1)
      ht5 = h(i)+h(i-1)+h(i-2)
      ht6 = h(i+1)+h(i)+h(i-1)+h(i-2)
      A(i) = -h(i)*h(i-1)*ht3/(h(i+1)*ht1*ht4*ht6)
      B(i) =  h(i-1)*ht1*ht3/(h(i)*h(i+1)*ht2*ht5)
!
      temp1 = 1.d0/ht1
      temp2 = 1.d0/h(i)
      temp3 = 1.d0/h(i-1)
      temp4 = 1.d0/ht3
      sum   = DIFAB(temp3,temp1)
      sum   = DIFAB(sum,temp2)
      C(i)  = DIFAB(sum,-temp4)
!     C(i)  = -1.d0/ht1-1.d0/h(i)+1.d0/h(i-1)+1.d0/ht3
!
      D(i) = -h(i)*ht1*ht3/(h(i-1)*h(i-2)*ht2*ht4)
      E(i) =  h(i)*h(i-1)*ht1/(h(i-2)*ht3*ht5*ht6)
    enddo
!-------------------------------------------------------------------------------
!   for i = 1
    ht1 = h(i2)+h(1)
    ht2 = h(1)+hm(2)
    ht3 = hm(2)+hm(1)
    ht4 = h(i2)+h(1)+hm(2)
    ht5 = h(1)+hm(2)+hm(1)
    ht6 = h(i2)+h(1)+hm(2)+hm(1)
    A(1) = -h(1)*hm(2)*ht3/(h(i2)*ht1*ht4*ht6)
    B(1) =  hm(2)*ht1*ht3/(h(1)*h(i2)*ht2*ht5)
!
    temp1 = 1.d0/ht1
    temp2 = 1.d0/h(1)
    temp3 = 1.d0/hm(2)
    temp4 = 1.d0/ht3
    sum   = DIFAB(temp3,temp1)
    sum   = DIFAB(sum,temp2)
    C(1)  = DIFAB(sum,-temp4)
!   C(1)  = -1.d0/ht1-1.d0/h(1)+1.d0/hm(2)+1.d0/ht3
!
    D(1) = -h(1)*ht1*ht3/(hm(2)*hm(1)*ht2*ht4)
    E(1) =  h(1)*hm(2)*ht1/(hm(1)*ht3*ht5*ht6)
!-------------------------------------------------------------------------------
!   for i = 2
    ht1 = h(i3)+h(i2)
    ht2 = h(i2)+h(1)      
    ht3 = h(1)+hm(2)
    ht4 = h(i3)+h(i2)+h(1)      
    ht5 = h(i2)+h(1)+hm(2)      
    ht6 = h(i3)+h(i2)+h(1)+hm(2)
    A(i2) = -h(i2)*h(1)*ht3/(h(i3)*ht1*ht4*ht6)
    B(i2) =  h(1)*ht1*ht3/(h(i2)*h(i3)*ht2*ht5)
!
    temp1 = 1.d0/ht1
    temp2 = 1.d0/h(i2)
    temp3 = 1.d0/h(1)
    temp4 = 1.d0/ht3
    sum   = DIFAB(temp3,temp1)
    sum   = DIFAB(sum,temp2)
    C(i2) = DIFAB(sum,-temp4)
!   C(i2) = -1.d0/ht1-1.d0/h(i2)+1.d0/h(1)+1.d0/ht3
!
    D(i2) = -h(i2)*ht1*ht3/(h(1)*hm(2)*ht2*ht4)      
    E(i2) =  h(i2)*h(1)*ht1/(hm(2)*ht3*ht5*ht6)
!-------------------------------------------------------------------------------
!   for i = N
    ht1 = hp(1)+h(N)
    ht2 = h(N)+h(N-1)
    ht3 = h(N-1)+h(N-2)
    ht4 = hp(1)+h(N)+h(N-1)
    ht5 = h(N)+h(N-1)+h(N-2)
    ht6 = hp(1)+h(N)+h(N-1)+h(N-2)
    A(N) = -h(N)*h(N-1)*ht3/(hp(1)*ht1*ht4*ht6)
    B(N) =  h(N-1)*ht1*ht3/(h(N)*hp(1)*ht2*ht5)
!
    temp1 = 1.d0/ht1
    temp2 = 1.d0/h(N)
    temp3 = 1.d0/h(N-1)
    temp4 = 1.d0/ht3
    sum   = DIFAB(temp3,temp1)
    sum   = DIFAB(sum,temp2)
    C(N)  = DIFAB(sum,-temp4)
!   C(N)  = -1.d0/ht1-1.d0/h(N)+1.d0/h(N-1)+1.d0/ht3
!
    D(N) = -h(N)*ht1*ht3/(h(N-1)*h(N-2)*ht2*ht4)
    E(N) =  h(N)*h(N-1)*ht1/(h(N-2)*ht3*ht5*ht6)
!-------------------------------------------------------------------------------
!   for i = N-1
    ht1 = h(N)+h(N-1)
    ht2 = h(N-1)+h(N-2)
    ht3 = h(N-2)+h(N-3)
    ht4 = h(N)+h(N-1)+h(N-2)
    ht5 = h(N-1)+h(N-2)+h(N-3)
    ht6 = h(N)+h(N-1)+h(N-2)+h(N-3)
    A(N-1) = -h(N-1)*h(N-2)*ht3/(h(N)*ht1*ht4*ht6)
    B(N-1) =  h(N-2)*ht1*ht3/(h(N-1)*h(N)*ht2*ht5)
!
    temp1  = 1.d0/ht1
    temp2  = 1.d0/h(N-1)
    temp3  = 1.d0/h(N-2)
    temp4  = 1.d0/ht3
    sum    = DIFAB(temp3,temp1)
    sum    = DIFAB(sum,temp2)
    C(N-1) = DIFAB(sum,-temp4)
!   C(N-1) = -1.d0/ht1-1.d0/h(N-1)+1.d0/h(N-2)+1.d0/ht3
!
    D(N-1) = -h(N-1)*ht1*ht3/(h(N-2)*h(N-3)*ht2*ht4)
    E(N-1) =  h(N-1)*h(N-2)*ht1/(h(N-3)*ht3*ht5*ht6)
!-------------------------------------------------------------------------------
    deallocate (hm, hp)
  end subroutine FCD3rd0_MPI  
!===============================================================================
!===============================================================================
!  This is the second step of third-order central difference
!    to calculate the first derivative fp along x-direction
!    for free boundary condition.
!
  subroutine FCDx3rd3D_MPI(f,A,B,C,D,E,fp,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
    implicit double precision (A-H,O-Z)
    dimension :: f(1), fp(1)
    dimension :: A(1), B(1), C(1), D(1), E(1)
    real*8, allocatable, dimension(:) :: fm2, fm1, fp1, fp2, fsend
    Nxy = Nx*Ny
    Nyz = Ny*Nz
    allocate (fm2(Nyz), fm1(Nyz), fp1(Nyz), fp2(Nyz), fsend(Nyz))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
    call getBCx3D(f(1),fsend(1),1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nyz,MPI_REAL8,iDm1,100, &
                      fp1(1)  ,Nyz,MPI_REAL8,iDp1,100, &
                      MPI_WORLD,ISTATUS,IERR)
!                      
    call getBCx3D(f(1),fsend(1),2,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nyz,MPI_REAL8,iDm1,105, &
                      fp2(1)  ,Nyz,MPI_REAL8,iDp1,105, &
                      MPI_WORLD,ISTATUS,IERR)                      
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
    call getBCx3D(f(1),fsend(1),Nx,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nyz,MPI_REAL8,iDp1,110, &
                      fm1(1)  ,Nyz,MPI_REAL8,iDm1,110, &
                      MPI_WORLD,ISTATUS,IERR)
!                      
    call getBCx3D(f(1),fsend(1),Nx-1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nyz,MPI_REAL8,iDp1,120, &
                      fm2(1)  ,Nyz,MPI_REAL8,iDm1,120, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
    if (iD == iDstart) then
      call getBCx3D(f(1),fm1(1),1,Nx,Ny,Nz)      
      call getBCx3D(f(1),fm2(1),1,Nx,Ny,Nz)      
    endif
!
    if (iD == iDend) then           
      call getBCx3D(f(1),fp1(1),Nx,Nx,Ny,Nz)      
      call getBCx3D(f(1),fp2(1),Nx,Nx,Ny,Nz)      
    endif                                         
!-------------------------------------------------------------------------------
    ii = 0
    do 30 k = 1, Nz
    do 30 j = 1, Ny
      ii = ii+1
      jk = (j-1)*Nx+(k-1)*Nxy
!      
      do i = 3, Nx-2
        ijkp2 = (i+2)+jk
        ijkp1 = (i+1)+jk
        ijk   = (i  )+jk
        ijkm1 = (i-1)+jk
        ijkm2 = (i-2)+jk       
!
        temp1 = A(i)*f(ijkp2)
        temp2 = B(i)*f(ijkp1)
        temp3 = C(i)*f(ijk)
        temp4 = D(i)*f(ijkm1)
        temp5 = E(i)*f(ijkm2)
        sum   = DIFAB(temp1,-temp2)
        sum   = DIFAB(sum,-temp3)
        sum   = DIFAB(sum,-temp4)
        fp(ijk) = DIFAB(sum,-temp5)
!       fp(ijk) = A(i)*f(ijkp2)+B(i)*f(ijkp1)+C(i)*f(ijk) &
!                +D(i)*f(ijkm1)+E(i)*f(ijkm2)
      enddo      
!-------------------------------------------------------------------------------
!     for i = 1
      ijkp2 = 3+jk
      ijkp1 = 2+jk
      ijk   = 1+jk      
!
      temp1 = A(1)*f(ijkp2)
      temp2 = B(1)*f(ijkp1)
      temp3 = C(1)*f(ijk)
      temp4 = D(1)*fm1(ii)
      temp5 = E(1)*fm2(ii)
      sum   = DIFAB(temp1,-temp2)
      sum   = DIFAB(sum,-temp3)
      sum   = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(1)*f(ijkp2)+B(1)*f(ijkp1)+C(1)*f(ijk) &
!              +D(1)*fm1(ii)+E(1)*fm2(ii)
! 
      if (iD == iDstart) fp(ijk) = 0.d0
!-------------------------------------------------------------------------------
!     for i = 2
      i2 = 2
      ijkp2 = 4+jk
      ijkp1 = 3+jk
      ijk   = 2+jk
      ijkm1 = 1+jk     
!
      temp1 = A(i2)*f(ijkp2)
      temp2 = B(i2)*f(ijkp1)
      temp3 = C(i2)*f(ijk)
      temp4 = D(i2)*f(ijkm1)
      temp5 = E(i2)*fm1(ii)
      sum = DIFAB(temp1,-temp2)
      sum = DIFAB(sum,-temp3)
      sum = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(i2)*f(ijkp2)+B(i2)*f(ijkp1)+C(i2)*f(ijk) &
!              +D(i2)*f(ijkm1)+E(i2)*fm1(ii)
!-------------------------------------------------------------------------------
!     for i = Nx-1
      ijkp1 = (Nx  )+jk
      ijk   = (Nx-1)+jk
      ijkm1 = (Nx-2)+jk
      ijkm2 = (Nx-3)+jk              
!
      temp1 = A(Nx-1)*fp1(ii)
      temp2 = B(Nx-1)*f(ijkp1)
      temp3 = C(Nx-1)*f(ijk)
      temp4 = D(Nx-1)*f(ijkm1)
      temp5 = E(Nx-1)*f(ijkm2)
      sum = DIFAB(temp1,-temp2)
      sum = DIFAB(sum,-temp3)
      sum = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(Nx-1)*fp1(ii)+B(Nx-1)*f(ijkp1)+C(Nx-1)*f(ijk) &
!              +D(Nx-1)*f(ijkm1)+E(Nx-1)*f(ijkm2)
!-------------------------------------------------------------------------------
!     for i = Nx
      ijk   = (Nx  )+jk
      ijkm1 = (Nx-1)+jk
      ijkm2 = (Nx-2)+jk              
!
      temp1 = A(Nx)*fp2(ii)
      temp2 = B(Nx)*fp1(ii)
      temp3 = C(Nx)*f(ijk)
      temp4 = D(Nx)*f(ijkm1)
      temp5 = E(Nx)*f(ijkm2)
      sum = DIFAB(temp1,-temp2)
      sum = DIFAB(sum,-temp3)
      sum = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(Nx)*fp2(ii)+B(Nx)*fp1(ii)+C(Nx)*f(ijk) &
!              +D(Nx)*f(ijkm1)+E(Nx)*f(ijkm2) 
!              
      if (iD == iDend) fp(ijk) = 0.d0
 30 continue
!-------------------------------------------------------------------------------
    deallocate (fm2, fm1, fp1, fp2, fsend)  
  end subroutine FCDx3rd3D_MPI
!===============================================================================
!===============================================================================
!  This is the second step of third-order central difference
!    to calculate the first derivative fp along y-direction
!    for free boundary condition.
!
  subroutine FCDy3rd3D_MPI(f,A,B,C,D,E,fp,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
    implicit double precision (A-H,O-Z)
    dimension :: f(1), fp(1)
    dimension :: A(1), B(1), C(1), D(1), E(1)
    real*8, allocatable, dimension(:) :: fm2, fm1, fp1, fp2, fsend
    Nxy = Nx*Ny
    Nxz = Nx*Nz
    allocate (fm2(Nxz), fm1(Nxz), fp1(Nxz), fp2(Nxz), fsend(Nxz))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
    call getBCy3D(f(1),fsend(1),1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxz,MPI_REAL8,iDm1,100, &
                      fp1(1)  ,Nxz,MPI_REAL8,iDp1,100, &
                      MPI_WORLD,ISTATUS,IERR)
!                      
    call getBCy3D(f(1),fsend(1),2,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxz,MPI_REAL8,iDm1,105, &
                      fp2(1)  ,Nxz,MPI_REAL8,iDp1,105, &
                      MPI_WORLD,ISTATUS,IERR)                      
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
    call getBCy3D(f(1),fsend(1),Ny,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxz,MPI_REAL8,iDp1,110, &
                      fm1(1)  ,Nxz,MPI_REAL8,iDm1,110, &
                      MPI_WORLD,ISTATUS,IERR)
!                      
    call getBCy3D(f(1),fsend(1),Ny-1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxz,MPI_REAL8,iDp1,120, &
                      fm2(1)  ,Nxz,MPI_REAL8,iDm1,120, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
    if (iD == iDstart) then
      call getBCy3D(f(1),fm1(1),1,Nx,Ny,Nz)      
      call getBCy3D(f(1),fm2(1),1,Nx,Ny,Nz)      
    endif
!
    if (iD == iDend) then           
      call getBCy3D(f(1),fp1(1),Ny,Nx,Ny,Nz)      
      call getBCy3D(f(1),fp2(1),Ny,Nx,Ny,Nz)      
    endif                                           
!-------------------------------------------------------------------------------
    ii = 0
    do 30 k = 1, Nz
    do 30 i = 1, Nx
      ii = ii+1
      ik = i+(k-1)*Nxy
!      
      do j = 3, Ny-2
        ijkp2 = ((j+2)-1)*Nx+ik
        ijkp1 = ((j+1)-1)*Nx+ik
        ijk   = ((j  )-1)*Nx+ik
        ijkm1 = ((j-1)-1)*Nx+ik
        ijkm2 = ((j-2)-1)*Nx+ik       
!
        temp1 = A(j)*f(ijkp2)
        temp2 = B(j)*f(ijkp1)
        temp3 = C(j)*f(ijk)
        temp4 = D(j)*f(ijkm1)
        temp5 = E(j)*f(ijkm2)
        sum   = DIFAB(temp1,-temp2)
        sum   = DIFAB(sum,-temp3)
        sum   = DIFAB(sum,-temp4)
        fp(ijk) = DIFAB(sum,-temp5)
!       fp(ijk) = A(j)*f(ijkp2)+B(j)*f(ijkp1)+C(j)*f(ijk) &
!                +D(j)*f(ijkm1)+E(j)*f(ijkm2)
      enddo      
!-------------------------------------------------------------------------------
!     for j = 1
      ijkp2 = (3-1)*Nx+ik
      ijkp1 = (2-1)*Nx+ik
      ijk   = (1-1)*Nx+ik  
!
      temp1 = A(1)*f(ijkp2)
      temp2 = B(1)*f(ijkp1)
      temp3 = C(1)*f(ijk)
      temp4 = D(1)*fm1(ii)
      temp5 = E(1)*fm2(ii)
      sum   = DIFAB(temp1,-temp2)
      sum   = DIFAB(sum,-temp3)
      sum   = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(1)*f(ijkp2)+B(1)*f(ijkp1)+C(1)*f(ijk) &
!              +D(1)*fm1(ii)+E(1)*fm2(ii)
!
      if (iD == iDstart) fp(ijk) = 0.d0
!-------------------------------------------------------------------------------
!     for j = 2
      j2 = 2
      ijkp2 = (4-1)*Nx+ik
      ijkp1 = (3-1)*Nx+ik
      ijk   = (2-1)*Nx+ik
      ijkm1 = (1-1)*Nx+ik    
!
      temp1 = A(j2)*f(ijkp2)
      temp2 = B(j2)*f(ijkp1)
      temp3 = C(j2)*f(ijk)
      temp4 = D(j2)*f(ijkm1)
      temp5 = E(j2)*fm1(ii)
      sum   = DIFAB(temp1,-temp2)
      sum   = DIFAB(sum,-temp3)
      sum   = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(j2)*f(ijkp2)+B(j2)*f(ijkp1)+C(j2)*f(ijk) &
!              +D(j2)*f(ijkm1)+E(j2)*fm1(ii)
!-------------------------------------------------------------------------------
!     for j = Ny-1
      ijkp1 = ( Ny   -1)*Nx+ik
      ijk   = ((Ny-1)-1)*Nx+ik
      ijkm1 = ((Ny-2)-1)*Nx+ik
      ijkm2 = ((Ny-3)-1)*Nx+ik             
!
      temp1 = A(Ny-1)*fp1(ii)
      temp2 = B(Ny-1)*f(ijkp1)
      temp3 = C(Ny-1)*f(ijk)
      temp4 = D(Ny-1)*f(ijkm1)
      temp5 = E(Ny-1)*f(ijkm2)
      sum   = DIFAB(temp1,-temp2)
      sum   = DIFAB(sum,-temp3)
      sum   = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(Ny-1)*fp1(ii)+B(Ny-1)*f(ijkp1)+C(Ny-1)*f(ijk) &
!              +D(Ny-1)*f(ijkm1)+E(Ny-1)*f(ijkm2)
!-------------------------------------------------------------------------------
!     for j = Ny
      ijk   = ((Ny  )-1)*Nx+ik
      ijkm1 = ((Ny-1)-1)*Nx+ik
      ijkm2 = ((Ny-2)-1)*Nx+ik             
!
      temp1 = A(Ny)*fp2(ii)
      temp2 = B(Ny)*fp1(ii)
      temp3 = C(Ny)*f(ijk)
      temp4 = D(Ny)*f(ijkm1)
      temp5 = E(Ny)*f(ijkm2) 
      sum   = DIFAB(temp1,-temp2)
      sum   = DIFAB(sum,-temp3)
      sum   = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(Ny)*fp2(ii)+B(Ny)*fp1(ii)+C(Ny)*f(ijk) &
!              +D(Ny)*f(ijkm1)+E(Ny)*f(ijkm2) 
!
      if (iD == iDend) fp(ijk) = 0.d0
 30 continue
!-------------------------------------------------------------------------------
    deallocate (fm2, fm1, fp1, fp2, fsend)  
  end subroutine FCDy3rd3D_MPI
!===============================================================================
!===============================================================================
!  This is the second step of third-order central difference
!    to calculate the first derivative fp along y-direction
!    for periodic boundary condition.
!
  subroutine FCDz3rd3D_MPI(f,A,B,C,D,E,fp,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
    implicit double precision (A-H,O-Z)
    dimension :: f(1), fp(1)
    dimension :: A(1), B(1), C(1), D(1), E(1)
    real*8, allocatable, dimension(:) :: fm2, fm1, fp1, fp2, fsend
    Nxy = Nx*Ny
    allocate (fm2(Nxy), fm1(Nxy), fp1(Nxy), fp2(Nxy), fsend(Nxy))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
    call getBCz3D(f(1),fsend(1),1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxy,MPI_REAL8,iDm1,100, &
                      fp1(1)  ,Nxy,MPI_REAL8,iDp1,100, &
                      MPI_WORLD,ISTATUS,IERR)
!                      
    call getBCz3D(f(1),fsend(1),2,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxy,MPI_REAL8,iDm1,105, &
                      fp2(1)  ,Nxy,MPI_REAL8,iDp1,105, &
                      MPI_WORLD,ISTATUS,IERR)                      
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
    call getBCz3D(f(1),fsend(1),Nz,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxy,MPI_REAL8,iDp1,110, &
                      fm1(1)  ,Nxy,MPI_REAL8,iDm1,110, &
                      MPI_WORLD,ISTATUS,IERR)
!                      
    call getBCz3D(f(1),fsend(1),Nz-1,Nx,Ny,Nz)
    call MPI_SENDRECV(fsend(1),Nxy,MPI_REAL8,iDp1,120, &
                      fm2(1)  ,Nxy,MPI_REAL8,iDm1,120, &
                      MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
    if (iD == iDstart) then
      call getBCz3D(f(1),fm1(1),1,Nx,Ny,Nz)      
      call getBCz3D(f(1),fm2(1),1,Nx,Ny,Nz)      
    endif
!
    if (iD == iDend) then           
      call getBCz3D(f(1),fp1(1),Nz,Nx,Ny,Nz)      
      call getBCz3D(f(1),fp2(1),Nz,Nx,Ny,Nz)      
    endif                                           
!-------------------------------------------------------------------------------
    ii = 0
    do 30 j = 1, Ny
    do 30 i = 1, Nx
      ii = ii+1
      !ijk = i+(j-1)*Nx+(k-1)*Nxy
      ij = i+(j-1)*Nx
!      
      do k = 3, Nz-2
        ijkp2 = ((k+2)-1)*Nxy+ij
        ijkp1 = ((k+1)-1)*Nxy+ij
        ijk   = ((k  )-1)*Nxy+ij
        ijkm1 = ((k-1)-1)*Nxy+ij
        ijkm2 = ((k-2)-1)*Nxy+ij       
!
        temp1 = A(k)*f(ijkp2)
        temp2 = B(k)*f(ijkp1)
        temp3 = C(k)*f(ijk)
        temp4 = D(k)*f(ijkm1)
        temp5 = E(k)*f(ijkm2)
        sum   = DIFAB(temp1,-temp2)
        sum   = DIFAB(sum,-temp3)
        sum   = DIFAB(sum,-temp4)
        fp(ijk) = DIFAB(sum,-temp5)
!       fp(ijk) = A(k)*f(ijkp2)+B(k)*f(ijkp1)+C(k)*f(ijk) &
!                +D(k)*f(ijkm1)+E(k)*f(ijkm2)
      enddo      
!-------------------------------------------------------------------------------
!     for k = 1
      ijkp2 = (3-1)*Nxy+ij
      ijkp1 = (2-1)*Nxy+ij
      ijk   = (1-1)*Nxy+ij  
!
      temp1 = A(1)*f(ijkp2)
      temp2 = B(1)*f(ijkp1)
      temp3 = C(1)*f(ijk)
      temp4 = D(1)*fm1(ii)
      temp5 = E(1)*fm2(ii)
      sum   = DIFAB(temp1,-temp2)
      sum   = DIFAB(sum,-temp3)
      sum   = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(1)*f(ijkp2)+B(1)*f(ijkp1)+C(1)*f(ijk) &
!              +D(1)*fm1(ii)+E(1)*fm2(ii)
!
      if (iD == iDstart) fp(ijk) = 0.d0
!-------------------------------------------------------------------------------
!     for k = 2
      k2 = 2
      ijkp2 = (4-1)*Nxy+ij
      ijkp1 = (3-1)*Nxy+ij
      ijk   = (2-1)*Nxy+ij
      ijkm1 = (1-1)*Nxy+ij    
!
      temp1 = A(k2)*f(ijkp2)
      temp2 = B(k2)*f(ijkp1)
      temp3 = C(k2)*f(ijk)
      temp4 = D(k2)*f(ijkm1)
      temp5 = E(k2)*fm1(ii)
      sum   = DIFAB(temp1,-temp2)
      sum   = DIFAB(sum,-temp3)
      sum   = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(k2)*f(ijkp2)+B(k2)*f(ijkp1)+C(k2)*f(ijk) &
!              +D(k2)*f(ijkm1)+E(k2)*fm1(ii)
!-------------------------------------------------------------------------------
!     for k = Nz-1
      ijkp1 = ((Nz  )-1)*Nxy+ij
      ijk   = ((Nz-1)-1)*Nxy+ij
      ijkm1 = ((Nz-2)-1)*Nxy+ij
      ijkm2 = ((Nz-3)-1)*Nxy+ij             
!
      temp1 = A(Nz-1)*fp1(ii)
      temp2 = B(Nz-1)*f(ijkp1)
      temp3 = C(Nz-1)*f(ijk)
      temp4 = D(Nz-1)*f(ijkm1)
      temp5 = E(Nz-1)*f(ijkm2)
      sum   = DIFAB(temp1,-temp2)
      sum   = DIFAB(sum,-temp3)
      sum   = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(Nz-1)*fp1(ii)+B(Nz-1)*f(ijkp1)+C(Nz-1)*f(ijk) &
!              +D(Nz-1)*f(ijkm1)+E(Nz-1)*f(ijkm2)
!-------------------------------------------------------------------------------
!     for k = Nz
      ijk   = ((Nz  )-1)*Nxy+ij
      ijkm1 = ((Nz-1)-1)*Nxy+ij
      ijkm2 = ((Nz-2)-1)*Nxy+ij             
!
      temp1 = A(Nz)*fp2(ii)
      temp2 = B(Nz)*fp1(ii)
      temp3 = C(Nz)*f(ijk)
      temp4 = D(Nz)*f(ijkm1)
      temp5 = E(Nz)*f(ijkm2)
      sum   = DIFAB(temp1,-temp2)
      sum   = DIFAB(sum,-temp3)
      sum   = DIFAB(sum,-temp4)
      fp(ijk) = DIFAB(sum,-temp5)
!     fp(ijk) = A(Nz)*fp2(ii)+B(Nz)*fp1(ii)+C(Nz)*f(ijk) &
!              +D(Nz)*f(ijkm1)+E(Nz)*f(ijkm2) 
!
      if (iD == iDend) fp(ijk) = 0.d0 
 30 continue
!-------------------------------------------------------------------------------
    deallocate (fm2, fm1, fp1, fp2, fsend)  
  end subroutine FCDz3rd3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fifth-order central difference
!    to calculate the first derivative fp along x-direction
!    for free boundary condition.
!-------------------------------------------------------------------------------
   subroutine FCDx5th3D_MPI(f,fp,dx,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy = Nx*Ny
     Nyz = Ny*Nz
     Nyz2 = 2*Nyz
     Nyz3 = 3*Nyz
     allocate (fm3(Nyz),fm2(Nyz),fm1(Nyz),fp1(Nyz),fp2(Nyz),fp3(Nyz))
     allocate (fsend(Nyz3),frecvp(Nyz3),frecvm(Nyz3))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCx3D(f(1),fsend(     1),1,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz +1),2,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nyz3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fp1(1:Nyz) = frecvp(1     :Nyz )
     fp2(1:Nyz) = frecvp(Nyz +1:Nyz2)
     fp3(1:Nyz) = frecvp(Nyz2+1:Nyz3)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCx3D(f(1),fsend(     1),Nx,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz +1),Nx-1,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz2+1),Nx-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nyz3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fm1(1:Nyz) = frecvm(1     :Nyz )
     fm2(1:Nyz) = frecvm(Nyz +1:Nyz2)
     fm3(1:Nyz) = frecvm(Nyz2+1:Nyz3)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
     if (iD == iDstart) then
       call getBCx3D(f(1),fm1(1),1,Nx,Ny,Nz)      
       call getBCx3D(f(1),fm2(1),1,Nx,Ny,Nz)
       call getBCx3D(f(1),fm3(1),1,Nx,Ny,Nz)      
     endif
!
     if (iD == iDend) then           
       call getBCx3D(f(1),fp1(1),Nx,Nx,Ny,Nz)      
       call getBCx3D(f(1),fp2(1),Nx,Nx,Ny,Nz)
       call getBCx3D(f(1),fp3(1),Nx,Nx,Ny,Nz)      
     endif                                         
!-------------------------------------------------------------------------------
     ii = 0
!$OMP parallel do private(ii,jk,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2) 
     do 30 k = 1, Nz
     do 30 j = 1, Ny
       ii = j+(k-1)*Ny
       jk = (j-1)*Nx+(k-1)*Nxy
!      
       do i = 4, Nx-3
         ijkp3 = (i+3)+jk
         ijkp2 = (i+2)+jk
         ijkp1 = (i+1)+jk
         ijk   = (i  )+jk
         ijkm1 = (i-1)+jk
         ijkm2 = (i-2)+jk
         ijkm3 = (i-3)+jk
         fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dx)
       enddo      
!-------------------------------------------------------------------------------
!     for i = 1
       ijkp3 = 4+jk
       ijkp2 = 3+jk
       ijkp1 = 2+jk
       ijk   = 1+jk
!
       fp(ijk) = FCD5th(fm3(ii),fm2(ii),fm1(ii),f(ijkp1),f(ijkp2),f(ijkp3),dx)
!-------------------------------------------------------------------------------
!      for i = 2
       ijkp3 = 5+jk
       ijkp2 = 4+jk
       ijkp1 = 3+jk
       ijk   = 2+jk
       ijkm1 = 1+jk     
!
       fp(ijk) = FCD5th(fm2(ii),fm1(ii),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dx)
!-------------------------------------------------------------------------------
!      for i = 3
       ijkp3 = 6+jk
       ijkp2 = 5+jk
       ijkp1 = 4+jk
       ijk   = 3+jk
       ijkm1 = 2+jk
       ijkm2 = 1+jk     
!
       fp(ijk) = FCD5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dx)
!-------------------------------------------------------------------------------
!      for i = Nx-2
       ijkp2 = (Nx  )+jk
       ijkp1 = (Nx-1)+jk
       ijk   = (Nx-2)+jk
       ijkm1 = (Nx-3)+jk
       ijkm2 = (Nx-4)+jk
       ijkm3 = (Nx-5)+jk              
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),fp1(ii),dx)
!-------------------------------------------------------------------------------
!      for i = Nx-1
       ijkp1 = (Nx  )+jk
       ijk   = (Nx-1)+jk
       ijkm1 = (Nx-2)+jk
       ijkm2 = (Nx-3)+jk
       ijkm3 = (Nx-4)+jk              
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),fp1(ii),fp2(ii),dx)
!-------------------------------------------------------------------------------
!      for i = Nx
       ijk   = (Nx  )+jk
       ijkm1 = (Nx-1)+jk
       ijkm2 = (Nx-2)+jk
       ijkm3 = (Nx-3)+jk              
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),fp1(ii),fp2(ii),fp3(ii),dx)
 30  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm)  
   end subroutine FCDx5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fifth-order central difference
!    to calculate the first derivative fp along x-direction
!    for periodic boundary condition.
!-------------------------------------------------------------------------------
   subroutine PCDx5th3D_MPI(f,fp,dx,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input data    : dx
!  input arrays  : f(Nxyz)
!  output arrays : fp(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy  = Nx*Ny
     Nyz  = Ny*Nz
     Nyz2 = 2*Nyz
     Nyz3 = 3*Nyz
     allocate (fm3(Nyz),fm2(Nyz),fm1(Nyz),fp1(Nyz),fp2(Nyz),fp3(Nyz))
     allocate (fsend(Nyz3),frecvp(Nyz3),frecvm(Nyz3))
!-------------------------------------------------------------------------------
!  if P = 1, it doesn't transfer data in MPI
!
     if (iDend == iDstart) then       
       call getBCx3D(f(1),fp3(1),4,Nx,Ny,Nz)
       call getBCx3D(f(1),fp2(1),3,Nx,Ny,Nz)
       call getBCx3D(f(1),fp1(1),2,Nx,Ny,Nz)
!       
       call getBCx3D(f(1),fm1(1),Nx-1,Nx,Ny,Nz)
       call getBCx3D(f(1),fm2(1),Nx-2,Nx,Ny,Nz)
       call getBCx3D(f(1),fm3(1),Nx-3,Nx,Ny,Nz)
!
       do 20 k = 1, Nz
       do 20 j = 1, Ny
         ii = j+(k-1)*Ny
         jk = (j-1)*Nx+(k-1)*Nxy
! 
         ijk1 = 1 +jk
         ijkN = Nx+jk
         f(ijkN) = f(ijk1)
  20   continue
       goto 123
     endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCx3D(f(1),fsend(1     ),1,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz +1),2,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz2+1),3,Nx,Ny,Nz)
!
     call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nyz3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR) 
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCx3D(f(1),fsend(1     ),Nx  ,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz +1),Nx-1,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz2+1),Nx-2,Nx,Ny,Nz)
!
     call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nyz3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR) 
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
     if (iD == iDstart) then
!
       call getBCx3D(f(1),fsend(1     ),2,Nx,Ny,Nz)
       call getBCx3D(f(1),fsend(Nyz +1),3,Nx,Ny,Nz)
       call getBCx3D(f(1),fsend(Nyz2+1),4,Nx,Ny,Nz)
!
       call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDend,120, &
                         frecvm(1),Nyz3,MPI_REAL8,iDend,120, &
                         MPI_WORLD,ISTATUS,IERR)						 
     endif
!
     if (iD == iDend) then      
!
       call getBCx3D(f(1),fsend(1     ),Nx-1,Nx,Ny,Nz)
       call getBCx3D(f(1),fsend(Nyz +1),Nx-2,Nx,Ny,Nz)
       call getBCx3D(f(1),fsend(Nyz2+1),Nx-3,Nx,Ny,Nz)
!
       call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDstart,120, &
                         frecvp(1),Nyz3,MPI_REAL8,iDstart,120, &
                         MPI_WORLD,ISTATUS,IERR)
     endif
!
     fp1(1:Nyz) = frecvp(1     :Nyz )
     fp2(1:Nyz) = frecvp(Nyz +1:Nyz2)
     fp3(1:Nyz) = frecvp(Nyz2+1:Nyz3)
     fm1(1:Nyz) = frecvm(1     :Nyz )
     fm2(1:Nyz) = frecvm(Nyz +1:Nyz2)
     fm3(1:Nyz) = frecvm(Nyz2+1:Nyz3)
!    
123  continue                                          
!-------------------------------------------------------------------------------
     ii = 0
!$OMP parallel do private(ii,jk,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 k = 1, Nz
     do 30 j = 1, Ny
       ii = j+(k-1)*Ny
       jk = (j-1)*Nx+(k-1)*Nxy
!      
       do i = 4, Nx-3
         ijkp3 = (i+3)+jk
         ijkp2 = (i+2)+jk
         ijkp1 = (i+1)+jk
         ijk   = (i  )+jk
         ijkm1 = (i-1)+jk
         ijkm2 = (i-2)+jk
         ijkm3 = (i-3)+jk
         fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dx)
       enddo      
!-------------------------------------------------------------------------------
!     for i = 1
       ijkp3 = 4+jk
       ijkp2 = 3+jk
       ijkp1 = 2+jk
       ijk   = 1+jk
!
       fp(ijk) = FCD5th(fm3(ii),fm2(ii),fm1(ii),f(ijkp1),f(ijkp2),f(ijkp3),dx)
!-------------------------------------------------------------------------------
!      for i = 2
       ijkp3 = 5+jk
       ijkp2 = 4+jk
       ijkp1 = 3+jk
       ijk   = 2+jk
       ijkm1 = 1+jk     
!
       fp(ijk) = FCD5th(fm2(ii),fm1(ii),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dx)
!-------------------------------------------------------------------------------
!      for i = 3
       ijkp3 = 6+jk
       ijkp2 = 5+jk
       ijkp1 = 4+jk
       ijk   = 3+jk
       ijkm1 = 2+jk
       ijkm2 = 1+jk     
!
       fp(ijk) = FCD5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dx)
!-------------------------------------------------------------------------------
!      for i = Nx-2
       ijkp2 = (Nx  )+jk
       ijkp1 = (Nx-1)+jk
       ijk   = (Nx-2)+jk
       ijkm1 = (Nx-3)+jk
       ijkm2 = (Nx-4)+jk
       ijkm3 = (Nx-5)+jk              
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),fp1(ii),dx)
!-------------------------------------------------------------------------------
!      for i = Nx-1
       ijkp1 = (Nx  )+jk
       ijk   = (Nx-1)+jk
       ijkm1 = (Nx-2)+jk
       ijkm2 = (Nx-3)+jk
       ijkm3 = (Nx-4)+jk              
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),fp1(ii),fp2(ii),dx)
!-------------------------------------------------------------------------------
!      for i = Nx
       ijk   = (Nx  )+jk
       ijkm1 = (Nx-1)+jk
       ijkm2 = (Nx-2)+jk
       ijkm3 = (Nx-3)+jk              
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),fp1(ii),fp2(ii),fp3(ii),dx)
 30  continue
!$OMP end parallel do 
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm) 
   end subroutine PCDx5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fifth-order central difference
!    to calculate the first derivative fp along y-direction
!    for free boundary condition.
!-------------------------------------------------------------------------------
   subroutine FCDy5th3D_MPI(f,fp,dy,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input data    : dy
!  input arrays  : f(Nxyz)
!  output arrays : fp(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy  = Nx*Ny
     Nxz  = Nx*Nz
     Nxz2 = 2*Nxz
     Nxz3 = 3*Nxz
     allocate (fm3(Nxz),fm2(Nxz),fm1(Nxz),fp1(Nxz),fp2(Nxz),fp3(Nxz))
     allocate (fsend(Nxz3),frecvp(Nxz3),frecvm(Nxz3))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCy3D(f(1),fsend(     1),1,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz +1),2,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxz3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nxz3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fp1(1:Nxz) = frecvp(1     :Nxz )
     fp2(1:Nxz) = frecvp(Nxz +1:Nxz2)
     fp3(1:Nxz) = frecvp(Nxz2+1:Nxz3)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCy3D(f(1),fsend(     1),Ny  ,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz +1),Ny-1,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz2+1),Ny-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxz3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nxz3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fm1(1:Nxz) = frecvm(1     :Nxz )
     fm2(1:Nxz) = frecvm(Nxz +1:Nxz2)
     fm3(1:Nxz) = frecvm(Nxz2+1:Nxz3)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
     if (iD == iDstart) then
       call getBCy3D(f(1),fm1(1),1,Nx,Ny,Nz)      
       call getBCy3D(f(1),fm2(1),1,Nx,Ny,Nz)
       call getBCy3D(f(1),fm3(1),1,Nx,Ny,Nz)     
     endif
!
     if (iD == iDend) then           
       call getBCy3D(f(1),fp1(1),Ny,Nx,Ny,Nz)      
       call getBCy3D(f(1),fp2(1),Ny,Nx,Ny,Nz)
       call getBCy3D(f(1),fp3(1),Ny,Nx,Ny,Nz)      
     endif                                           
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,ik,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 k = 1, Nz
     do 30 i = 1, Nx
       ii = i+(k-1)*Nx
       ik = i+(k-1)*Nxy
!      
       do j = 4, Ny-3
         ijkp3 = ((j+3)-1)*Nx+ik
         ijkp2 = ((j+2)-1)*Nx+ik
         ijkp1 = ((j+1)-1)*Nx+ik
         ijk   = ((j  )-1)*Nx+ik
         ijkm1 = ((j-1)-1)*Nx+ik
         ijkm2 = ((j-2)-1)*Nx+ik
         ijkm3 = ((j-3)-1)*Nx+ik       
!
         fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dy)
       enddo      
!-------------------------------------------------------------------------------
!      for j = 1
       ijkp3 = (4-1)*Nx+ik
       ijkp2 = (3-1)*Nx+ik
       ijkp1 = (2-1)*Nx+ik
       ijk   = (1-1)*Nx+ik  
!
       fp(ijk) = FCD5th(fm3(ii),fm2(ii),fm1(ii),f(ijkp1),f(ijkp2),f(ijkp3),dy)
!-------------------------------------------------------------------------------
!      for j = 2
       ijkp3 = (5-1)*Nx+ik
       ijkp2 = (4-1)*Nx+ik
       ijkp1 = (3-1)*Nx+ik
       ijk   = (2-1)*Nx+ik
       ijkm1 = (1-1)*Nx+ik    
!
       fp(ijk) = FCD5th(fm2(ii),fm1(ii),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dy)
!-------------------------------------------------------------------------------
!      for j = 3
       ijkp3 = (6-1)*Nx+ik
       ijkp2 = (5-1)*Nx+ik
       ijkp1 = (4-1)*Nx+ik
       ijk   = (3-1)*Nx+ik
       ijkm1 = (2-1)*Nx+ik
       ijkm2 = (1-1)*Nx+ik    
!
       fp(ijk) = FCD5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dy)
!-------------------------------------------------------------------------------
!      for j = Ny-2
       ijkp2 = ((Ny  )-1)*Nx+ik
       ijkp1 = ((Ny-1)-1)*Nx+ik
       ijk   = ((Ny-2)-1)*Nx+ik
       ijkm1 = ((Ny-3)-1)*Nx+ik
       ijkm2 = ((Ny-4)-1)*Nx+ik
       ijkm3 = ((Ny-5)-1)*Nx+ik             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),fp1(ii),dy)
!-------------------------------------------------------------------------------
!      for j = Ny-1
       ijkp1 = ((Ny  )-1)*Nx+ik
       ijk   = ((Ny-1)-1)*Nx+ik
       ijkm1 = ((Ny-2)-1)*Nx+ik
       ijkm2 = ((Ny-3)-1)*Nx+ik
       ijkm3 = ((Ny-4)-1)*Nx+ik             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),fp1(ii),fp2(ii),dy)
!-------------------------------------------------------------------------------
!      for j = Ny
       ijk   = ((Ny  )-1)*Nx+ik
       ijkm1 = ((Ny-1)-1)*Nx+ik
       ijkm2 = ((Ny-2)-1)*Nx+ik
       ijkm3 = ((Ny-3)-1)*Nx+ik             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),fp1(ii),fp2(ii),fp3(ii),dy)
 30  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm)
   end subroutine FCDy5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fifth-order central difference
!    to calculate the first derivative fp along y-direction
!    for periodic boundary condition.
!-------------------------------------------------------------------------------
   subroutine PCDy5th3D_MPI(f,fp,dy,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input data    : dy
!  input arrays  : f(Nxyz)
!  output arrays : fp(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy  = Nx*Ny
     Nxz  = Nx*Nz
     Nxz2 = Nxz*2
     Nxz3 = Nxz*3
     allocate (fm3(Nxz),fm2(Nxz),fm1(Nxz),fp1(Nxz),fp2(Nxz),fp3(Nxz))
     allocate (fsend(Nxz3),frecvp(Nxz3),frecvm(Nxz3))
!-------------------------------------------------------------------------------
!  if Q = 1, it doesn't transfer data in MPI
!
     if (iDend == iDstart) then
       call getBCy3D(f(1),fp3(1),4,Nx,Ny,Nz)
       call getBCy3D(f(1),fp2(1),3,Nx,Ny,Nz)
       call getBCy3D(f(1),fp1(1),2,Nx,Ny,Nz)
!
       call getBCy3D(f(1),fm1(1),Ny-1,Nx,Ny,Nz)
       call getBCy3D(f(1),fm2(1),Ny-2,Nx,Ny,Nz)
       call getBCy3D(f(1),fm3(1),Ny-3,Nx,Ny,Nz)
       goto 123
     endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCy3D(f(1),fsend(1     ),1,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz +1),2,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend (1),Nxz3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nxz3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCy3D(f(1),fsend(1     ),Ny  ,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz +1),Ny-1,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz2+1),Ny-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend (1),Nxz3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nxz3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
     if (iD == iDstart) then
!
       call getBCy3D(f(1),fsend(1     ),2,Nx,Ny,Nz)
       call getBCy3D(f(1),fsend(Nxz +1),3,Nx,Ny,Nz)
       call getBCy3D(f(1),fsend(Nxz2+1),4,Nx,Ny,Nz)
       call MPI_SENDRECV(fsend (1),Nxz3,MPI_REAL8,iDend,120, &
                         frecvm(1),Nxz3,MPI_REAL8,iDend,120, &
                         MPI_WORLD,ISTATUS,IERR)
     endif
!
     if (iD == iDend) then 
!
       call getBCy3D(f(1),fsend(1     ),Ny-1,Nx,Ny,Nz)
       call getBCy3D(f(1),fsend(Nxz +1),Ny-2,Nx,Ny,Nz)
       call getBCy3D(f(1),fsend(Nxz2+1),Ny-3,Nx,Ny,Nz)
       call MPI_SENDRECV(fsend (1),Nxz3,MPI_REAL8,iDstart,120, &
                         frecvp(1),Nxz3,MPI_REAL8,iDstart,120, &
                         MPI_WORLD,ISTATUS,IERR) 
     endif
!-------------------------------------------------------------------------------
     fp1(1:Nxz) = frecvp(1     :Nxz )
     fp2(1:Nxz) = frecvp(Nxz +1:Nxz2)
     fp3(1:Nxz) = frecvp(Nxz2+1:Nxz3)
     fm1(1:Nxz) = frecvm(1     :Nxz )
     fm2(1:Nxz) = frecvm(Nxz +1:Nxz2)
     fm3(1:Nxz) = frecvm(Nxz2+1:Nxz3)
!   
123  continue                                          
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,ik,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 k = 1, Nz
     do 30 i = 1, Nx
       ii = i+(k-1)*Nx
       ik = i+(k-1)*Nxy
!      
       do j = 4, Ny-3
         ijkp3 = ((j+3)-1)*Nx+ik
         ijkp2 = ((j+2)-1)*Nx+ik
         ijkp1 = ((j+1)-1)*Nx+ik
         ijk   = ((j  )-1)*Nx+ik
         ijkm1 = ((j-1)-1)*Nx+ik
         ijkm2 = ((j-2)-1)*Nx+ik
         ijkm3 = ((j-3)-1)*Nx+ik       
!
         fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dy)
       enddo      
!-------------------------------------------------------------------------------
!      for j = 1
       ijkp3 = (4-1)*Nx+ik
       ijkp2 = (3-1)*Nx+ik
       ijkp1 = (2-1)*Nx+ik
       ijk   = (1-1)*Nx+ik  
!
       fp(ijk) = FCD5th(fm3(ii),fm2(ii),fm1(ii),f(ijkp1),f(ijkp2),f(ijkp3),dy)
!-------------------------------------------------------------------------------
!      for j = 2
       ijkp3 = (5-1)*Nx+ik
       ijkp2 = (4-1)*Nx+ik
       ijkp1 = (3-1)*Nx+ik
       ijk   = (2-1)*Nx+ik
       ijkm1 = (1-1)*Nx+ik    
!
       fp(ijk) = FCD5th(fm2(ii),fm1(ii),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dy)
!-------------------------------------------------------------------------------
!      for j = 3
       ijkp3 = (6-1)*Nx+ik
       ijkp2 = (5-1)*Nx+ik
       ijkp1 = (4-1)*Nx+ik
       ijk   = (3-1)*Nx+ik
       ijkm1 = (2-1)*Nx+ik
       ijkm2 = (1-1)*Nx+ik    
!
       fp(ijk) = FCD5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dy)
!-------------------------------------------------------------------------------
!      for j = Ny-2
       ijkp2 = ((Ny  )-1)*Nx+ik
       ijkp1 = ((Ny-1)-1)*Nx+ik
       ijk   = ((Ny-2)-1)*Nx+ik
       ijkm1 = ((Ny-3)-1)*Nx+ik
       ijkm2 = ((Ny-4)-1)*Nx+ik
       ijkm3 = ((Ny-5)-1)*Nx+ik             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),fp1(ii),dy)
!-------------------------------------------------------------------------------
!      for j = Ny-1
       ijkp1 = ((Ny  )-1)*Nx+ik
       ijk   = ((Ny-1)-1)*Nx+ik
       ijkm1 = ((Ny-2)-1)*Nx+ik
       ijkm2 = ((Ny-3)-1)*Nx+ik
       ijkm3 = ((Ny-4)-1)*Nx+ik             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),fp1(ii),fp2(ii),dy)
!-------------------------------------------------------------------------------
!      for j = Ny
       ijk   = ((Ny  )-1)*Nx+ik
       ijkm1 = ((Ny-1)-1)*Nx+ik
       ijkm2 = ((Ny-2)-1)*Nx+ik
       ijkm3 = ((Ny-3)-1)*Nx+ik             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),fp1(ii),fp2(ii),fp3(ii),dy)
 30  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm) 
   end subroutine PCDy5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of the fifth-order central difference
!    to calculate the first derivative fp along z-direction
!    for free boundary condition.
!-------------------------------------------------------------------------------
   subroutine FCDz5th3D_MPI(f,fp,dz,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input data    : dz
!  input arrays  : f(Nxyz)
!  output arrays : fp(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy  = Nx*Ny
     Nxy2 = 2*Nxy
     Nxy3 = 3*Nxy
     allocate (fm3(Nxy),fm2(Nxy),fm1(Nxy),fp1(Nxy),fp2(Nxy),fp3(Nxy))
     allocate (fsend(Nxy3),frecvp(Nxy3),frecvm(Nxy3))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCz3D(f(1),fsend(     1),1,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy +1),2,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nxy3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fp1(1:Nxy) = frecvp(1     :Nxy )
     fp2(1:Nxy) = frecvp(Nxy +1:Nxy2)
     fp3(1:Nxy) = frecvp(Nxy2+1:Nxy3)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCz3D(f(1),fsend(     1),Nz  ,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy +1),Nz-1,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy2+1),Nz-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nxy3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fm1(1:Nxy) = frecvm(1     :Nxy )
     fm2(1:Nxy) = frecvm(Nxy +1:Nxy2)
     fm3(1:Nxy) = frecvm(Nxy2+1:Nxy3)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
     if (iD == iDstart) then
       call getBCz3D(f(1),fm1(1),1,Nx,Ny,Nz)      
       call getBCz3D(f(1),fm2(1),1,Nx,Ny,Nz)
       call getBCz3D(f(1),fm3(1),1,Nx,Ny,Nz)      
     endif
!
     if (iD == iDend) then           
       call getBCz3D(f(1),fp1(1),Nz,Nx,Ny,Nz)      
       call getBCz3D(f(1),fp2(1),Nz,Nx,Ny,Nz)
       call getBCz3D(f(1),fp3(1),Nz,Nx,Ny,Nz)      
     endif                                           
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,ij,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 j = 1, Ny
     do 30 i = 1, Nx
       ii = i+(j-1)*Nx
       ij = ii !i+(j-1)*Nx
!      
       do k = 4, Nz-3
         ijkp3 = ((k+3)-1)*Nxy+ij
         ijkp2 = ((k+2)-1)*Nxy+ij
         ijkp1 = ((k+1)-1)*Nxy+ij
         ijk   = ((k  )-1)*Nxy+ij
         ijkm1 = ((k-1)-1)*Nxy+ij
         ijkm2 = ((k-2)-1)*Nxy+ij
         ijkm3 = ((k-3)-1)*Nxy+ij       
!
         fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dz)
       enddo      
!-------------------------------------------------------------------------------
!      for k = 1
       ijkp3 = (4-1)*Nxy+ij
       ijkp2 = (3-1)*Nxy+ij
       ijkp1 = (2-1)*Nxy+ij
       ijk   = (1-1)*Nxy+ij  
!
       fp(ijk) = FCD5th(fm3(ii),fm2(ii),fm1(ii),f(ijkp1),f(ijkp2),f(ijkp3),dz)
!-------------------------------------------------------------------------------
!      for k = 2
       ijkp3 = (5-1)*Nxy+ij
       ijkp2 = (4-1)*Nxy+ij
       ijkp1 = (3-1)*Nxy+ij
       ijk   = (2-1)*Nxy+ij
       ijkm1 = (1-1)*Nxy+ij    
!
       fp(ijk) = FCD5th(fm2(ii),fm1(ii),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dz)
!-------------------------------------------------------------------------------
!      for k = 3
       ijkp3 = (6-1)*Nxy+ij
       ijkp2 = (5-1)*Nxy+ij
       ijkp1 = (4-1)*Nxy+ij
       ijk   = (3-1)*Nxy+ij
       ijkm1 = (2-1)*Nxy+ij 
       ijkm2 = (1-1)*Nxy+ij   
!
       fp(ijk) = FCD5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dz)
!-------------------------------------------------------------------------------
!      for k = Nz-2
       ijkp2 = ((Nz  )-1)*Nxy+ij
       ijkp1 = ((Nz-1)-1)*Nxy+ij
       ijk   = ((Nz-2)-1)*Nxy+ij
       ijkm1 = ((Nz-3)-1)*Nxy+ij
       ijkm2 = ((Nz-4)-1)*Nxy+ij
       ijkm3 = ((Nz-5)-1)*Nxy+ij             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),fp1(ii),dz)
!-------------------------------------------------------------------------------
!      for k = Nz-1
       ijkp1 = ((Nz  )-1)*Nxy+ij
       ijk   = ((Nz-1)-1)*Nxy+ij
       ijkm1 = ((Nz-2)-1)*Nxy+ij
       ijkm2 = ((Nz-3)-1)*Nxy+ij
       ijkm3 = ((Nz-4)-1)*Nxy+ij             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),fp1(ii),fp2(ii),dz)
!-------------------------------------------------------------------------------
!      for k = Nz
       ijk   = ((Nz  )-1)*Nxy+ij
       ijkm1 = ((Nz-1)-1)*Nxy+ij
       ijkm2 = ((Nz-2)-1)*Nxy+ij
       ijkm3 = ((Nz-3)-1)*Nxy+ij             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),fp1(ii),fp2(ii),fp3(ii),dz) 
 30  continue
!$OMP end parallel do     
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm) 
   end subroutine FCDz5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of the fifth-order central difference
!    to calculate the first derivative fp along z-direction
!    for periodic boundary condition.
!-------------------------------------------------------------------------------
   subroutine PCDz5th3D_MPI(f,fp,dz,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input data    : dz
!  input arrays  : f(Nxyz)
!  output arrays : fp(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm 
     Nxy  = Nx*Ny
     Nxy2 = 2*Nxy
     Nxy3 = 3*Nxy
     allocate (fm3(Nxy),fm2(Nxy),fm1(Nxy),fp1(Nxy),fp2(Nxy),fp3(Nxy))
     allocate (fsend(Nxy3),frecvp(Nxy3),frecvm(Nxy3))
!-------------------------------------------------------------------------------
!  if R = 1, it doesn't transfer data in MPI
!
     if (iDend == iDstart) then
       call getBCz3D(f(1),fp1(1),2,Nx,Ny,Nz)
       call getBCz3D(f(1),fp2(1),3,Nx,Ny,Nz)
       call getBCz3D(f(1),fp3(1),4,Nx,Ny,Nz)
!       
       call getBCz3D(f(1),fm1(1),Nz-1,Nx,Ny,Nz)
       call getBCz3D(f(1),fm2(1),Nz-2,Nx,Ny,Nz)
       call getBCz3D(f(1),fm3(1),Nz-3,Nx,Ny,Nz)
       goto 123
     endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCz3D(f(1),fsend(     1),1,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy +1),2,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nxy3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCz3D(f(1),fsend(     1),Nz  ,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy +1),Nz-1,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy2+1),Nz-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nxy3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
     if (iD == iDstart) then
!
       call getBCz3D(f(1),fsend(     1),2,Nx,Ny,Nz)
       call getBCz3D(f(1),fsend(Nxy +1),3,Nx,Ny,Nz)
       call getBCz3D(f(1),fsend(Nxy2+1),4,Nx,Ny,Nz)
       call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDend,120, &
                         frecvm(1),Nxy3,MPI_REAL8,iDend,120, &
                         MPI_WORLD,ISTATUS,IERR)
     endif
!
     if (iD == iDend) then      
!
       call getBCz3D(f(1),fsend(     1),Nz-1,Nx,Ny,Nz)
       call getBCz3D(f(1),fsend(Nxy +1),Nz-2,Nx,Ny,Nz)
       call getBCz3D(f(1),fsend(Nxy2+1),Nz-3,Nx,Ny,Nz)
       call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDstart,120, &
                         frecvp(1),Nxy3,MPI_REAL8,iDstart,120, &
                         MPI_WORLD,ISTATUS,IERR) 
     endif
!-------------------------------------------------------------------------------
     fp1(1:Nxy) = frecvp(1     :Nxy )
     fp2(1:Nxy) = frecvp(Nxy +1:Nxy2)
     fp3(1:Nxy) = frecvp(Nxy2+1:Nxy3)
     fm1(1:Nxy) = frecvm(1     :Nxy )
     fm2(1:Nxy) = frecvm(Nxy +1:Nxy2)
     fm3(1:Nxy) = frecvm(Nxy2+1:Nxy3)
!    
123  continue                                          
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,ij,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 j = 1, Ny
     do 30 i = 1, Nx
       ii = i+(j-1)*Nx
       ij = ii !i+(j-1)*Nx
!      
       do k = 3, Nz-3
         ijkp3 = ((k+3)-1)*Nxy+ij
         ijkp2 = ((k+2)-1)*Nxy+ij
         ijkp1 = ((k+1)-1)*Nxy+ij
         ijk   = ((k  )-1)*Nxy+ij
         ijkm1 = ((k-1)-1)*Nxy+ij
         ijkm2 = ((k-2)-1)*Nxy+ij
         ijkm3 = ((k-3)-1)*Nxy+ij       
!
         fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dz)
       enddo
!-------------------------------------------------------------------------------
!      for k = 1
       ijkp3 = (4-1)*Nxy+ij
       ijkp2 = (3-1)*Nxy+ij
       ijkp1 = (2-1)*Nxy+ij
       ijk   = (1-1)*Nxy+ij  
!
       fp(ijk) = FCD5th(fm3(ii),fm2(ii),fm1(ii),f(ijkp1),f(ijkp2),f(ijkp3),dz)
!-------------------------------------------------------------------------------
!      for k = 2
       ijkp3 = (5-1)*Nxy+ij
       ijkp2 = (4-1)*Nxy+ij
       ijkp1 = (3-1)*Nxy+ij
       ijk   = (2-1)*Nxy+ij
       ijkm1 = (1-1)*Nxy+ij    
!
       fp(ijk) = FCD5th(fm2(ii),fm1(ii),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dz)
!-------------------------------------------------------------------------------
!      for k = 3
       ijkp3 = (6-1)*Nxy+ij
       ijkp2 = (5-1)*Nxy+ij
       ijkp1 = (4-1)*Nxy+ij
       ijk   = (3-1)*Nxy+ij
       ijkm1 = (2-1)*Nxy+ij 
       ijkm2 = (1-1)*Nxy+ij   
!
       fp(ijk) = FCD5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dz)
!-------------------------------------------------------------------------------
!      for k = Nz-2
       ijkp2 = ((Nz  )-1)*Nxy+ij
       ijkp1 = ((Nz-1)-1)*Nxy+ij
       ijk   = ((Nz-2)-1)*Nxy+ij
       ijkm1 = ((Nz-3)-1)*Nxy+ij
       ijkm2 = ((Nz-4)-1)*Nxy+ij
       ijkm3 = ((Nz-5)-1)*Nxy+ij             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),fp1(ii),dz)
!-------------------------------------------------------------------------------
!      for k = Nz-1
       ijkp1 = ((Nz  )-1)*Nxy+ij
       ijk   = ((Nz-1)-1)*Nxy+ij
       ijkm1 = ((Nz-2)-1)*Nxy+ij
       ijkm2 = ((Nz-3)-1)*Nxy+ij
       ijkm3 = ((Nz-4)-1)*Nxy+ij             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),fp1(ii),fp2(ii),dz)
!-------------------------------------------------------------------------------
!      for k = Nz
       ijk   = ((Nz  )-1)*Nxy+ij
       ijkm1 = ((Nz-1)-1)*Nxy+ij
       ijkm2 = ((Nz-2)-1)*Nxy+ij
       ijkm3 = ((Nz-3)-1)*Nxy+ij             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),fp1(ii),fp2(ii),fp3(ii),dz) 
 30  continue
!$OMP end parallel do     
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm)  
   end subroutine PCDz5th3D_MPI
!===============================================================================
!===============================================================================
   function FCD5th(fm3,fm2,fm1,fp1,fp2,fp3,dx)
     implicit double precision (A-H,O-Z)
     a160 = 1.d0/60.d0
     a320 = 3.d0/20.d0
     a34  = 0.75d0
!
     sumA = sumABC(a160*fp3,-a320*fp2,a34*fp1)
     sumB = sumABC(-a34*fm1,a320*fm2,-a160*fm3)
     FCD5th = DIFAB(sumA,-sumB)/dx
   end function FCD5th
!===============================================================================
   function FCD3rd(fm2,fm1,fp1,fp2,dx)
     implicit double precision (A-H,O-Z)
     a112 = 1.d0/12.d0
     a23  = 2.d0/3.d0
!
     sumA = DIFAB(a23*fp1,a112*fp2)
     sumB = DIFAB(a112*fm2,a23*fm1)
     FCD3rd = DIFAB(sumA,-sumB)/dx
   end function FCD3rd
!===============================================================================
!===============================================================================
   function FCD1st(fm1,fp1,dx)
     implicit double precision (A-H,O-Z)
     FCD1st = DIFAB(fp1,fm1)/dx
   end function FCD1st
!===============================================================================
!===============================================================================
!  This is the step of fifth-order central difference
!    to calculate the second derivative fpp along x-direction
!    for free boundary condition.
!-------------------------------------------------------------------------------
   subroutine FCDxx5th3D_MPI(f,fpp,dx,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fpp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy = Nx*Ny
     Nyz = Ny*Nz
     Nyz2 = 2*Nyz
     Nyz3 = 3*Nyz
     dx2 = dx*dx
     allocate (fm3(Nyz),fm2(Nyz),fm1(Nyz),fp1(Nyz),fp2(Nyz),fp3(Nyz))
     allocate (fsend(Nyz3),frecvp(Nyz3),frecvm(Nyz3))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCx3D(f(1),fsend(     1),1,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz +1),2,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nyz3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nyz3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fp1(1:Nyz) = frecvp(1     :Nyz )
     fp2(1:Nyz) = frecvp(Nyz +1:Nyz2)
     fp3(1:Nyz) = frecvp(Nyz2+1:Nyz3)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCx3D(f(1),fsend(     1),Nx  ,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz +1),Nx-1,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz2+1),Nx-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nyz3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nyz3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fm1(1:Nyz) = frecvm(1     :Nyz )
     fm2(1:Nyz) = frecvm(Nyz +1:Nyz2)
     fm3(1:Nyz) = frecvm(Nyz2+1:Nyz3)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
     if (iD == iDstart) then
       call getBCx3D(f(1),fm1(1),1,Nx,Ny,Nz)      
       call getBCx3D(f(1),fm2(1),1,Nx,Ny,Nz)
       call getBCx3D(f(1),fm3(1),1,Nx,Ny,Nz)      
     endif
!
     if (iD == iDend) then           
       call getBCx3D(f(1),fp1(1),Nx,Nx,Ny,Nz)      
       call getBCx3D(f(1),fp2(1),Nx,Nx,Ny,Nz)
       call getBCx3D(f(1),fp3(1),Nx,Nx,Ny,Nz)      
     endif                                         
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,jk,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 k = 1, Nz
     do 30 j = 1, Ny
       ii = j+(k-1)*Ny
       jk = (j-1)*Nx+(k-1)*Nxy
!      
       do i = 4, Nx-3
         ijkp3 = (i+3)+jk
         ijkp2 = (i+2)+jk
         ijkp1 = (i+1)+jk
         ijk   = (i  )+jk
         ijkm1 = (i-1)+jk
         ijkm2 = (i-2)+jk
         ijkm3 = (i-3)+jk
         fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                           f(ijkp1),f(ijkp2),f(ijkp3),dx2)
       enddo      
!-------------------------------------------------------------------------------
!     for i = 1
       ijkp3 = 4+jk
       ijkp2 = 3+jk
       ijkp1 = 2+jk
       ijk   = 1+jk
!
       fpp(ijk) = FSD5th(fm3(ii),fm2(ii),fm1(ii),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dx2)
!-------------------------------------------------------------------------------
!      for i = 2
       ijkp3 = 5+jk
       ijkp2 = 4+jk
       ijkp1 = 3+jk
       ijk   = 2+jk
       ijkm1 = 1+jk     
!
       fpp(ijk) = FSD5th(fm2(ii),fm1(ii),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dx2)
!-------------------------------------------------------------------------------
!      for i = 3
       ijkp3 = 6+jk
       ijkp2 = 5+jk
       ijkp1 = 4+jk
       ijk   = 3+jk
       ijkm1 = 2+jk
       ijkm2 = 1+jk     
!
       fpp(ijk) = FSD5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dx2)
!-------------------------------------------------------------------------------
!      for i = Nx-2
       ijkp2 = (Nx  )+jk
       ijkp1 = (Nx-1)+jk
       ijk   = (Nx-2)+jk
       ijkm1 = (Nx-3)+jk
       ijkm2 = (Nx-4)+jk
       ijkm3 = (Nx-5)+jk              
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),fp1(ii),dx2)
!-------------------------------------------------------------------------------
!      for i = Nx-1
       ijkp1 = (Nx  )+jk
       ijk   = (Nx-1)+jk
       ijkm1 = (Nx-2)+jk
       ijkm2 = (Nx-3)+jk
       ijkm3 = (Nx-4)+jk              
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),fp1(ii),fp2(ii),dx2)
!-------------------------------------------------------------------------------
!      for i = Nx
       ijk   = (Nx  )+jk
       ijkm1 = (Nx-1)+jk
       ijkm2 = (Nx-2)+jk
       ijkm3 = (Nx-3)+jk              
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         fp1(ii),fp2(ii),fp3(ii),dx2)
 30  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend)  
   end subroutine FCDxx5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fifth-order central difference
!    to calculate the second derivative fpp along x-direction
!    for periodic boundary condition.
!-------------------------------------------------------------------------------
   subroutine PCDxx5th3D_MPI(f,fpp,dx,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input data    : dx
!  input arrays  : f(Nxyz)
!  output arrays : fpp(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fpp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy = Nx*Ny
     Nyz = Ny*Nz
     Nyz2 = 2*Nyz
     Nyz3 = 3*Nyz
     dx2 = dx*dx
     allocate (fm3(Nyz),fm2(Nyz),fm1(Nyz),fp1(Nyz),fp2(Nyz),fp3(Nyz))
     allocate (fsend(Nyz3),frecvp(Nyz3),frecvm(Nyz3))
!-------------------------------------------------------------------------------
!  if P = 1, it doesn't transfer data in MPI
!
     if (iDend == iDstart) then       
       call getBCx3D(f(1),fp3(1),4,Nx,Ny,Nz)
       call getBCx3D(f(1),fp2(1),3,Nx,Ny,Nz)
       call getBCx3D(f(1),fp1(1),2,Nx,Ny,Nz)
!       
       call getBCx3D(f(1),fm1(1),Nx-1,Nx,Ny,Nz)
       call getBCx3D(f(1),fm2(1),Nx-2,Nx,Ny,Nz)
       call getBCx3D(f(1),fm3(1),Nx-3,Nx,Ny,Nz)
       goto 123
     endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCx3D(f(1),fsend(1     ),1,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz +1),2,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz2+1),3,Nx,Ny,Nz)
!
     call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nyz3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCx3D(f(1),fsend(1     ),Nx  ,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz +1),Nx-1,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz2+1),Nx-2,Nx,Ny,Nz)
!
     call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nyz3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
     if (iD == iDstart) then
!
       call getBCx3D(f(1),fsend(1     ),2,Nx,Ny,Nz)
       call getBCx3D(f(1),fsend(Nyz +1),3,Nx,Ny,Nz)
       call getBCx3D(f(1),fsend(Nyz2+1),4,Nx,Ny,Nz)
!
       call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDend,120, &
                         frecvm(1),Nyz3,MPI_REAL8,iDend,120, &
                         MPI_WORLD,ISTATUS,IERR)       
     endif
!
     if (iD == iDend) then      
!
       call getBCx3D(f(1),fsend(1     ),Nx-1,Nx,Ny,Nz)
       call getBCx3D(f(1),fsend(Nyz +1),Nx-2,Nx,Ny,Nz)
       call getBCx3D(f(1),fsend(Nyz2+1),Nx-3,Nx,Ny,Nz)
!
       call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDstart,120, &
                         frecvp(1),Nyz3,MPI_REAL8,iDstart,120, &
                         MPI_WORLD,ISTATUS,IERR)
     endif
!-------------------------------------------------------------------------------
     fp1(1:Nyz) = frecvp(1     :Nyz )
     fp2(1:Nyz) = frecvp(Nyz +1:Nyz2)
     fp3(1:Nyz) = frecvp(Nyz2+1:Nyz3)
     fm1(1:Nyz) = frecvm(1     :Nyz )
     fm2(1:Nyz) = frecvm(Nyz +1:Nyz2)
     fm3(1:Nyz) = frecvm(Nyz2+1:Nyz3)
!    
123  continue                                          
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,jk,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 k = 1, Nz
     do 30 j = 1, Ny
       ii = j+(k-1)*Ny
       jk = (j-1)*Nx+(k-1)*Nxy
!      
       do i = 4, Nx-3
         ijkp3 = (i+3)+jk
         ijkp2 = (i+2)+jk
         ijkp1 = (i+1)+jk
         ijk   = (i  )+jk
         ijkm1 = (i-1)+jk
         ijkm2 = (i-2)+jk
         ijkm3 = (i-3)+jk
         fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                           f(ijkp1),f(ijkp2),f(ijkp3),dx2)
       enddo      
!-------------------------------------------------------------------------------
!     for i = 1
       ijkp3 = 4+jk
       ijkp2 = 3+jk
       ijkp1 = 2+jk
       ijk   = 1+jk
!
       fpp(ijk) = FSD5th(fm3(ii),fm2(ii),fm1(ii),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dx2)
!-------------------------------------------------------------------------------
!      for i = 2
       ijkp3 = 5+jk
       ijkp2 = 4+jk
       ijkp1 = 3+jk
       ijk   = 2+jk
       ijkm1 = 1+jk     
!
       fpp(ijk) = FSD5th(fm2(ii),fm1(ii),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dx2)
!-------------------------------------------------------------------------------
!      for i = 3
       ijkp3 = 6+jk
       ijkp2 = 5+jk
       ijkp1 = 4+jk
       ijk   = 3+jk
       ijkm1 = 2+jk
       ijkm2 = 1+jk     
!
       fpp(ijk) = FSD5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dx2)
!-------------------------------------------------------------------------------
!      for i = Nx-2
       ijkp2 = (Nx  )+jk
       ijkp1 = (Nx-1)+jk
       ijk   = (Nx-2)+jk
       ijkm1 = (Nx-3)+jk
       ijkm2 = (Nx-4)+jk
       ijkm3 = (Nx-5)+jk              
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),fp1(ii),dx2)
!-------------------------------------------------------------------------------
!      for i = Nx-1
       ijkp1 = (Nx  )+jk
       ijk   = (Nx-1)+jk
       ijkm1 = (Nx-2)+jk
       ijkm2 = (Nx-3)+jk
       ijkm3 = (Nx-4)+jk              
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),fp1(ii),fp2(ii),dx2)
!-------------------------------------------------------------------------------
!      for i = Nx
       ijk   = (Nx  )+jk
       ijkm1 = (Nx-1)+jk
       ijkm2 = (Nx-2)+jk
       ijkm3 = (Nx-3)+jk              
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         fp1(ii),fp2(ii),fp3(ii),dx2)
 30  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm) 
   end subroutine PCDxx5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fifth-order central difference
!    to calculate the second derivative fpp along y-direction
!    for free boundary condition.
!-------------------------------------------------------------------------------
   subroutine FCDyy5th3D_MPI(f,fpp,dy,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input data    : dy
!  input arrays  : f(Nxyz)
!  output arrays : fpp(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fpp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm 
     Nxy = Nx*Ny
     Nxz = Nx*Nz
     Nxz2 = 2*Nxz
     Nxz3 = 3*Nxz
     dy2 = dy*dy
     allocate (fm3(Nxz),fm2(Nxz),fm1(Nxz),fp1(Nxz),fp2(Nxz),fp3(Nxz))
     allocate (fsend(Nxz3),frecvp(Nxz3),frecvm(Nxz3))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCy3D(f(1),fsend(     1),1,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz +1),2,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxz3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nxz3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)

     fp1(1:Nxz) = frecvp(1     :Nxz )
     fp2(1:Nxz) = frecvp(Nxz +1:Nxz2)
     fp3(1:Nxz) = frecvp(Nxz2+1:Nxz3)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCy3D(f(1),fsend(     1),Ny  ,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz +1),Ny-1,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz2+1),Ny-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxz3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nxz3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fm1(1:Nxz) = frecvm(1     :Nxz )
     fm2(1:Nxz) = frecvm(Nxz +1:Nxz2)
     fm3(1:Nxz) = frecvm(Nxz2+1:Nxz3)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
     if (iD == iDstart) then
       call getBCy3D(f(1),fm1(1),1,Nx,Ny,Nz)      
       call getBCy3D(f(1),fm2(1),1,Nx,Ny,Nz)
       call getBCy3D(f(1),fm3(1),1,Nx,Ny,Nz)     
     endif
!
     if (iD == iDend) then           
       call getBCy3D(f(1),fp1(1),Ny,Nx,Ny,Nz)      
       call getBCy3D(f(1),fp2(1),Ny,Nx,Ny,Nz)
       call getBCy3D(f(1),fp3(1),Ny,Nx,Ny,Nz)      
     endif                                           
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,ik,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 k = 1, Nz
     do 30 i = 1, Nx
       ii = i+(k-1)*Nx
       ik = i+(k-1)*Nxy
!      
       do j = 4, Ny-3
         ijkp3 = ((j+3)-1)*Nx+ik
         ijkp2 = ((j+2)-1)*Nx+ik
         ijkp1 = ((j+1)-1)*Nx+ik
         ijk   = ((j  )-1)*Nx+ik
         ijkm1 = ((j-1)-1)*Nx+ik
         ijkm2 = ((j-2)-1)*Nx+ik
         ijkm3 = ((j-3)-1)*Nx+ik       
!
         fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                           f(ijkp1),f(ijkp2),f(ijkp3),dy2)
       enddo      
!-------------------------------------------------------------------------------
!      for j = 1
       ijkp3 = (4-1)*Nx+ik
       ijkp2 = (3-1)*Nx+ik
       ijkp1 = (2-1)*Nx+ik
       ijk   = (1-1)*Nx+ik  
!
       fpp(ijk) = FSD5th(fm3(ii),fm2(ii),fm1(ii),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dy2)
!-------------------------------------------------------------------------------
!      for j = 2
       ijkp3 = (5-1)*Nx+ik
       ijkp2 = (4-1)*Nx+ik
       ijkp1 = (3-1)*Nx+ik
       ijk   = (2-1)*Nx+ik
       ijkm1 = (1-1)*Nx+ik    
!
       fpp(ijk) = FSD5th(fm2(ii),fm1(ii),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dy2)
!-------------------------------------------------------------------------------
!      for j = 3
       ijkp3 = (6-1)*Nx+ik
       ijkp2 = (5-1)*Nx+ik
       ijkp1 = (4-1)*Nx+ik
       ijk   = (3-1)*Nx+ik
       ijkm1 = (2-1)*Nx+ik
       ijkm2 = (1-1)*Nx+ik    
!
       fpp(ijk) = FSD5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dy2)
!-------------------------------------------------------------------------------
!      for j = Ny-2
       ijkp2 = ((Ny  )-1)*Nx+ik
       ijkp1 = ((Ny-1)-1)*Nx+ik
       ijk   = ((Ny-2)-1)*Nx+ik
       ijkm1 = ((Ny-3)-1)*Nx+ik
       ijkm2 = ((Ny-4)-1)*Nx+ik
       ijkm3 = ((Ny-5)-1)*Nx+ik             
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),fp1(ii),dy2)
!-------------------------------------------------------------------------------
!      for j = Ny-1
       ijkp1 = ((Ny  )-1)*Nx+ik
       ijk   = ((Ny-1)-1)*Nx+ik
       ijkm1 = ((Ny-2)-1)*Nx+ik
       ijkm2 = ((Ny-3)-1)*Nx+ik
       ijkm3 = ((Ny-4)-1)*Nx+ik             
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),fp1(ii),fp2(ii),dy2)
!-------------------------------------------------------------------------------
!      for j = Ny
       ijk   = ((Ny  )-1)*Nx+ik
       ijkm1 = ((Ny-1)-1)*Nx+ik
       ijkm2 = ((Ny-2)-1)*Nx+ik
       ijkm3 = ((Ny-3)-1)*Nx+ik             
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         fp1(ii),fp2(ii),fp3(ii),dy2)
 30  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend)
   end subroutine FCDyy5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fifth-order central difference
!    to calculate the second derivative fpp along y-direction
!    for periodic boundary condition.
!-------------------------------------------------------------------------------
   subroutine PCDyy5th3D_MPI(f,fpp,dy,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input data    : dy
!  input arrays  : f(Nxyz)
!  output arrays : fpp(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fpp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm 
     Nxy = Nx*Ny
     Nxz = Nx*Nz
     Nxz2 = 2*Nxz
     Nxz3 = 3*Nxz
     dy2 = dy*dy
     allocate (fm3(Nxz),fm2(Nxz),fm1(Nxz),fp1(Nxz),fp2(Nxz),fp3(Nxz))
     allocate (fsend(Nxz3),frecvp(Nxz3),frecvm(Nxz3))
!-------------------------------------------------------------------------------
!  if Q = 1, it doesn't transfer data in MPI
!
     if (iDend == iDstart) then
       call getBCy3D(f(1),fp3(1),4,Nx,Ny,Nz)
       call getBCy3D(f(1),fp2(1),3,Nx,Ny,Nz)
       call getBCy3D(f(1),fp1(1),2,Nx,Ny,Nz)
!
       call getBCy3D(f(1),fm1(1),Ny-1,Nx,Ny,Nz)
       call getBCy3D(f(1),fm2(1),Ny-2,Nx,Ny,Nz)
       call getBCy3D(f(1),fm3(1),Ny-3,Nx,Ny,Nz)
       goto 123
     endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCy3D(f(1),fsend(     1),1,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz +1),2,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxz3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nxz3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!

     call getBCy3D(f(1),fsend(     1),Ny  ,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz +1),Ny-1,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz2+1),Ny-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxz3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nxz3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
     if (iD == iDstart) then
!
       call getBCy3D(f(1),fsend(1     ),2,Nx,Ny,Nz)
       call getBCy3D(f(1),fsend(Nxz +1),3,Nx,Ny,Nz)
       call getBCy3D(f(1),fsend(Nxz2+1),4,Nx,Ny,Nz)
       call MPI_SENDRECV(fsend (1),Nxz3,MPI_REAL8,iDend,120, &
                         frecvm(1),Nxz3,MPI_REAL8,iDend,120, &
                         MPI_WORLD,ISTATUS,IERR)
     endif
!
     if (iD == iDend) then 
!
       call getBCy3D(f(1),fsend(1     ),Ny-1,Nx,Ny,Nz)
       call getBCy3D(f(1),fsend(Nxz +1),Ny-2,Nx,Ny,Nz)
       call getBCy3D(f(1),fsend(Nxz2+1),Ny-3,Nx,Ny,Nz)
       call MPI_SENDRECV(fsend (1),Nxz3,MPI_REAL8,iDstart,120, &
                         frecvp(1),Nxz3,MPI_REAL8,iDstart,120, &
                         MPI_WORLD,ISTATUS,IERR) 
     endif
!-------------------------------------------------------------------------------
     fp1(1:Nxz) = frecvp(1     :Nxz )
     fp2(1:Nxz) = frecvp(Nxz +1:Nxz2)
     fp3(1:Nxz) = frecvp(Nxz2+1:Nxz3)
     fm1(1:Nxz) = frecvm(1     :Nxz )
     fm2(1:Nxz) = frecvm(Nxz +1:Nxz2)
     fm3(1:Nxz) = frecvm(Nxz2+1:Nxz3)
!    
123  continue                                          
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,ik,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 k = 1, Nz
     do 30 i = 1, Nx
       ii = i+(k-1)*Nx
       ik = i+(k-1)*Nxy
!      
       do j = 4, Ny-3
         ijkp3 = ((j+3)-1)*Nx+ik
         ijkp2 = ((j+2)-1)*Nx+ik
         ijkp1 = ((j+1)-1)*Nx+ik
         ijk   = ((j  )-1)*Nx+ik
         ijkm1 = ((j-1)-1)*Nx+ik
         ijkm2 = ((j-2)-1)*Nx+ik
         ijkm3 = ((j-3)-1)*Nx+ik       
!
         fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                           f(ijkp1),f(ijkp2),f(ijkp3),dy2)
       enddo      
!-------------------------------------------------------------------------------
!      for j = 1
       ijkp3 = (4-1)*Nx+ik
       ijkp2 = (3-1)*Nx+ik
       ijkp1 = (2-1)*Nx+ik
       ijk   = (1-1)*Nx+ik  
!
       fpp(ijk) = FSD5th(fm3(ii),fm2(ii),fm1(ii),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dy2)
!-------------------------------------------------------------------------------
!      for j = 2
       ijkp3 = (5-1)*Nx+ik
       ijkp2 = (4-1)*Nx+ik
       ijkp1 = (3-1)*Nx+ik
       ijk   = (2-1)*Nx+ik
       ijkm1 = (1-1)*Nx+ik    
!
       fpp(ijk) = FSD5th(fm2(ii),fm1(ii),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dy2)
!-------------------------------------------------------------------------------
!      for j = 3
       ijkp3 = (6-1)*Nx+ik
       ijkp2 = (5-1)*Nx+ik
       ijkp1 = (4-1)*Nx+ik
       ijk   = (3-1)*Nx+ik
       ijkm1 = (2-1)*Nx+ik
       ijkm2 = (1-1)*Nx+ik    
!
       fpp(ijk) = FSD5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dy2)
!-------------------------------------------------------------------------------
!      for j = Ny-2
       ijkp2 = ((Ny  )-1)*Nx+ik
       ijkp1 = ((Ny-1)-1)*Nx+ik
       ijk   = ((Ny-2)-1)*Nx+ik
       ijkm1 = ((Ny-3)-1)*Nx+ik
       ijkm2 = ((Ny-4)-1)*Nx+ik
       ijkm3 = ((Ny-5)-1)*Nx+ik             
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),fp1(ii),dy2)
!-------------------------------------------------------------------------------
!      for j = Ny-1
       ijkp1 = ((Ny  )-1)*Nx+ik
       ijk   = ((Ny-1)-1)*Nx+ik
       ijkm1 = ((Ny-2)-1)*Nx+ik
       ijkm2 = ((Ny-3)-1)*Nx+ik
       ijkm3 = ((Ny-4)-1)*Nx+ik             
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),fp1(ii),fp2(ii),dy2)
!-------------------------------------------------------------------------------
!      for j = Ny
       ijk   = ((Ny  )-1)*Nx+ik
       ijkm1 = ((Ny-1)-1)*Nx+ik
       ijkm2 = ((Ny-2)-1)*Nx+ik
       ijkm3 = ((Ny-3)-1)*Nx+ik             
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         fp1(ii),fp2(ii),fp3(ii),dy2)
 30  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm) 
   end subroutine PCDyy5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of the fifth-order central difference
!    to calculate the second derivative fpp along z-direction
!    for free boundary condition.
!-------------------------------------------------------------------------------
   subroutine FCDzz5th3D_MPI(f,fpp,dz,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input data    : dz
!  input arrays  : f(Nxyz)
!  output arrays : fpp(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fpp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy = Nx*Ny
     Nxy2 = 2*Nxy
     Nxy3 = 2*Nxy
     dz2 = dz*dz
     allocate (fm3(Nxy),fm2(Nxy),fm1(Nxy),fp1(Nxy),fp2(Nxy),fp3(Nxy))
     allocate (fsend(Nxy3),frecvp(Nxy3),frecvm(Nxy3))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCz3D(f(1),fsend(     1),1,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy +1),2,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nxy3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fp1(1:Nxy) = frecvp(1     :Nxy )
     fp2(1:Nxy) = frecvp(Nxy +1:Nxy2)
     fp3(1:Nxy) = frecvp(Nxy2+1:Nxy3)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCz3D(f(1),fsend(     1),Nz  ,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy +1),Nz-1,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy2+1),Nz-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nxy3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fm1(1:Nxy) = frecvm(1     :Nxy )
     fm2(1:Nxy) = frecvm(Nxy +1:Nxy2)
     fm3(1:Nxy) = frecvm(Nxy2+1:Nxy3)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
     if (iD == iDstart) then
       call getBCz3D(f(1),fm1(1),1,Nx,Ny,Nz)      
       call getBCz3D(f(1),fm2(1),1,Nx,Ny,Nz)
       call getBCz3D(f(1),fm3(1),1,Nx,Ny,Nz)      
     endif
!
     if (iD == iDend) then           
       call getBCz3D(f(1),fp1(1),Nz,Nx,Ny,Nz)      
       call getBCz3D(f(1),fp2(1),Nz,Nx,Ny,Nz)
       call getBCz3D(f(1),fp3(1),Nz,Nx,Ny,Nz)      
     endif                                           
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,ij,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 j = 1, Ny
     do 30 i = 1, Nx
       ii = i+(j-1)*Nx
       ij = ii !i+(j-1)*Nx
!      
       do k = 4, Nz-3
         ijkp3 = ((k+3)-1)*Nxy+ij
         ijkp2 = ((k+2)-1)*Nxy+ij
         ijkp1 = ((k+1)-1)*Nxy+ij
         ijk   = ((k  )-1)*Nxy+ij
         ijkm1 = ((k-1)-1)*Nxy+ij
         ijkm2 = ((k-2)-1)*Nxy+ij
         ijkm3 = ((k-3)-1)*Nxy+ij       
!
         fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                           f(ijkp1),f(ijkp2),f(ijkp3),dz2)
       enddo      
!-------------------------------------------------------------------------------
!      for k = 1
       ijkp3 = (4-1)*Nxy+ij
       ijkp2 = (3-1)*Nxy+ij
       ijkp1 = (2-1)*Nxy+ij
       ijk   = (1-1)*Nxy+ij  
!
       fpp(ijk) = FSD5th(fm3(ii),fm2(ii),fm1(ii),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dz2)
!-------------------------------------------------------------------------------
!      for k = 2
       ijkp3 = (5-1)*Nxy+ij
       ijkp2 = (4-1)*Nxy+ij
       ijkp1 = (3-1)*Nxy+ij
       ijk   = (2-1)*Nxy+ij
       ijkm1 = (1-1)*Nxy+ij    
!
       fpp(ijk) = FSD5th(fm2(ii),fm1(ii),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dz2)
!-------------------------------------------------------------------------------
!      for k = 3
       ijkp3 = (6-1)*Nxy+ij
       ijkp2 = (5-1)*Nxy+ij
       ijkp1 = (4-1)*Nxy+ij
       ijk   = (3-1)*Nxy+ij
       ijkm1 = (2-1)*Nxy+ij 
       ijkm2 = (1-1)*Nxy+ij   
!
       fpp(ijk) = FSD5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dz2)
!-------------------------------------------------------------------------------
!      for k = Nz-2
       ijkp2 = ((Nz  )-1)*Nxy+ij
       ijkp1 = ((Nz-1)-1)*Nxy+ij
       ijk   = ((Nz-2)-1)*Nxy+ij
       ijkm1 = ((Nz-3)-1)*Nxy+ij
       ijkm2 = ((Nz-4)-1)*Nxy+ij
       ijkm3 = ((Nz-5)-1)*Nxy+ij             
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),fp1(ii),dz2)
!-------------------------------------------------------------------------------
!      for k = Nz-1
       ijkp1 = ((Nz  )-1)*Nxy+ij
       ijk   = ((Nz-1)-1)*Nxy+ij
       ijkm1 = ((Nz-2)-1)*Nxy+ij
       ijkm2 = ((Nz-3)-1)*Nxy+ij
       ijkm3 = ((Nz-4)-1)*Nxy+ij             
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),fp1(ii),fp2(ii),dz2)
!-------------------------------------------------------------------------------
!      for k = Nz
       ijk   = ((Nz  )-1)*Nxy+ij
       ijkm1 = ((Nz-1)-1)*Nxy+ij
       ijkm2 = ((Nz-2)-1)*Nxy+ij
       ijkm3 = ((Nz-3)-1)*Nxy+ij             
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         fp1(ii),fp2(ii),fp3(ii),dz2) 
 30  continue
!$OMP end parallel do     
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm) 
   end subroutine FCDzz5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of the fifth-order central difference
!    to calculate the second derivative fpp along z-direction
!    for periodic boundary condition.
!-------------------------------------------------------------------------------
   subroutine PCDzz5th3D_MPI(f,fpp,dz,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input data    : dz
!  input arrays  : f(Nxyz)
!  output arrays : fpp(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fpp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm 
     Nxy = Nx*Ny
     Nxy2 = 2*Nxy
     Nxy3 = 3*Nxy
     dz2 = dz*dz
     allocate (fm3(Nxy),fm2(Nxy),fm1(Nxy),fp1(Nxy),fp2(Nxy),fp3(Nxy))
     allocate (fsend(Nxy3),frecvp(Nxy3),frecvm(Nxy3))
!-------------------------------------------------------------------------------
!  if R = 1, it doesn't transfer data in MPI
!
     if (iDend == iDstart) then
       call getBCz3D(f(1),fp1(1),2,Nx,Ny,Nz)
       call getBCz3D(f(1),fp2(1),3,Nx,Ny,Nz)
       call getBCz3D(f(1),fp3(1),4,Nx,Ny,Nz)
!       
       call getBCz3D(f(1),fm1(1),Nz-1,Nx,Ny,Nz)
       call getBCz3D(f(1),fm2(1),Nz-2,Nx,Ny,Nz)
       call getBCz3D(f(1),fm3(1),Nz-3,Nx,Ny,Nz)
       goto 123
     endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCz3D(f(1),fsend(     1),1,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy +1),2,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nxy3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCz3D(f(1),fsend(     1),Nz  ,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy +1),Nz-1,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy2+1),Nz-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nxy3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
     if (iD == iDstart) then
!
       call getBCz3D(f(1),fsend(     1),2,Nx,Ny,Nz)
       call getBCz3D(f(1),fsend(Nxy +1),3,Nx,Ny,Nz)
       call getBCz3D(f(1),fsend(Nxy2+1),4,Nx,Ny,Nz)
       call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDend,120, &
                         frecvm(1),Nxy3,MPI_REAL8,iDend,120, &
                         MPI_WORLD,ISTATUS,IERR)
     endif
!
     if (iD == iDend) then      
!
       call getBCz3D(f(1),fsend(     1),Nz-1,Nx,Ny,Nz)
       call getBCz3D(f(1),fsend(Nxy +1),Nz-2,Nx,Ny,Nz)
       call getBCz3D(f(1),fsend(Nxy2+1),Nz-3,Nx,Ny,Nz)
       call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDstart,120, &
                         frecvp(1),Nxy3,MPI_REAL8,iDstart,120, &
                         MPI_WORLD,ISTATUS,IERR)
     endif
!-------------------------------------------------------------------------------
     fp1(1:Nxy) = frecvp(1     :Nxy )
     fp2(1:Nxy) = frecvp(Nxy +1:Nxy2)
     fp3(1:Nxy) = frecvp(Nxy2+1:Nxy3)
     fm1(1:Nxy) = frecvm(1     :Nxy )
     fm2(1:Nxy) = frecvm(Nxy +1:Nxy2)
     fm3(1:Nxy) = frecvm(Nxy2+1:Nxy3)
!    
123  continue                                          
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,ij,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 j = 1, Ny
     do 30 i = 1, Nx
       ii = i+(j-1)*Nx
       ij = ii !i+(j-1)*Nx
!      
       do k = 3, Nz-3
         ijkp3 = ((k+3)-1)*Nxy+ij
         ijkp2 = ((k+2)-1)*Nxy+ij
         ijkp1 = ((k+1)-1)*Nxy+ij
         ijk   = ((k  )-1)*Nxy+ij
         ijkm1 = ((k-1)-1)*Nxy+ij
         ijkm2 = ((k-2)-1)*Nxy+ij
         ijkm3 = ((k-3)-1)*Nxy+ij       
!
         fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                           f(ijkp1),f(ijkp2),f(ijkp3),dz2)
       enddo
!-------------------------------------------------------------------------------
!      for k = 1
       ijkp3 = (4-1)*Nxy+ij
       ijkp2 = (3-1)*Nxy+ij
       ijkp1 = (2-1)*Nxy+ij
       ijk   = (1-1)*Nxy+ij  
!
       fpp(ijk) = FSD5th(fm3(ii),fm2(ii),fm1(ii),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dz2)
!-------------------------------------------------------------------------------
!      for k = 2
       ijkp3 = (5-1)*Nxy+ij
       ijkp2 = (4-1)*Nxy+ij
       ijkp1 = (3-1)*Nxy+ij
       ijk   = (2-1)*Nxy+ij
       ijkm1 = (1-1)*Nxy+ij    
!
       fpp(ijk) = FSD5th(fm2(ii),fm1(ii),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dz2)
!-------------------------------------------------------------------------------
!      for k = 3
       ijkp3 = (6-1)*Nxy+ij
       ijkp2 = (5-1)*Nxy+ij
       ijkp1 = (4-1)*Nxy+ij
       ijk   = (3-1)*Nxy+ij
       ijkm1 = (2-1)*Nxy+ij 
       ijkm2 = (1-1)*Nxy+ij   
!
       fpp(ijk) = FSD5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),f(ijkp3),dz2)
!-------------------------------------------------------------------------------
!      for k = Nz-2
       ijkp2 = ((Nz  )-1)*Nxy+ij
       ijkp1 = ((Nz-1)-1)*Nxy+ij
       ijk   = ((Nz-2)-1)*Nxy+ij
       ijkm1 = ((Nz-3)-1)*Nxy+ij
       ijkm2 = ((Nz-4)-1)*Nxy+ij
       ijkm3 = ((Nz-5)-1)*Nxy+ij             
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),f(ijkp2),fp1(ii),dz2)
!-------------------------------------------------------------------------------
!      for k = Nz-1
       ijkp1 = ((Nz  )-1)*Nxy+ij
       ijk   = ((Nz-1)-1)*Nxy+ij
       ijkm1 = ((Nz-2)-1)*Nxy+ij
       ijkm2 = ((Nz-3)-1)*Nxy+ij
       ijkm3 = ((Nz-4)-1)*Nxy+ij             
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         f(ijkp1),fp1(ii),fp2(ii),dz2)
!-------------------------------------------------------------------------------
!      for k = Nz
       ijk   = ((Nz  )-1)*Nxy+ij
       ijkm1 = ((Nz-1)-1)*Nxy+ij
       ijkm2 = ((Nz-2)-1)*Nxy+ij
       ijkm3 = ((Nz-3)-1)*Nxy+ij             
!
       fpp(ijk) = FSD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk), &
                         fp1(ii),fp2(ii),fp3(ii),dz2) 
 30  continue
!$OMP end parallel do     
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm)  
   end subroutine PCDzz5th3D_MPI
!===============================================================================
!===============================================================================
   function FSD5th(fm3,fm2,fm1,f,fp1,fp2,fp3,dx2)
!
!  5th-order central finite difference method for the second order of derivative
!
     implicit double precision (A-H,O-Z)
     a190 = 1.d0/90.d0
     a320 = 3.d0/20.d0
     a32  = 1.5d0
     ac   = 49.d0/18.d0
!
     sumA = sumABC(+a190*fp3,-a320*fp2,+a32 *fp1)
     sumB = sumABC(+a32 *fm1,-a320*fm2,+a190*fm3)
     temp = -ac*f
     FSD5th = sumABC(sumA,sumB,temp)/dx2
   end function FSD5th
!===============================================================================
!===============================================================================
!  This is the step of fifth-order central difference
!    to calculate the first derivative fp along x-direction
!    for free boundary condition.
!    "No MPI version"
!-------------------------------------------------------------------------------
!  input data    : dx
!  input arrays  : f(Nxyz)
!  output arrays : fp(Nxyz)
!-------------------------------------------------------------------------------
   subroutine FCDx5th3D(f,fp,dx,Nx,Ny,Nz)
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3
     Nxy = Nx*Ny
     Nyz = Ny*Nz
     Nyz2 = 2*Nyz
     Nyz3 = 3*Nyz
     allocate (fm3(Nyz),fm2(Nyz),fm1(Nyz),fp1(Nyz),fp2(Nyz),fp3(Nyz))
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
     call getBCx3D(f(1),fm1(1),1,Nx,Ny,Nz)      
     call getBCx3D(f(1),fm2(1),1,Nx,Ny,Nz)
     call getBCx3D(f(1),fm3(1),1,Nx,Ny,Nz)
!          
     call getBCx3D(f(1),fp1(1),Nx,Nx,Ny,Nz)      
     call getBCx3D(f(1),fp2(1),Nx,Nx,Ny,Nz)
     call getBCx3D(f(1),fp3(1),Nx,Nx,Ny,Nz)                                         
!-------------------------------------------------------------------------------
     do 30 k = 1, Nz
     do 30 j = 1, Ny
       ii = j+(k-1)*Ny
       jk = (j-1)*Nx+(k-1)*Nxy
!      
       do i = 4, Nx-3
         ijkp3 = (i+3)+jk
         ijkp2 = (i+2)+jk
         ijkp1 = (i+1)+jk
         ijk   = (i  )+jk
         ijkm1 = (i-1)+jk
         ijkm2 = (i-2)+jk
         ijkm3 = (i-3)+jk
         fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dx)
       enddo      
!-------------------------------------------------------------------------------
!     for i = 1
       ijkp3 = 4+jk
       ijkp2 = 3+jk
       ijkp1 = 2+jk
       ijk   = 1+jk
!
       fp(ijk) = FCD5th(fm3(ii),fm2(ii),fm1(ii),f(ijkp1),f(ijkp2),f(ijkp3),dx)
!-------------------------------------------------------------------------------
!      for i = 2
       ijkp3 = 5+jk
       ijkp2 = 4+jk
       ijkp1 = 3+jk
       ijk   = 2+jk
       ijkm1 = 1+jk     
!
       fp(ijk) = FCD5th(fm2(ii),fm1(ii),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dx)
!-------------------------------------------------------------------------------
!      for i = 3
       ijkp3 = 6+jk
       ijkp2 = 5+jk
       ijkp1 = 4+jk
       ijk   = 3+jk
       ijkm1 = 2+jk
       ijkm2 = 1+jk     
!
       fp(ijk) = FCD5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dx)
!-------------------------------------------------------------------------------
!      for i = Nx-2
       ijkp2 = (Nx  )+jk
       ijkp1 = (Nx-1)+jk
       ijk   = (Nx-2)+jk
       ijkm1 = (Nx-3)+jk
       ijkm2 = (Nx-4)+jk
       ijkm3 = (Nx-5)+jk              
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),fp1(ii),dx)
!-------------------------------------------------------------------------------
!      for i = Nx-1
       ijkp1 = (Nx  )+jk
       ijk   = (Nx-1)+jk
       ijkm1 = (Nx-2)+jk
       ijkm2 = (Nx-3)+jk
       ijkm3 = (Nx-4)+jk              
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),fp1(ii),fp2(ii),dx)
!-------------------------------------------------------------------------------
!      for i = Nx
       ijk   = (Nx  )+jk
       ijkm1 = (Nx-1)+jk
       ijkm2 = (Nx-2)+jk
       ijkm3 = (Nx-3)+jk              
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),fp1(ii),fp2(ii),fp3(ii),dx)
 30  continue
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3)  
   end subroutine FCDx5th3D
!===============================================================================
!===============================================================================
!  This is the step of fifth-order central difference
!    to calculate the first derivative fp along y-direction
!    for free boundary condition.
!    "No MPI version"
!-------------------------------------------------------------------------------
   subroutine FCDy5th3D(f,fp,dy,Nx,Ny,Nz)
!-------------------------------------------------------------------------------
!  input data    : dy
!  input arrays  : f(Nxyz)
!  output arrays : fp(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3
     Nxy  = Nx*Ny
     Nxz  = Nx*Nz
!     Nxz2 = 2*Nxz
!     Nxz3 = 3*Nxz
     allocate (fm3(Nxz),fm2(Nxz),fm1(Nxz),fp1(Nxz),fp2(Nxz),fp3(Nxz))
!-------------------------------------------------------------------------------
!  set the data for free boundary condition

     call getBCy3D(f(1),fm1(1),1,Nx,Ny,Nz)      
     call getBCy3D(f(1),fm2(1),1,Nx,Ny,Nz)
     call getBCy3D(f(1),fm3(1),1,Nx,Ny,Nz)
!           
     call getBCy3D(f(1),fp1(1),Ny,Nx,Ny,Nz)      
     call getBCy3D(f(1),fp2(1),Ny,Nx,Ny,Nz)
     call getBCy3D(f(1),fp3(1),Ny,Nx,Ny,Nz)                                            
!-------------------------------------------------------------------------------
     do 30 k = 1, Nz
     do 30 i = 1, Nx
       ii = i+(k-1)*Nx
       ik = i+(k-1)*Nxy
!      
       do j = 4, Ny-3
         ijkp3 = ((j+3)-1)*Nx+ik
         ijkp2 = ((j+2)-1)*Nx+ik
         ijkp1 = ((j+1)-1)*Nx+ik
         ijk   = ((j  )-1)*Nx+ik
         ijkm1 = ((j-1)-1)*Nx+ik
         ijkm2 = ((j-2)-1)*Nx+ik
         ijkm3 = ((j-3)-1)*Nx+ik       
!
         fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dy)
       enddo      
!-------------------------------------------------------------------------------
!      for j = 1
       ijkp3 = (4-1)*Nx+ik
       ijkp2 = (3-1)*Nx+ik
       ijkp1 = (2-1)*Nx+ik
       ijk   = (1-1)*Nx+ik  
!
       fp(ijk) = FCD5th(fm3(ii),fm2(ii),fm1(ii),f(ijkp1),f(ijkp2),f(ijkp3),dy)
!-------------------------------------------------------------------------------
!      for j = 2
       ijkp3 = (5-1)*Nx+ik
       ijkp2 = (4-1)*Nx+ik
       ijkp1 = (3-1)*Nx+ik
       ijk   = (2-1)*Nx+ik
       ijkm1 = (1-1)*Nx+ik    
!
       fp(ijk) = FCD5th(fm2(ii),fm1(ii),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dy)
!-------------------------------------------------------------------------------
!      for j = 3
       ijkp3 = (6-1)*Nx+ik
       ijkp2 = (5-1)*Nx+ik
       ijkp1 = (4-1)*Nx+ik
       ijk   = (3-1)*Nx+ik
       ijkm1 = (2-1)*Nx+ik
       ijkm2 = (1-1)*Nx+ik    
!
       fp(ijk) = FCD5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dy)
!-------------------------------------------------------------------------------
!      for j = Ny-2
       ijkp2 = ((Ny  )-1)*Nx+ik
       ijkp1 = ((Ny-1)-1)*Nx+ik
       ijk   = ((Ny-2)-1)*Nx+ik
       ijkm1 = ((Ny-3)-1)*Nx+ik
       ijkm2 = ((Ny-4)-1)*Nx+ik
       ijkm3 = ((Ny-5)-1)*Nx+ik             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),fp1(ii),dy)
!-------------------------------------------------------------------------------
!      for j = Ny-1
       ijkp1 = ((Ny  )-1)*Nx+ik
       ijk   = ((Ny-1)-1)*Nx+ik
       ijkm1 = ((Ny-2)-1)*Nx+ik
       ijkm2 = ((Ny-3)-1)*Nx+ik
       ijkm3 = ((Ny-4)-1)*Nx+ik             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),fp1(ii),fp2(ii),dy)
!-------------------------------------------------------------------------------
!      for j = Ny
       ijk   = ((Ny  )-1)*Nx+ik
       ijkm1 = ((Ny-1)-1)*Nx+ik
       ijkm2 = ((Ny-2)-1)*Nx+ik
       ijkm3 = ((Ny-3)-1)*Nx+ik             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),fp1(ii),fp2(ii),fp3(ii),dy)
 30  continue
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3)
   end subroutine FCDy5th3D
!===============================================================================
!===============================================================================
!  This is the step of the fifth-order central difference
!    to calculate the first derivative fp along z-direction
!    for free boundary condition.
!    "No MPI version"
!-------------------------------------------------------------------------------
   subroutine FCDz5th3D(f,fp,dz,Nx,Ny,Nz)
!-------------------------------------------------------------------------------
!  input data    : dz
!  input arrays  : f(Nxyz)
!  output arrays : fp(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3
     Nxy  = Nx*Ny
     allocate (fm3(Nxy),fm2(Nxy),fm1(Nxy),fp1(Nxy),fp2(Nxy),fp3(Nxy))
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
     call getBCz3D(f(1),fm1(1),1,Nx,Ny,Nz)      
     call getBCz3D(f(1),fm2(1),1,Nx,Ny,Nz)
     call getBCz3D(f(1),fm3(1),1,Nx,Ny,Nz)
!           
     call getBCz3D(f(1),fp1(1),Nz,Nx,Ny,Nz)      
     call getBCz3D(f(1),fp2(1),Nz,Nx,Ny,Nz)
     call getBCz3D(f(1),fp3(1),Nz,Nx,Ny,Nz)                                      
!-------------------------------------------------------------------------------
     do 30 j = 1, Ny
     do 30 i = 1, Nx
       ii = i+(j-1)*Nx
       ij = ii !i+(j-1)*Nx
!      
       do k = 4, Nz-3
         ijkp3 = ((k+3)-1)*Nxy+ij
         ijkp2 = ((k+2)-1)*Nxy+ij
         ijkp1 = ((k+1)-1)*Nxy+ij
         ijk   = ((k  )-1)*Nxy+ij
         ijkm1 = ((k-1)-1)*Nxy+ij
         ijkm2 = ((k-2)-1)*Nxy+ij
         ijkm3 = ((k-3)-1)*Nxy+ij       
!
         fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dz)
       enddo      
!-------------------------------------------------------------------------------
!      for k = 1
       ijkp3 = (4-1)*Nxy+ij
       ijkp2 = (3-1)*Nxy+ij
       ijkp1 = (2-1)*Nxy+ij
       ijk   = (1-1)*Nxy+ij  
!
       fp(ijk) = FCD5th(fm3(ii),fm2(ii),fm1(ii),f(ijkp1),f(ijkp2),f(ijkp3),dz)
!-------------------------------------------------------------------------------
!      for k = 2
       ijkp3 = (5-1)*Nxy+ij
       ijkp2 = (4-1)*Nxy+ij
       ijkp1 = (3-1)*Nxy+ij
       ijk   = (2-1)*Nxy+ij
       ijkm1 = (1-1)*Nxy+ij    
!
       fp(ijk) = FCD5th(fm2(ii),fm1(ii),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dz)
!-------------------------------------------------------------------------------
!      for k = 3
       ijkp3 = (6-1)*Nxy+ij
       ijkp2 = (5-1)*Nxy+ij
       ijkp1 = (4-1)*Nxy+ij
       ijk   = (3-1)*Nxy+ij
       ijkm1 = (2-1)*Nxy+ij 
       ijkm2 = (1-1)*Nxy+ij   
!
       fp(ijk) = FCD5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),f(ijkp3),dz)
!-------------------------------------------------------------------------------
!      for k = Nz-2
       ijkp2 = ((Nz  )-1)*Nxy+ij
       ijkp1 = ((Nz-1)-1)*Nxy+ij
       ijk   = ((Nz-2)-1)*Nxy+ij
       ijkm1 = ((Nz-3)-1)*Nxy+ij
       ijkm2 = ((Nz-4)-1)*Nxy+ij
       ijkm3 = ((Nz-5)-1)*Nxy+ij             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2),fp1(ii),dz)
!-------------------------------------------------------------------------------
!      for k = Nz-1
       ijkp1 = ((Nz  )-1)*Nxy+ij
       ijk   = ((Nz-1)-1)*Nxy+ij
       ijkm1 = ((Nz-2)-1)*Nxy+ij
       ijkm2 = ((Nz-3)-1)*Nxy+ij
       ijkm3 = ((Nz-4)-1)*Nxy+ij             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),fp1(ii),fp2(ii),dz)
!-------------------------------------------------------------------------------
!      for k = Nz
       ijk   = ((Nz  )-1)*Nxy+ij
       ijkm1 = ((Nz-1)-1)*Nxy+ij
       ijkm2 = ((Nz-2)-1)*Nxy+ij
       ijkm3 = ((Nz-3)-1)*Nxy+ij             
!
       fp(ijk) = FCD5th(f(ijkm3),f(ijkm2),f(ijkm1),fp1(ii),fp2(ii),fp3(ii),dz) 
 30  continue     
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3) 
   end subroutine FCDz5th3D
!===============================================================================
!===============================================================================
!  This is the step of fifth-order central difference
!    to calculate the first derivative fp and second derivative fpp
!    along x-direction for free boundary condition.
!-------------------------------------------------------------------------------
   subroutine FCDx5th6D_MPI(f,fp,fpp,dx,Nx,Ny,Nz,Nux,Nuy,Nuz,iDm1,iD,iDp1,iDstart,iDend)
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1), fpp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy = Nx*Ny
     Nyz = Ny*Nz
     Nuxy  = Nux*Nuy
     Nuxyz = Nux*Nuy*Nuz
     NN  = Nyz*Nuxyz
     NN2 = NN*2
     NN3 = NN*3
     dx2 = dx*dx
     allocate (fm3(NN),fm2(NN),fm1(NN),fp1(NN),fp2(NN),fp3(NN))
     allocate (fsend(NN3),frecvp(NN3),frecvm(NN3))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCx6D(f(1),fsend(    1),1,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCx6D(f(1),fsend(NN +1),2,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCx6D(f(1),fsend(NN2+1),3,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call MPI_SENDRECV(fsend (1),NN3,MPI_REAL8,iDm1,100, &
                       frecvp(1),NN3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fp1(1:NN) = frecvp(1    :NN )
     fp2(1:NN) = frecvp(NN +1:NN2)
     fp3(1:NN) = frecvp(NN2+1:NN3)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCx6D(f(1),fsend(    1),Nx  ,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCx6D(f(1),fsend(NN +1),Nx-1,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCx6D(f(1),fsend(NN2+1),Nx-2,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call MPI_SENDRECV(fsend (1),NN3,MPI_REAL8,iDp1,110, &
                       frecvm(1),NN3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fm1(1:NN) = frecvm(1    :NN )
     fm2(1:NN) = frecvm(NN +1:NN2)
     fm3(1:NN) = frecvm(NN2+1:NN3)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
     if (iD == iDstart) then
       call getBCx6D(f(1),fm1(1),1,Nx,Ny,Nz,Nux,Nuy,Nuz)      
       call getBCx6D(f(1),fm2(1),1,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCx6D(f(1),fm3(1),1,Nx,Ny,Nz,Nux,Nuy,Nuz)      
     endif
!
     if (iD == iDend) then           
       call getBCx6D(f(1),fp1(1),Nx,Nx,Ny,Nz,Nux,Nuy,Nuz)      
       call getBCx6D(f(1),fp2(1),Nx,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCx6D(f(1),fp3(1),Nx,Nx,Ny,Nz,Nux,Nuy,Nuz)      
     endif                                         
!-------------------------------------------------------------------------------
!$OMP parallel do private(iim3,iim2,iim1,ii,iip1,iip2,iip3,ijku,jk1,jk2,ii1) collapse(3) 
     do 35 ku = 1, Nuz
     do 35 ju = 1, Nuy
     do 35 iu = 1, Nux
       ijku = iu+(ju-1)*Nux+(ku-1)*Nuxy
!
       do 30 k = 1, Nz
       do 30 j = 1, Ny
         jk1 = j+(k-1)*Ny
         jk2 = (j-1)*Nx+(k-1)*Nxy
         ii1 = ijku+(jk1-1)*Nuxyz
!      
         do i = 4, Nx-3
           iip3 = ijku+((i+3)+jk2-1)*Nuxyz
           iip2 = ijku+((i+2)+jk2-1)*Nuxyz
           iip1 = ijku+((i+1)+jk2-1)*Nuxyz
           ii   = ijku+((i  )+jk2-1)*Nuxyz
           iim1 = ijku+((i-1)+jk2-1)*Nuxyz
           iim2 = ijku+((i-2)+jk2-1)*Nuxyz
           iim3 = ijku+((i-3)+jk2-1)*Nuxyz
           fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),f(iip2),f(iip3),dx)
           fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
         enddo      
!-------------------------------------------------------------------------------
!  for i = 1
         iip3 = ijku+(4+jk2-1)*Nuxyz
         iip2 = ijku+(3+jk2-1)*Nuxyz
         iip1 = ijku+(2+jk2-1)*Nuxyz
         ii   = ijku+(1+jk2-1)*Nuxyz
!
         fp(ii)  = FCD5th(fm3(ii1),fm2(ii1),fm1(ii1),f(iip1),f(iip2),f(iip3),dx)
         fpp(ii) = FSD5th(fm3(ii1),fm2(ii1),fm1(ii1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = 2
         iip3 = ijku+(5+jk2-1)*Nuxyz
         iip2 = ijku+(4+jk2-1)*Nuxyz
         iip1 = ijku+(3+jk2-1)*Nuxyz
         ii   = ijku+(2+jk2-1)*Nuxyz
         iim1 = ijku+(1+jk2-1)*Nuxyz     
!
         fp(ii)  = FCD5th(fm2(ii1),fm1(ii1),f(iim1),f(iip1),f(iip2),f(iip3),dx)
         fpp(ii) = FSD5th(fm2(ii1),fm1(ii1),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = 3
         iip3 = ijku+(6+jk2-1)*Nuxyz
         iip2 = ijku+(5+jk2-1)*Nuxyz
         iip1 = ijku+(4+jk2-1)*Nuxyz
         ii   = ijku+(3+jk2-1)*Nuxyz
         iim1 = ijku+(2+jk2-1)*Nuxyz
         iim2 = ijku+(1+jk2-1)*Nuxyz     
!
         fp(ii)  = FCD5th(fm1(ii1),f(iim2),f(iim1),f(iip1),f(iip2),f(iip3),dx)
         fpp(ii) = FSD5th(fm1(ii1),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx-2
         iip2 = ijku+((Nx  )+jk2-1)*Nuxyz
         iip1 = ijku+((Nx-1)+jk2-1)*Nuxyz
         ii   = ijku+((Nx-2)+jk2-1)*Nuxyz
         iim1 = ijku+((Nx-3)+jk2-1)*Nuxyz
         iim2 = ijku+((Nx-4)+jk2-1)*Nuxyz
         iim3 = ijku+((Nx-5)+jk2-1)*Nuxyz              
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),f(iip2),fp1(ii1),dx)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),fp1(ii1),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx-1
         iip1 = ijku+((Nx  )+jk2-1)*Nuxyz
         ii   = ijku+((Nx-1)+jk2-1)*Nuxyz
         iim1 = ijku+((Nx-2)+jk2-1)*Nuxyz
         iim2 = ijku+((Nx-3)+jk2-1)*Nuxyz
         iim3 = ijku+((Nx-4)+jk2-1)*Nuxyz              
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),fp1(ii1),fp2(ii1),dx)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),fp1(ii1),fp2(ii1),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx
         ii   = ijku+((Nx  )+jk2-1)*Nuxyz
         iim1 = ijku+((Nx-1)+jk2-1)*Nuxyz
         iim2 = ijku+((Nx-2)+jk2-1)*Nuxyz
         iim3 = ijku+((Nx-3)+jk2-1)*Nuxyz              
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),fp1(ii1),fp2(ii1),fp3(ii1),dx)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),fp1(ii1),fp2(ii1),fp3(ii1),dx2)
 30    continue
!
 35  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm)  
   end subroutine FCDx5th6D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fifth-order central difference
!    to calculate the first derivative fp and second derivative fpp
!    along x-direction for free boundary condition.
!-------------------------------------------------------------------------------
   subroutine PCDx5th6D_MPI(f,fp,fpp,dx,Nx,Ny,Nz,Nux,Nuy,Nuz,iDm1,iD,iDp1,iDstart,iDend)
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1), fpp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy = Nx*Ny
     Nyz = Ny*Nz
     Nuxy  = Nux*Nuy
     Nuxyz = Nux*Nuy*Nuz
     NN  = Nyz*Nuxyz
     NN2 = NN*2
     NN3 = NN*3
     dx2 = dx*dx
     allocate (fm3(NN),fm2(NN),fm1(NN),fp1(NN),fp2(NN),fp3(NN))
     allocate (fsend(NN3),frecvp(NN3),frecvm(NN3))
!-------------------------------------------------------------------------------
!  if P = 1, it doesn't transfer data in MPI
!
     if (iDend == iDstart) then       
       call getBCx6D(f(1),fp3(1),4,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCx6D(f(1),fp2(1),3,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCx6D(f(1),fp1(1),2,Nx,Ny,Nz,Nux,Nuy,Nuz)
!       
       call getBCx6D(f(1),fm1(1),Nx-1,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCx6D(f(1),fm2(1),Nx-2,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCx6D(f(1),fm3(1),Nx-3,Nx,Ny,Nz,Nux,Nuy,Nuz)
!
       do 110 ku = 1, Nuz
       do 110 ju = 1, Nuy
       do 110 iu = 1, Nux
!
           ijku = iu+(ju-1)*Nux+(ku-1)*Nuxy
!
           ijk1 = 1
           ii1  = ijku+(ijk1-1)*Nuxyz
!
           ijkn = Nx
           iin  = ijku+(ijkn-1)*Nuxyz
!
           f(iin) = f(ii1)
 110   continue
!
       goto 123
     endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCx6D(f(1),fsend(1    ),1,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCx6D(f(1),fsend(NN +1),2,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCx6D(f(1),fsend(NN2+1),3,Nx,Ny,Nz,Nux,Nuy,Nuz)
!
     call MPI_SENDRECV(fsend (1),NN3,MPI_REAL8,iDm1,100, &
                       frecvp(1),NN3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR) 
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCx6D(f(1),fsend(1    ),Nx  ,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCx6D(f(1),fsend(NN +1),Nx-1,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCx6D(f(1),fsend(NN2+1),Nx-2,Nx,Ny,Nz,Nux,Nuy,Nuz)
!
     call MPI_SENDRECV(fsend (1),NN3,MPI_REAL8,iDp1,110, &
                       frecvm(1),NN3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR) 
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
     if (iD == iDstart) then
!
       call getBCx6D(f(1),fsend(1    ),2,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCx6D(f(1),fsend(NN +1),3,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCx6D(f(1),fsend(NN2+1),4,Nx,Ny,Nz,Nux,Nuy,Nuz)
!
       call MPI_SENDRECV(fsend (1),NN3,MPI_REAL8,iDend,120, &
                         frecvm(1),NN3,MPI_REAL8,iDend,120, &
                         MPI_WORLD,ISTATUS,IERR)      
     endif
!
     if (iD == iDend) then      
!
       call getBCx6D(f(1),fsend(1    ),Nx-1,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCx6D(f(1),fsend(NN +1),Nx-2,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCx6D(f(1),fsend(NN2+1),Nx-3,Nx,Ny,Nz,Nux,Nuy,Nuz)
!
       call MPI_SENDRECV(fsend (1),NN3,MPI_REAL8,iDstart,120, &
                         frecvp(1),NN3,MPI_REAL8,iDstart,120, &
                         MPI_WORLD,ISTATUS,IERR)
     endif
!
     fp1(1:NN) = frecvp(1    :NN )
     fp2(1:NN) = frecvp(NN +1:NN2)
     fp3(1:NN) = frecvp(NN2+1:NN3)
     fm1(1:NN) = frecvm(1    :NN )
     fm2(1:NN) = frecvm(NN +1:NN2)
     fm3(1:NN) = frecvm(NN2+1:NN3)
!    
123  continue                                         
!-------------------------------------------------------------------------------
!$OMP parallel do private(iim3,iim2,iim1,ii,iip1,iip2,iip3,ijku,jk1,jk2,ii1) collapse(3) 
     do 35 ku = 1, Nuz
     do 35 ju = 1, Nuy
     do 35 iu = 1, Nux
       ijku = iu+(ju-1)*Nux+(ku-1)*Nuxy
!
       do 30 k = 1, Nz
       do 30 j = 1, Ny
         jk1 = j+(k-1)*Ny
         jk2 = (j-1)*Nx+(k-1)*Nxy
         ii1 = ijku+(jk1-1)*Nuxyz
!      
         do i = 4, Nx-3
           iip3 = ijku+((i+3)+jk2-1)*Nuxyz
           iip2 = ijku+((i+2)+jk2-1)*Nuxyz
           iip1 = ijku+((i+1)+jk2-1)*Nuxyz
           ii   = ijku+((i  )+jk2-1)*Nuxyz
           iim1 = ijku+((i-1)+jk2-1)*Nuxyz
           iim2 = ijku+((i-2)+jk2-1)*Nuxyz
           iim3 = ijku+((i-3)+jk2-1)*Nuxyz
           fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),f(iip2),f(iip3),dx)
           fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
         enddo      
!-------------------------------------------------------------------------------
!  for i = 1
         iip3 = ijku+(4+jk2-1)*Nuxyz
         iip2 = ijku+(3+jk2-1)*Nuxyz
         iip1 = ijku+(2+jk2-1)*Nuxyz
         ii   = ijku+(1+jk2-1)*Nuxyz
!
         fp(ii)  = FCD5th(fm3(ii1),fm2(ii1),fm1(ii1),f(iip1),f(iip2),f(iip3),dx)
         fpp(ii) = FSD5th(fm3(ii1),fm2(ii1),fm1(ii1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = 2
         iip3 = ijku+(5+jk2-1)*Nuxyz
         iip2 = ijku+(4+jk2-1)*Nuxyz
         iip1 = ijku+(3+jk2-1)*Nuxyz
         ii   = ijku+(2+jk2-1)*Nuxyz
         iim1 = ijku+(1+jk2-1)*Nuxyz     
!
         fp(ii)  = FCD5th(fm2(ii1),fm1(ii1),f(iim1),f(iip1),f(iip2),f(iip3),dx)
         fpp(ii) = FSD5th(fm2(ii1),fm1(ii1),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = 3
         iip3 = ijku+(6+jk2-1)*Nuxyz
         iip2 = ijku+(5+jk2-1)*Nuxyz
         iip1 = ijku+(4+jk2-1)*Nuxyz
         ii   = ijku+(3+jk2-1)*Nuxyz
         iim1 = ijku+(2+jk2-1)*Nuxyz
         iim2 = ijku+(1+jk2-1)*Nuxyz     
!
         fp(ii)  = FCD5th(fm1(ii1),f(iim2),f(iim1),f(iip1),f(iip2),f(iip3),dx)
         fpp(ii) = FSD5th(fm1(ii1),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx-2
         iip2 = ijku+((Nx  )+jk2-1)*Nuxyz
         iip1 = ijku+((Nx-1)+jk2-1)*Nuxyz
         ii   = ijku+((Nx-2)+jk2-1)*Nuxyz
         iim1 = ijku+((Nx-3)+jk2-1)*Nuxyz
         iim2 = ijku+((Nx-4)+jk2-1)*Nuxyz
         iim3 = ijku+((Nx-5)+jk2-1)*Nuxyz              
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),f(iip2),fp1(ii1),dx)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),fp1(ii1),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx-1
         iip1 = ijku+((Nx  )+jk2-1)*Nuxyz
         ii   = ijku+((Nx-1)+jk2-1)*Nuxyz
         iim1 = ijku+((Nx-2)+jk2-1)*Nuxyz
         iim2 = ijku+((Nx-3)+jk2-1)*Nuxyz
         iim3 = ijku+((Nx-4)+jk2-1)*Nuxyz              
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),fp1(ii1),fp2(ii1),dx)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),fp1(ii1),fp2(ii1),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx
         ii   = ijku+((Nx  )+jk2-1)*Nuxyz
         iim1 = ijku+((Nx-1)+jk2-1)*Nuxyz
         iim2 = ijku+((Nx-2)+jk2-1)*Nuxyz
         iim3 = ijku+((Nx-3)+jk2-1)*Nuxyz              
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),fp1(ii1),fp2(ii1),fp3(ii1),dx)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),fp1(ii1),fp2(ii1),fp3(ii1),dx2)
 30    continue
!
 35  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm)  
   end subroutine PCDx5th6D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fifth-order central difference
!    to calculate the first derivative fp along y-direction
!    for free boundary condition.
!-------------------------------------------------------------------------------
   subroutine FCDy5th6D_MPI(f,fp,fpp,dy,Nx,Ny,Nz,Nux,Nuy,Nuz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input data    : dy
!  input arrays  : f(Nxyz*Nuxyz)
!  output arrays : fp(Nxyz*Nuxyz), fpp(Nxyz*Nuxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1), fpp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nuxy  = Nux*Nuy
     Nuxyz = Nux*Nuy*Nuz
     Nxy   = Nx*Ny
     Nxz   = Nx*Nz
     NN    = Nxz*Nuxyz
     NN2   = 2*NN
     NN3   = 3*NN
     dy2   = dy*dy
     allocate (fm3(NN),fm2(NN),fm1(NN),fp1(NN),fp2(NN),fp3(NN))
     allocate (fsend(NN3),frecvp(NN3),frecvm(NN3))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCy6D(f(1),fsend(    1),1,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCy6D(f(1),fsend(NN +1),2,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCy6D(f(1),fsend(NN2+1),3,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call MPI_SENDRECV(fsend(1) ,NN3,MPI_REAL8,iDm1,100, &
                       frecvp(1),NN3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fp1(1:NN) = frecvp(1    :NN )
     fp2(1:NN) = frecvp(NN +1:NN2)
     fp3(1:NN) = frecvp(NN2+1:NN3)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCy6D(f(1),fsend(    1),Ny  ,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCy6D(f(1),fsend(NN +1),Ny-1,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCy6D(f(1),fsend(NN2+1),Ny-2,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call MPI_SENDRECV(fsend(1) ,NN3,MPI_REAL8,iDp1,110, &
                       frecvm(1),NN3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fm1(1:NN) = frecvm(1    :NN )
     fm2(1:NN) = frecvm(NN +1:NN2)
     fm3(1:NN) = frecvm(NN2+1:NN3)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
     if (iD == iDstart) then
       call getBCy6D(f(1),fm1(1),1,Nx,Ny,Nz,Nux,Nuy,Nuz)      
       call getBCy6D(f(1),fm2(1),1,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCy6D(f(1),fm3(1),1,Nx,Ny,Nz,Nux,Nuy,Nuz)     
     endif
!
     if (iD == iDend) then           
       call getBCy6D(f(1),fp1(1),Ny,Nx,Ny,Nz,Nux,Nuy,Nuz)      
       call getBCy6D(f(1),fp2(1),Ny,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCy6D(f(1),fp3(1),Ny,Nx,Ny,Nz,Nux,Nuy,Nuz)      
     endif                                           
!-------------------------------------------------------------------------------
!$OMP parallel do private(iim3,iim2,iim1,ii,iip1,iip2,iip3,ijku,ik1,ik2,ii1) collapse(3)
     do 35 ku = 1, Nuz
     do 35 ju = 1, Nuy
     do 35 iu = 1, Nux
       ijku = iu+(ju-1)*Nux+(ku-1)*Nuxy
!
       do 30 k = 1, Nz
       do 30 i = 1, Nx
         ik1 = i+(k-1)*Nx
         ik2 = i+(k-1)*Nxy
         ii1 = ijku+(ik1-1)*Nuxyz
!      
         do j = 4, Ny-3
           iip3 = ijku+(((j+3)-1)*Nx+ik2-1)*Nuxyz
           iip2 = ijku+(((j+2)-1)*Nx+ik2-1)*Nuxyz
           iip1 = ijku+(((j+1)-1)*Nx+ik2-1)*Nuxyz
           ii   = ijku+(((j  )-1)*Nx+ik2-1)*Nuxyz
           iim1 = ijku+(((j-1)-1)*Nx+ik2-1)*Nuxyz
           iim2 = ijku+(((j-2)-1)*Nx+ik2-1)*Nuxyz
           iim3 = ijku+(((j-3)-1)*Nx+ik2-1)*Nuxyz       
!        
           fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),f(iip2),f(iip3),dy)
           fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dy2)
         enddo      
!-------------------------------------------------------------------------------
!  for j = 1
         iip3 = ijku+((4-1)*Nx+ik2-1)*Nuxyz
         iip2 = ijku+((3-1)*Nx+ik2-1)*Nuxyz
         iip1 = ijku+((2-1)*Nx+ik2-1)*Nuxyz
         ii   = ijku+((1-1)*Nx+ik2-1)*Nuxyz  
!
         fp(ii)  = FCD5th(fm3(ii1),fm2(ii1),fm1(ii1),f(iip1),f(iip2),f(iip3),dy)
         fpp(ii) = FSD5th(fm3(ii1),fm2(ii1),fm1(ii1),f(ii),f(iip1),f(iip2),f(iip3),dy2)
!-------------------------------------------------------------------------------
!  for j = 2
         iip3 = ijku+((5-1)*Nx+ik2-1)*Nuxyz
         iip2 = ijku+((4-1)*Nx+ik2-1)*Nuxyz
         iip1 = ijku+((3-1)*Nx+ik2-1)*Nuxyz
         ii   = ijku+((2-1)*Nx+ik2-1)*Nuxyz
         iim1 = ijku+((1-1)*Nx+ik2-1)*Nuxyz    
!
         fp(ii)  = FCD5th(fm2(ii1),fm1(ii1),f(iim1),f(iip1),f(iip2),f(iip3),dy)
         fpp(ii) = FSD5th(fm2(ii1),fm1(ii1),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dy2)
!-------------------------------------------------------------------------------
!  for j = 3
         iip3 = ijku+((6-1)*Nx+ik2-1)*Nuxyz
         iip2 = ijku+((5-1)*Nx+ik2-1)*Nuxyz
         iip1 = ijku+((4-1)*Nx+ik2-1)*Nuxyz
         ii   = ijku+((3-1)*Nx+ik2-1)*Nuxyz
         iim1 = ijku+((2-1)*Nx+ik2-1)*Nuxyz
         iim2 = ijku+((1-1)*Nx+ik2-1)*Nuxyz    
!
         fp(ii)  = FCD5th(fm1(ii1),f(iim2),f(iim1),f(iip1),f(iip2),f(iip3),dy)
         fpp(ii) = FSD5th(fm1(ii1),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dy2)
!-------------------------------------------------------------------------------
!  for j = Ny-2
         iip2 = ijku+(((Ny  )-1)*Nx+ik2-1)*Nuxyz
         iip1 = ijku+(((Ny-1)-1)*Nx+ik2-1)*Nuxyz
         ii   = ijku+(((Ny-2)-1)*Nx+ik2-1)*Nuxyz
         iim1 = ijku+(((Ny-3)-1)*Nx+ik2-1)*Nuxyz
         iim2 = ijku+(((Ny-4)-1)*Nx+ik2-1)*Nuxyz
         iim3 = ijku+(((Ny-5)-1)*Nx+ik2-1)*Nuxyz             
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),f(iip2),fp1(ii1),dy)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),fp1(ii1),dy2)
!-------------------------------------------------------------------------------
!  for j = Ny-1
         iip1 = ijku+(((Ny  )-1)*Nx+ik2-1)*Nuxyz
         ii   = ijku+(((Ny-1)-1)*Nx+ik2-1)*Nuxyz
         iim1 = ijku+(((Ny-2)-1)*Nx+ik2-1)*Nuxyz
         iim2 = ijku+(((Ny-3)-1)*Nx+ik2-1)*Nuxyz
         iim3 = ijku+(((Ny-4)-1)*Nx+ik2-1)*Nuxyz             
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),fp1(ii1),fp2(ii1),dy)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),fp1(ii1),fp2(ii1),dy2)
!-------------------------------------------------------------------------------
!  for j = Ny
         ii   = ijku+(((Ny  )-1)*Nx+ik2-1)*Nuxyz 
         iim1 = ijku+(((Ny-1)-1)*Nx+ik2-1)*Nuxyz 
         iim2 = ijku+(((Ny-2)-1)*Nx+ik2-1)*Nuxyz 
         iim3 = ijku+(((Ny-3)-1)*Nx+ik2-1)*Nuxyz              
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),fp1(ii1),fp2(ii1),fp3(ii1),dy)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),fp1(ii1),fp2(ii1),fp3(ii1),dy2)
 30    continue
!
 35  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm)
   end subroutine FCDy5th6D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fifth-order central difference
!    to calculate the first derivative fp along y-direction
!    for periodic boundary condition.
!-------------------------------------------------------------------------------
   subroutine PCDy5th6D_MPI(f,fp,fpp,dy,Nx,Ny,Nz,Nux,Nuy,Nuz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input data    : dy
!  input arrays  : f(Nxyz*Nuxyz)
!  output arrays : fp(Nxyz*Nuxyz), fpp(Nxyz*Nuxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1), fpp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nuxy  = Nux*Nuy
     Nuxyz = Nux*Nuy*Nuz
     Nxy   = Nx*Ny
     Nxz   = Nx*Nz
     NN    = Nxz*Nuxyz
     NN2   = 2*NN
     NN3   = 3*NN
     dy2   = dy*dy
     allocate (fm3(NN),fm2(NN),fm1(NN),fp1(NN),fp2(NN),fp3(NN))
     allocate (fsend(NN3),frecvp(NN3),frecvm(NN3))
!-------------------------------------------------------------------------------
!  if Q = 1, it doesn't transfer data in MPI
!
     if (iDend == iDstart) then
       call getBCy6D(f(1),fp3(1),4,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCy6D(f(1),fp2(1),3,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCy6D(f(1),fp1(1),2,Nx,Ny,Nz,Nux,Nuy,Nuz)
!
       call getBCy6D(f(1),fm1(1),Ny-1,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCy6D(f(1),fm2(1),Ny-2,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCy6D(f(1),fm3(1),Ny-3,Nx,Ny,Nz,Nux,Nuy,Nuz)
       goto 123
     endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCy6D(f(1),fsend(    1),1,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCy6D(f(1),fsend(NN +1),2,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCy6D(f(1),fsend(NN2+1),3,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call MPI_SENDRECV(fsend (1),NN3,MPI_REAL8,iDm1,100, &
                       frecvp(1),NN3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCy6D(f(1),fsend(    1),Ny  ,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCy6D(f(1),fsend(NN +1),Ny-1,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCy6D(f(1),fsend(NN2+1),Ny-2,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call MPI_SENDRECV(fsend (1),NN3,MPI_REAL8,iDp1,110, &
                       frecvm(1),NN3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
     if (iD == iDstart) then
!
       call getBCy6D(f(1),fsend(    1),2,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCy6D(f(1),fsend(NN +1),3,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCy6D(f(1),fsend(NN2+1),4,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call MPI_SENDRECV(fsend (1),NN3,MPI_REAL8,iDend,120, &
                         frecvm(1),NN3,MPI_REAL8,iDend,120, &
                         MPI_WORLD,ISTATUS,IERR)
     endif
!
     if (iD == iDend) then 
!
       call getBCy6D(f(1),fsend(    1),Ny-1,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCy6D(f(1),fsend(NN +1),Ny-2,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCy6D(f(1),fsend(NN2+1),Ny-3,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call MPI_SENDRECV(fsend (1),NN3,MPI_REAL8,iDstart,120, &
                         frecvp(1),NN3,MPI_REAL8,iDstart,120, &
                         MPI_WORLD,ISTATUS,IERR) 
     endif
!-------------------------------------------------------------------------------
     fp1(1:NN) = frecvp(1    :NN )
     fp2(1:NN) = frecvp(NN +1:NN2)
     fp3(1:NN) = frecvp(NN2+1:NN3)
     fm1(1:NN) = frecvm(1    :NN )
     fm2(1:NN) = frecvm(NN +1:NN2)
     fm3(1:NN) = frecvm(NN2+1:NN3)
!   
123  continue                                          
!-------------------------------------------------------------------------------
!$OMP parallel do private(iim3,iim2,iim1,ii,iip1,iip2,iip3,ijku,ik1,ik2,ii1) collapse(3)
     do 35 ku = 1, Nuz
     do 35 ju = 1, Nuy
     do 35 iu = 1, Nux
       ijku = iu+(ju-1)*Nux+(ku-1)*Nuxy
!
       do 30 k = 1, Nz
       do 30 i = 1, Nx
         ik1 = i+(k-1)*Nx
         ik2 = i+(k-1)*Nxy
         ii1 = ijku+(ik1-1)*Nuxyz
!      
         do j = 4, Ny-3
           iip3 = ijku+(((j+3)-1)*Nx+ik2-1)*Nuxyz
           iip2 = ijku+(((j+2)-1)*Nx+ik2-1)*Nuxyz
           iip1 = ijku+(((j+1)-1)*Nx+ik2-1)*Nuxyz
           ii   = ijku+(((j  )-1)*Nx+ik2-1)*Nuxyz
           iim1 = ijku+(((j-1)-1)*Nx+ik2-1)*Nuxyz
           iim2 = ijku+(((j-2)-1)*Nx+ik2-1)*Nuxyz
           iim3 = ijku+(((j-3)-1)*Nx+ik2-1)*Nuxyz       
!        
           fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),f(iip2),f(iip3),dy)
           fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dy2)
         enddo      
!-------------------------------------------------------------------------------
!  for j = 1
         iip3 = ijku+((4-1)*Nx+ik2-1)*Nuxyz
         iip2 = ijku+((3-1)*Nx+ik2-1)*Nuxyz
         iip1 = ijku+((2-1)*Nx+ik2-1)*Nuxyz
         ii   = ijku+((1-1)*Nx+ik2-1)*Nuxyz  
!
         fp(ii)  = FCD5th(fm3(ii1),fm2(ii1),fm1(ii1),f(iip1),f(iip2),f(iip3),dy)
         fpp(ii) = FSD5th(fm3(ii1),fm2(ii1),fm1(ii1),f(ii),f(iip1),f(iip2),f(iip3),dy2)
!-------------------------------------------------------------------------------
!  for j = 2
         iip3 = ijku+((5-1)*Nx+ik2-1)*Nuxyz
         iip2 = ijku+((4-1)*Nx+ik2-1)*Nuxyz
         iip1 = ijku+((3-1)*Nx+ik2-1)*Nuxyz
         ii   = ijku+((2-1)*Nx+ik2-1)*Nuxyz
         iim1 = ijku+((1-1)*Nx+ik2-1)*Nuxyz    
!
         fp(ii)  = FCD5th(fm2(ii1),fm1(ii1),f(iim1),f(iip1),f(iip2),f(iip3),dy)
         fpp(ii) = FSD5th(fm2(ii1),fm1(ii1),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dy2)
!-------------------------------------------------------------------------------
!  for j = 3
         iip3 = ijku+((6-1)*Nx+ik2-1)*Nuxyz
         iip2 = ijku+((5-1)*Nx+ik2-1)*Nuxyz
         iip1 = ijku+((4-1)*Nx+ik2-1)*Nuxyz
         ii   = ijku+((3-1)*Nx+ik2-1)*Nuxyz
         iim1 = ijku+((2-1)*Nx+ik2-1)*Nuxyz
         iim2 = ijku+((1-1)*Nx+ik2-1)*Nuxyz    
!
         fp(ii)  = FCD5th(fm1(ii1),f(iim2),f(iim1),f(iip1),f(iip2),f(iip3),dy)
         fpp(ii) = FSD5th(fm1(ii1),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dy2)
!-------------------------------------------------------------------------------
!  for j = Ny-2
         iip2 = ijku+(((Ny  )-1)*Nx+ik2-1)*Nuxyz
         iip1 = ijku+(((Ny-1)-1)*Nx+ik2-1)*Nuxyz
         ii   = ijku+(((Ny-2)-1)*Nx+ik2-1)*Nuxyz
         iim1 = ijku+(((Ny-3)-1)*Nx+ik2-1)*Nuxyz
         iim2 = ijku+(((Ny-4)-1)*Nx+ik2-1)*Nuxyz
         iim3 = ijku+(((Ny-5)-1)*Nx+ik2-1)*Nuxyz             
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),f(iip2),fp1(ii1),dy)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),fp1(ii1),dy2)
!-------------------------------------------------------------------------------
!  for j = Ny-1
         iip1 = ijku+(((Ny  )-1)*Nx+ik2-1)*Nuxyz
         ii   = ijku+(((Ny-1)-1)*Nx+ik2-1)*Nuxyz
         iim1 = ijku+(((Ny-2)-1)*Nx+ik2-1)*Nuxyz
         iim2 = ijku+(((Ny-3)-1)*Nx+ik2-1)*Nuxyz
         iim3 = ijku+(((Ny-4)-1)*Nx+ik2-1)*Nuxyz             
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),fp1(ii1),fp2(ii1),dy)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),fp1(ii1),fp2(ii1),dy2)
!-------------------------------------------------------------------------------
!  for j = Ny
         ii   = ijku+(((Ny  )-1)*Nx+ik2-1)*Nuxyz 
         iim1 = ijku+(((Ny-1)-1)*Nx+ik2-1)*Nuxyz 
         iim2 = ijku+(((Ny-2)-1)*Nx+ik2-1)*Nuxyz 
         iim3 = ijku+(((Ny-3)-1)*Nx+ik2-1)*Nuxyz              
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),fp1(ii1),fp2(ii1),fp3(ii1),dy)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),fp1(ii1),fp2(ii1),fp3(ii1),dy2)
 30    continue
!
 35  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm) 
   end subroutine PCDy5th6D_MPI
!===============================================================================
!===============================================================================
!  This is the step of the fifth-order central difference
!    to calculate the first derivative fp along z-direction
!    for free boundary condition.
!-------------------------------------------------------------------------------
   subroutine FCDz5th6D_MPI(f,fp,fpp,dz,Nx,Ny,Nz,Nux,Nuy,Nuz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input data    : dz
!  input arrays  : f(Nxyz*Nuxyz)
!  output arrays : fp(Nxyz*Nuxyz), fpp(Nxyz*Nuxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1), fpp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nuxy  = Nux*Nuy
     Nuxyz = Nux*Nuy*Nuz
     Nxy   = Nx*Ny
     NN    = Nxy*Nuxyz     
     NN2   = 2*NN
     NN3   = 3*NN
     dz2   = dz*dz
     allocate (fm3(NN),fm2(NN),fm1(NN),fp1(NN),fp2(NN),fp3(NN))
     allocate (fsend(NN3),frecvp(NN3),frecvm(NN3))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCz6D(f(1),fsend(    1),1,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCz6D(f(1),fsend(NN +1),2,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCz6D(f(1),fsend(NN2+1),3,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call MPI_SENDRECV(fsend(1) ,NN3,MPI_REAL8,iDm1,100, &
                       frecvp(1),NN3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fp1(1:NN) = frecvp(1    :NN )
     fp2(1:NN) = frecvp(NN +1:NN2)
     fp3(1:NN) = frecvp(NN2+1:NN3)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCz6D(f(1),fsend(    1),Nz  ,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCz6D(f(1),fsend(NN +1),Nz-1,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCz6D(f(1),fsend(NN2+1),Nz-2,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call MPI_SENDRECV(fsend(1) ,NN3,MPI_REAL8,iDp1,110, &
                       frecvm(1),NN3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fm1(1:NN) = frecvm(1    :NN )
     fm2(1:NN) = frecvm(NN +1:NN2)
     fm3(1:NN) = frecvm(NN2+1:NN3)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
     if (iD == iDstart) then
       call getBCz6D(f(1),fm1(1),1,Nx,Ny,Nz,Nux,Nuy,Nuz)      
       call getBCz6D(f(1),fm2(1),1,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCz6D(f(1),fm3(1),1,Nx,Ny,Nz,Nux,Nuy,Nuz)      
     endif
!
     if (iD == iDend) then           
       call getBCz6D(f(1),fp1(1),Nz,Nx,Ny,Nz,Nux,Nuy,Nuz)      
       call getBCz6D(f(1),fp2(1),Nz,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCz6D(f(1),fp3(1),Nz,Nx,Ny,Nz,Nux,Nuy,Nuz)      
     endif                                           
!-------------------------------------------------------------------------------
!$OMP parallel do private(iim3,iim2,iim1,ii,iip1,iip2,iip3,ijku,ij,ii1) collapse(3)
     do 35 ku = 1, Nuz
     do 35 ju = 1, Nuy
     do 35 iu = 1, Nux
       ijku = iu+(ju-1)*Nux+(ku-1)*Nuxy
!
       do 30 j = 1, Ny
       do 30 i = 1, Nx
         ij = i+(j-1)*Nx
         ii1 = ijku+(ij-1)*Nuxyz
!      
         do k = 4, Nz-3
           iip3 = ijku+(((k+3)-1)*Nxy+ij-1)*Nuxyz
           iip2 = ijku+(((k+2)-1)*Nxy+ij-1)*Nuxyz
           iip1 = ijku+(((k+1)-1)*Nxy+ij-1)*Nuxyz
           ii   = ijku+(((k  )-1)*Nxy+ij-1)*Nuxyz
           iim1 = ijku+(((k-1)-1)*Nxy+ij-1)*Nuxyz
           iim2 = ijku+(((k-2)-1)*Nxy+ij-1)*Nuxyz
           iim3 = ijku+(((k-3)-1)*Nxy+ij-1)*Nuxyz       
!
           fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),f(iip2),f(iip3),dz)
           fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dz2)
         enddo      
!-------------------------------------------------------------------------------
!  for k = 1
         iip3 = ijku+((4-1)*Nxy+ij-1)*Nuxyz
         iip2 = ijku+((3-1)*Nxy+ij-1)*Nuxyz
         iip1 = ijku+((2-1)*Nxy+ij-1)*Nuxyz
         ii   = ijku+((1-1)*Nxy+ij-1)*Nuxyz  
!
         fp(ii)  = FCD5th(fm3(ii1),fm2(ii1),fm1(ii1),f(iip1),f(iip2),f(iip3),dz)
         fpp(ii) = FSD5th(fm3(ii1),fm2(ii1),fm1(ii1),f(ii),f(iip1),f(iip2),f(iip3),dz2)
!-------------------------------------------------------------------------------
!  for k = 2
         iip3 = ijku+((5-1)*Nxy+ij-1)*Nuxyz
         iip2 = ijku+((4-1)*Nxy+ij-1)*Nuxyz
         iip1 = ijku+((3-1)*Nxy+ij-1)*Nuxyz
         ii   = ijku+((2-1)*Nxy+ij-1)*Nuxyz
         iim1 = ijku+((1-1)*Nxy+ij-1)*Nuxyz    
!
         fp(ii)  = FCD5th(fm2(ii1),fm1(ii1),f(iim1),f(iip1),f(iip2),f(iip3),dz)
         fpp(ii) = FSD5th(fm2(ii1),fm1(ii1),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dz2)
!-------------------------------------------------------------------------------
!  for k = 3
         iip3 = ijku+((6-1)*Nxy+ij-1)*Nuxyz
         iip2 = ijku+((5-1)*Nxy+ij-1)*Nuxyz
         iip1 = ijku+((4-1)*Nxy+ij-1)*Nuxyz
         ii   = ijku+((3-1)*Nxy+ij-1)*Nuxyz
         iim1 = ijku+((2-1)*Nxy+ij-1)*Nuxyz 
         iim2 = ijku+((1-1)*Nxy+ij-1)*Nuxyz   
!
         fp(ii)  = FCD5th(fm1(ii1),f(iim2),f(iim1),f(iip1),f(iip2),f(iip3),dz)
         fpp(ii) = FSD5th(fm1(ii1),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dz2)
!-------------------------------------------------------------------------------
!  for k = Nz-2
         iip2 = ijku+(((Nz  )-1)*Nxy+ij-1)*Nuxyz
         iip1 = ijku+(((Nz-1)-1)*Nxy+ij-1)*Nuxyz
         ii   = ijku+(((Nz-2)-1)*Nxy+ij-1)*Nuxyz
         iim1 = ijku+(((Nz-3)-1)*Nxy+ij-1)*Nuxyz
         iim2 = ijku+(((Nz-4)-1)*Nxy+ij-1)*Nuxyz
         iim3 = ijku+(((Nz-5)-1)*Nxy+ij-1)*Nuxyz             
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),f(iip2),fp1(ii1),dz)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),fp1(ii1),dz2)
!-------------------------------------------------------------------------------
!  for k = Nz-1
         iip1 = ijku+(((Nz  )-1)*Nxy+ij-1)*Nuxyz
         ii   = ijku+(((Nz-1)-1)*Nxy+ij-1)*Nuxyz
         iim1 = ijku+(((Nz-2)-1)*Nxy+ij-1)*Nuxyz
         iim2 = ijku+(((Nz-3)-1)*Nxy+ij-1)*Nuxyz
         iim3 = ijku+(((Nz-4)-1)*Nxy+ij-1)*Nuxyz             
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),fp1(ii1),fp2(ii1),dz)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),fp1(ii1),fp2(ii1),dz2)
!-------------------------------------------------------------------------------
!  for k = Nz
         ii   = ijku+(((Nz  )-1)*Nxy+ij-1)*Nuxyz
         iim1 = ijku+(((Nz-1)-1)*Nxy+ij-1)*Nuxyz
         iim2 = ijku+(((Nz-2)-1)*Nxy+ij-1)*Nuxyz
         iim3 = ijku+(((Nz-3)-1)*Nxy+ij-1)*Nuxyz             
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),fp1(ii1),fp2(ii1),fp3(ii1),dz)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),fp1(ii1),fp2(ii1),fp3(ii1),dz2) 
 30    continue
!
 35  continue
!$OMP end parallel do     
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm) 
   end subroutine FCDz5th6D_MPI
!===============================================================================
!===============================================================================
!  This is the step of the fifth-order central difference
!    to calculate the first derivative fp along z-direction
!    for periodic boundary condition.
!-------------------------------------------------------------------------------
   subroutine PCDz5th6D_MPI(f,fp,fpp,dz,Nx,Ny,Nz,Nux,Nuy,Nuz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input data    : dz
!  input arrays  : f(Nxyz*Nuxyz)
!  output arrays : fp(Nxyz*Nuxyz), fpp(Nxyz*Nuxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1), fpp(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm 
     Nuxy  = Nux*Nuy
     Nuxyz = Nux*Nuy*Nuz
     Nxy   = Nx*Ny
     NN    = Nxy*Nuxyz     
     NN2   = 2*NN
     NN3   = 3*NN
     dz2   = dz*dz
     allocate (fm3(NN),fm2(NN),fm1(NN),fp1(NN),fp2(NN),fp3(NN))
     allocate (fsend(NN3),frecvp(NN3),frecvm(NN3))
!-------------------------------------------------------------------------------
!  if R = 1, it doesn't transfer data in MPI
!
     if (iDend == iDstart) then
       call getBCz6D(f(1),fp1(1),2,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCz6D(f(1),fp2(1),3,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCz6D(f(1),fp3(1),4,Nx,Ny,Nz,Nux,Nuy,Nuz)
!       
       call getBCz6D(f(1),fm1(1),Nz-1,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCz6D(f(1),fm2(1),Nz-2,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCz6D(f(1),fm3(1),Nz-3,Nx,Ny,Nz,Nux,Nuy,Nuz)
       goto 123
     endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCz6D(f(1),fsend(    1),1,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCz6D(f(1),fsend(NN +1),2,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCz6D(f(1),fsend(NN2+1),3,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call MPI_SENDRECV(fsend(1) ,NN3,MPI_REAL8,iDm1,100, &
                       frecvp(1),NN3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCz6D(f(1),fsend(    1),Nz  ,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCz6D(f(1),fsend(NN +1),Nz-1,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call getBCz6D(f(1),fsend(NN2+1),Nz-2,Nx,Ny,Nz,Nux,Nuy,Nuz)
     call MPI_SENDRECV(fsend(1) ,NN3,MPI_REAL8,iDp1,110, &
                       frecvm(1),NN3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
     if (iD == iDstart) then
!
       call getBCz6D(f(1),fsend(    1),2,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCz6D(f(1),fsend(NN +1),3,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCz6D(f(1),fsend(NN2+1),4,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call MPI_SENDRECV(fsend(1) ,NN3,MPI_REAL8,iDend,120, &
                         frecvm(1),NN3,MPI_REAL8,iDend,120, &
                         MPI_WORLD,ISTATUS,IERR)
     endif
!
     if (iD == iDend) then      
!
       call getBCz6D(f(1),fsend(    1),Nz-1,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCz6D(f(1),fsend(NN +1),Nz-2,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call getBCz6D(f(1),fsend(NN2+1),Nz-3,Nx,Ny,Nz,Nux,Nuy,Nuz)
       call MPI_SENDRECV(fsend(1) ,NN3,MPI_REAL8,iDstart,120, &
                         frecvp(1),NN3,MPI_REAL8,iDstart,120, &
                         MPI_WORLD,ISTATUS,IERR) 
     endif
!-------------------------------------------------------------------------------
     fp1(1:NN) = frecvp(1    :NN )
     fp2(1:NN) = frecvp(NN +1:NN2)
     fp3(1:NN) = frecvp(NN2+1:NN3)
     fm1(1:NN) = frecvm(1    :NN )
     fm2(1:NN) = frecvm(NN +1:NN2)
     fm3(1:NN) = frecvm(NN2+1:NN3)
!    
123  continue                                          
!-------------------------------------------------------------------------------
!$OMP parallel do private(iim3,iim2,iim1,ii,iip1,iip2,iip3,ijku,ij,ii1) collapse(3)
     do 35 ku = 1, Nuz
     do 35 ju = 1, Nuy
     do 35 iu = 1, Nux
       ijku = iu+(ju-1)*Nux+(ku-1)*Nuxy
!
       do 30 j = 1, Ny
       do 30 i = 1, Nx
         ij = i+(j-1)*Nx
         ii1 = ijku+(ij-1)*Nuxyz
!      
         do k = 4, Nz-3
           iip3 = ijku+(((k+3)-1)*Nxy+ij-1)*Nuxyz
           iip2 = ijku+(((k+2)-1)*Nxy+ij-1)*Nuxyz
           iip1 = ijku+(((k+1)-1)*Nxy+ij-1)*Nuxyz
           ii   = ijku+(((k  )-1)*Nxy+ij-1)*Nuxyz
           iim1 = ijku+(((k-1)-1)*Nxy+ij-1)*Nuxyz
           iim2 = ijku+(((k-2)-1)*Nxy+ij-1)*Nuxyz
           iim3 = ijku+(((k-3)-1)*Nxy+ij-1)*Nuxyz       
!
           fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),f(iip2),f(iip3),dz)
           fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dz2)
         enddo      
!-------------------------------------------------------------------------------
!  for k = 1
         iip3 = ijku+((4-1)*Nxy+ij-1)*Nuxyz
         iip2 = ijku+((3-1)*Nxy+ij-1)*Nuxyz
         iip1 = ijku+((2-1)*Nxy+ij-1)*Nuxyz
         ii   = ijku+((1-1)*Nxy+ij-1)*Nuxyz  
!
         fp(ii)  = FCD5th(fm3(ii1),fm2(ii1),fm1(ii1),f(iip1),f(iip2),f(iip3),dz)
         fpp(ii) = FSD5th(fm3(ii1),fm2(ii1),fm1(ii1),f(ii),f(iip1),f(iip2),f(iip3),dz2)
!-------------------------------------------------------------------------------
!  for k = 2
         iip3 = ijku+((5-1)*Nxy+ij-1)*Nuxyz
         iip2 = ijku+((4-1)*Nxy+ij-1)*Nuxyz
         iip1 = ijku+((3-1)*Nxy+ij-1)*Nuxyz
         ii   = ijku+((2-1)*Nxy+ij-1)*Nuxyz
         iim1 = ijku+((1-1)*Nxy+ij-1)*Nuxyz    
!
         fp(ii)  = FCD5th(fm2(ii1),fm1(ii1),f(iim1),f(iip1),f(iip2),f(iip3),dz)
         fpp(ii) = FSD5th(fm2(ii1),fm1(ii1),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dz2)
!-------------------------------------------------------------------------------
!  for k = 3
         iip3 = ijku+((6-1)*Nxy+ij-1)*Nuxyz
         iip2 = ijku+((5-1)*Nxy+ij-1)*Nuxyz
         iip1 = ijku+((4-1)*Nxy+ij-1)*Nuxyz
         ii   = ijku+((3-1)*Nxy+ij-1)*Nuxyz
         iim1 = ijku+((2-1)*Nxy+ij-1)*Nuxyz 
         iim2 = ijku+((1-1)*Nxy+ij-1)*Nuxyz   
!
         fp(ii)  = FCD5th(fm1(ii1),f(iim2),f(iim1),f(iip1),f(iip2),f(iip3),dz)
         fpp(ii) = FSD5th(fm1(ii1),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dz2)
!-------------------------------------------------------------------------------
!  for k = Nz-2
         iip2 = ijku+(((Nz  )-1)*Nxy+ij-1)*Nuxyz
         iip1 = ijku+(((Nz-1)-1)*Nxy+ij-1)*Nuxyz
         ii   = ijku+(((Nz-2)-1)*Nxy+ij-1)*Nuxyz
         iim1 = ijku+(((Nz-3)-1)*Nxy+ij-1)*Nuxyz
         iim2 = ijku+(((Nz-4)-1)*Nxy+ij-1)*Nuxyz
         iim3 = ijku+(((Nz-5)-1)*Nxy+ij-1)*Nuxyz             
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),f(iip2),fp1(ii1),dz)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),fp1(ii1),dz2)
!-------------------------------------------------------------------------------
!  for k = Nz-1
         iip1 = ijku+(((Nz  )-1)*Nxy+ij-1)*Nuxyz
         ii   = ijku+(((Nz-1)-1)*Nxy+ij-1)*Nuxyz
         iim1 = ijku+(((Nz-2)-1)*Nxy+ij-1)*Nuxyz
         iim2 = ijku+(((Nz-3)-1)*Nxy+ij-1)*Nuxyz
         iim3 = ijku+(((Nz-4)-1)*Nxy+ij-1)*Nuxyz             
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),f(iip1),fp1(ii1),fp2(ii1),dz)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),fp1(ii1),fp2(ii1),dz2)
!-------------------------------------------------------------------------------
!  for k = Nz
         ii   = ijku+(((Nz  )-1)*Nxy+ij-1)*Nuxyz
         iim1 = ijku+(((Nz-1)-1)*Nxy+ij-1)*Nuxyz
         iim2 = ijku+(((Nz-2)-1)*Nxy+ij-1)*Nuxyz
         iim3 = ijku+(((Nz-3)-1)*Nxy+ij-1)*Nuxyz             
!
         fp(ii)  = FCD5th(f(iim3),f(iim2),f(iim1),fp1(ii1),fp2(ii1),fp3(ii1),dz)
         fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),fp1(ii1),fp2(ii1),fp3(ii1),dz2) 
 30    continue
!
 35  continue
!$OMP end parallel do     
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm)  
   end subroutine PCDz5th6D_MPI
!===============================================================================
!===============================================================================
!  return the Boundary values along the x-direction
!
   subroutine getBCx6D(f,fBC,ix,Nx,Ny,Nz,Nux,Nuy,Nuz)
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fBC(1)
     Nxy = Nx*Ny
     Nuxyz = Nux*Nuy*Nuz
!
     do 10 k = 1, Nz
     do 10 j = 1, Ny
       ii = j+(k-1)*Ny
       ijk = ix+(j-1)*Nx+(k-1)*Nxy
!
       iis1 = (ijk-1)*Nuxyz+1
       iie1 = (ijk  )*Nuxyz
!
       iis2 = (ii-1)*Nuxyz+1
       iie2 = (ii  )*Nuxyz      
!
       fBC(iis2:iie2) = f(iis1:iie1)        
 10  continue
   end subroutine getBCx6D
!===============================================================================
!===============================================================================
!  return the Boundary values along the y-direction
!
   subroutine getBCy6D(f,fBC,jy,Nx,Ny,Nz,Nux,Nuy,Nuz)
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fBC(1)
     Nxy = Nx*Ny
     Nuxyz = Nux*Nuy*Nuz
!
     do 10 k = 1, Nz
     do 10 i = 1, Nx
       ii = i+(k-1)*Nx
       ijk = i+(jy-1)*Nx+(k-1)*Nxy
!
       iis1 = (ijk-1)*Nuxyz+1
       iie1 = (ijk  )*Nuxyz
!
       iis2 = (ii-1)*Nuxyz+1
       iie2 = (ii  )*Nuxyz
!
       fBC(iis2:iie2) = f(iis1:iie1)        
 10  continue
   end subroutine getBCy6D
!===============================================================================
!===============================================================================
!  return the Boundary values along the z-direction
!
   subroutine getBCz6D(f,fBC,kz,Nx,Ny,Nz,Nux,Nuy,Nuz)
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fBC(1)
     Nxy = Nx*Ny
     Nuxyz = Nux*Nuy*Nuz
!
     do 10 j = 1, Ny
     do 10 i = 1, Nx
       ii = i+(j-1)*Nx
       ijk = i+(j-1)*Nx+(kz-1)*Nxy
!
       iis1 = (ijk-1)*Nuxyz+1
       iie1 = (ijk  )*Nuxyz
!
       iis2 = (ii-1)*Nuxyz+1
       iie2 = (ii  )*Nuxyz
!
       fBC(iis2:iie2) = f(iis1:iie1)       
 10  continue
   end subroutine getBCz6D
!===============================================================================
!===============================================================================
!  This is the step of tenth-order central finite difference
!    to calculate the first derivative fp along x-direction
!    for free boundary condition.
!    "No MPI version"
!-------------------------------------------------------------------------------
!  input data    : dx
!  input arrays  : f(Nxyz)
!  output arrays : fp(Nxyz)
!-------------------------------------------------------------------------------
   subroutine FCDx10th3D(f,fp,dx,Nx,Ny,Nz)
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1)
     Nxy = Nx*Ny                                       
!-------------------------------------------------------------------------------
     do 30 k = 1, Nz
     do 30 j = 1, Ny
       ii = j+(k-1)*Ny
       jk = (j-1)*Nx+(k-1)*Nxy
!      
       do i = 6, Nx-5
         ijkp5 = (i+5)+jk
         ijkp4 = (i+4)+jk
         ijkp3 = (i+3)+jk
         ijkp2 = (i+2)+jk
         ijkp1 = (i+1)+jk
         ijk   = (i  )+jk
         ijkm1 = (i-1)+jk
         ijkm2 = (i-2)+jk
         ijkm3 = (i-3)+jk
         ijkm4 = (i-4)+jk
         ijkm5 = (i-5)+jk
         fp(ijk) = FCD10th(f(ijkm5),f(ijkm4),f(ijkm3),f(ijkm2),f(ijkm1), &
                           f(ijkp1),f(ijkp2),f(ijkp3),f(ijkp4),f(ijkp5),dx)
       enddo      
!-------------------------------------------------------------------------------
       ijk11 = 11+jk
       ijk10 = 10+jk
       ijk9  =  9+jk
       ijk8  =  8+jk
       ijk7  =  7+jk
       ijk6  =  6+jk
       ijk5  =  5+jk
       ijk4  =  4+jk
       ijk3  =  3+jk
       ijk2  =  2+jk
       ijk1  =  1+jk
!
!      for i = 1
       fp(ijk1) = FCD10th_B1(f(ijk1),f(ijk2),f(ijk3),f(ijk4),f(ijk5),f(ijk6), &
                             f(ijk7),f(ijk8),f(ijk9),f(ijk10),f(ijk11),dx)
!
!      for i = 2
       fp(ijk2) = FCD10th_B2(f(ijk1),f(ijk2),f(ijk3),f(ijk4),f(ijk5),f(ijk6), &
                             f(ijk7),f(ijk8),f(ijk9),f(ijk10),f(ijk11),dx)
!
!      for i = 3
       fp(ijk3) = FCD10th_B3(f(ijk1),f(ijk2),f(ijk3),f(ijk4),f(ijk5),f(ijk6), &
                             f(ijk7),f(ijk8),f(ijk9),f(ijk10),f(ijk11),dx)
!
!      for i = 4
       fp(ijk4) = FCD10th_B4(f(ijk1),f(ijk2),f(ijk3),f(ijk4),f(ijk5),f(ijk6), &
                             f(ijk7),f(ijk8),f(ijk9),f(ijk10),f(ijk11),dx)
!
!      for i = 5
       fp(ijk5) = FCD10th_B5(f(ijk1),f(ijk2),f(ijk3),f(ijk4),f(ijk5),f(ijk6), &
                             f(ijk7),f(ijk8),f(ijk9),f(ijk10),f(ijk11),dx)
!-------------------------------------------------------------------------------
       ijkN    = (Nx   )+jk
       ijkNm1  = (Nx-1 )+jk
       ijkNm2  = (Nx-2 )+jk
       ijkNm3  = (Nx-3 )+jk
       ijkNm4  = (Nx-4 )+jk
       ijkNm5  = (Nx-5 )+jk
       ijkNm6  = (Nx-6 )+jk
       ijkNm7  = (Nx-7 )+jk
       ijkNm8  = (Nx-8 )+jk
       ijkNm9  = (Nx-9 )+jk
       ijkNm10 = (Nx-10)+jk              
!
!      for i = Nx
       fp(ijkN) = FCD10th_BN(f(ijkNm10),f(ijkNm9),f(ijkNm8),f(ijkNm7),f(ijkNm6), &
                             f(ijkNm5),f(ijkNm4),f(ijkNm3),f(ijkNm2),f(ijkNm1),f(ijkN),dx)
!
!      for i = Nx-1
       fp(ijkNm1) = FCD10th_BN1(f(ijkNm10),f(ijkNm9),f(ijkNm8),f(ijkNm7),f(ijkNm6), &
                                f(ijkNm5),f(ijkNm4),f(ijkNm3),f(ijkNm2),f(ijkNm1),f(ijkN),dx)
!
!      for i = Nx-2
       fp(ijkNm2) = FCD10th_BN2(f(ijkNm10),f(ijkNm9),f(ijkNm8),f(ijkNm7),f(ijkNm6), &
                                f(ijkNm5),f(ijkNm4),f(ijkNm3),f(ijkNm2),f(ijkNm1),f(ijkN),dx)
!
!      for i = Nx-3
       fp(ijkNm3) = FCD10th_BN3(f(ijkNm10),f(ijkNm9),f(ijkNm8),f(ijkNm7),f(ijkNm6), &
                                f(ijkNm5),f(ijkNm4),f(ijkNm3),f(ijkNm2),f(ijkNm1),f(ijkN),dx)
!
!      for i = Nx-4
       fp(ijkNm4) = FCD10th_BN4(f(ijkNm10),f(ijkNm9),f(ijkNm8),f(ijkNm7),f(ijkNm6), &
                                f(ijkNm5),f(ijkNm4),f(ijkNm3),f(ijkNm2),f(ijkNm1),f(ijkN),dx)
!-------------------------------------------------------------------------------
 30  continue
!-------------------------------------------------------------------------------
   end subroutine FCDx10th3D
!===============================================================================
!===============================================================================
!  This is the step of tenth-order central finite difference
!    to calculate the first derivative fp along y-direction
!    for free boundary condition.
!    "No MPI version"
!-------------------------------------------------------------------------------
   subroutine FCDy10th3D(f,fp,dy,Nx,Ny,Nz)
!-------------------------------------------------------------------------------
!  input data    : dy
!  input arrays  : f(Nxyz)
!  output arrays : fp(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1)
     Nxy  = Nx*Ny
     Nxz  = Nx*Nz
!-------------------------------------------------------------------------------
     do 30 k = 1, Nz
     do 30 i = 1, Nx
       ii = i+(k-1)*Nx
       ik = i+(k-1)*Nxy
!      
       do j = 6, Ny-5
         ijkp5 = ((j+5)-1)*Nx+ik
         ijkp4 = ((j+4)-1)*Nx+ik
         ijkp3 = ((j+3)-1)*Nx+ik
         ijkp2 = ((j+2)-1)*Nx+ik
         ijkp1 = ((j+1)-1)*Nx+ik
         ijk   = ((j  )-1)*Nx+ik
         ijkm1 = ((j-1)-1)*Nx+ik
         ijkm2 = ((j-2)-1)*Nx+ik
         ijkm3 = ((j-3)-1)*Nx+ik
         ijkm4 = ((j-4)-1)*Nx+ik
         ijkm5 = ((j-5)-1)*Nx+ik       
!
         fp(ijk) = FCD10th(f(ijkm5),f(ijkm4),f(ijkm3),f(ijkm2),f(ijkm1), &
                           f(ijkp1),f(ijkp2),f(ijkp3),f(ijkp4),f(ijkp5),dy)
       enddo      
!-------------------------------------------------------------------------------
       ijk11 = (11-1)*Nx+ik
       ijk10 = (10-1)*Nx+ik
       ijk9  = ( 9-1)*Nx+ik
       ijk8  = ( 8-1)*Nx+ik
       ijk7  = ( 7-1)*Nx+ik
       ijk6  = ( 6-1)*Nx+ik
       ijk5  = ( 5-1)*Nx+ik
       ijk4  = ( 4-1)*Nx+ik
       ijk3  = ( 3-1)*Nx+ik
       ijk2  = ( 2-1)*Nx+ik
       ijk1  = ( 1-1)*Nx+ik  
!
!      for j = 1
       fp(ijk1) = FCD10th_B1(f(ijk1),f(ijk2),f(ijk3),f(ijk4),f(ijk5),f(ijk6), &
                             f(ijk7),f(ijk8),f(ijk9),f(ijk10),f(ijk11),dy)
!
!      for j = 2
       fp(ijk2) = FCD10th_B2(f(ijk1),f(ijk2),f(ijk3),f(ijk4),f(ijk5),f(ijk6), &
                             f(ijk7),f(ijk8),f(ijk9),f(ijk10),f(ijk11),dy)
!
!      for j = 3
       fp(ijk3) = FCD10th_B3(f(ijk1),f(ijk2),f(ijk3),f(ijk4),f(ijk5),f(ijk6), &
                             f(ijk7),f(ijk8),f(ijk9),f(ijk10),f(ijk11),dy)
!
!      for j = 4
       fp(ijk4) = FCD10th_B4(f(ijk1),f(ijk2),f(ijk3),f(ijk4),f(ijk5),f(ijk6), &
                             f(ijk7),f(ijk8),f(ijk9),f(ijk10),f(ijk11),dy)
!
!      for j = 5
       fp(ijk5) = FCD10th_B5(f(ijk1),f(ijk2),f(ijk3),f(ijk4),f(ijk5),f(ijk6), &
                             f(ijk7),f(ijk8),f(ijk9),f(ijk10),f(ijk11),dy)
!-------------------------------------------------------------------------------
!      for j = Ny
       ijkN    = ((Ny   )-1)*Nx+ik
       ijkNm1  = ((Ny-1 )-1)*Nx+ik
       ijkNm2  = ((Ny-2 )-1)*Nx+ik
       ijkNm3  = ((Ny-3 )-1)*Nx+ik
       ijkNm4  = ((Ny-4 )-1)*Nx+ik
       ijkNm5  = ((Ny-5 )-1)*Nx+ik
       ijkNm6  = ((Ny-6 )-1)*Nx+ik
       ijkNm7  = ((Ny-7 )-1)*Nx+ik
       ijkNm8  = ((Ny-8 )-1)*Nx+ik
       ijkNm9  = ((Ny-9 )-1)*Nx+ik
       ijkNm10 = ((Ny-10)-1)*Nx+ik             
!
!      for j = Ny
       fp(ijkN) = FCD10th_BN(f(ijkNm10),f(ijkNm9),f(ijkNm8),f(ijkNm7),f(ijkNm6), &
                             f(ijkNm5),f(ijkNm4),f(ijkNm3),f(ijkNm2),f(ijkNm1),f(ijkN),dy)
!
!      for j = Ny-1
       fp(ijkNm1) = FCD10th_BN1(f(ijkNm10),f(ijkNm9),f(ijkNm8),f(ijkNm7),f(ijkNm6), &
                                f(ijkNm5),f(ijkNm4),f(ijkNm3),f(ijkNm2),f(ijkNm1),f(ijkN),dy)
!
!      for j = Ny-2
       fp(ijkNm2) = FCD10th_BN2(f(ijkNm10),f(ijkNm9),f(ijkNm8),f(ijkNm7),f(ijkNm6), &
                                f(ijkNm5),f(ijkNm4),f(ijkNm3),f(ijkNm2),f(ijkNm1),f(ijkN),dy)
!
!      for j = Ny-3
       fp(ijkNm3) = FCD10th_BN3(f(ijkNm10),f(ijkNm9),f(ijkNm8),f(ijkNm7),f(ijkNm6), &
                                f(ijkNm5),f(ijkNm4),f(ijkNm3),f(ijkNm2),f(ijkNm1),f(ijkN),dy)
!
!      for j = Ny-4
       fp(ijkNm4) = FCD10th_BN4(f(ijkNm10),f(ijkNm9),f(ijkNm8),f(ijkNm7),f(ijkNm6), &
                                f(ijkNm5),f(ijkNm4),f(ijkNm3),f(ijkNm2),f(ijkNm1),f(ijkN),dy)
!-------------------------------------------------------------------------------
 30  continue
!-------------------------------------------------------------------------------
   end subroutine FCDy10th3D
!===============================================================================
!===============================================================================
!  This is the step of the fifth-order central difference
!    to calculate the first derivative fp along z-direction
!    for free boundary condition.
!    "No MPI version"
!-------------------------------------------------------------------------------
   subroutine FCDz10th3D(f,fp,dz,Nx,Ny,Nz)
!-------------------------------------------------------------------------------
!  input data    : dz
!  input arrays  : f(Nxyz)
!  output arrays : fp(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1)
     Nxy  = Nx*Ny                                     
!-------------------------------------------------------------------------------
     do 30 j = 1, Ny
     do 30 i = 1, Nx
       ii = i+(j-1)*Nx
       ij = ii !i+(j-1)*Nx
!      
       do k = 4, Nz-3
         ijkp5 = ((k+5)-1)*Nxy+ij
         ijkp4 = ((k+4)-1)*Nxy+ij
         ijkp3 = ((k+3)-1)*Nxy+ij
         ijkp2 = ((k+2)-1)*Nxy+ij
         ijkp1 = ((k+1)-1)*Nxy+ij
         ijk   = ((k  )-1)*Nxy+ij
         ijkm1 = ((k-1)-1)*Nxy+ij
         ijkm2 = ((k-2)-1)*Nxy+ij
         ijkm3 = ((k-3)-1)*Nxy+ij
         ijkm4 = ((k-4)-1)*Nxy+ij
         ijkm5 = ((k-5)-1)*Nxy+ij       
!
         fp(ijk) = FCD10th(f(ijkm5),f(ijkm4),f(ijkm3),f(ijkm2),f(ijkm1), &
                           f(ijkp1),f(ijkp2),f(ijkp3),f(ijkp4),f(ijkp5),dz)
       enddo      
!-------------------------------------------------------------------------------
!      for k = 1
       ijk11 = (11-1)*Nxy+ij
       ijk10 = (10-1)*Nxy+ij
       ijk9  = (9 -1)*Nxy+ij
       ijk8  = (8 -1)*Nxy+ij
       ijk7  = (7 -1)*Nxy+ij
       ijk6  = (6 -1)*Nxy+ij
       ijk5  = (5 -1)*Nxy+ij
       ijk4  = (4 -1)*Nxy+ij
       ijk3  = (3 -1)*Nxy+ij
       ijk2  = (2 -1)*Nxy+ij
       ijk1  = (1 -1)*Nxy+ij  
!
!
!      for k = 1
       fp(ijk1) = FCD10th_B1(f(ijk1),f(ijk2),f(ijk3),f(ijk4),f(ijk5),f(ijk6), &
                             f(ijk7),f(ijk8),f(ijk9),f(ijk10),f(ijk11),dz)
!
!      for k = 2
       fp(ijk2) = FCD10th_B2(f(ijk1),f(ijk2),f(ijk3),f(ijk4),f(ijk5),f(ijk6), &
                             f(ijk7),f(ijk8),f(ijk9),f(ijk10),f(ijk11),dz)
!
!      for k = 3
       fp(ijk3) = FCD10th_B3(f(ijk1),f(ijk2),f(ijk3),f(ijk4),f(ijk5),f(ijk6), &
                             f(ijk7),f(ijk8),f(ijk9),f(ijk10),f(ijk11),dz)
!
!      for k = 4
       fp(ijk4) = FCD10th_B4(f(ijk1),f(ijk2),f(ijk3),f(ijk4),f(ijk5),f(ijk6), &
                             f(ijk7),f(ijk8),f(ijk9),f(ijk10),f(ijk11),dz)
!
!      for k = 5
       fp(ijk5) = FCD10th_B5(f(ijk1),f(ijk2),f(ijk3),f(ijk4),f(ijk5),f(ijk6), &
                             f(ijk7),f(ijk8),f(ijk9),f(ijk10),f(ijk11),dz)
!-------------------------------------------------------------------------------
!      for k = Nz
       ijkN    = ((Nz   )-1)*Nxy+ij
       ijkNm1  = ((Nz- 1)-1)*Nxy+ij
       ijkNm2  = ((Nz- 2)-1)*Nxy+ij
       ijkNm3  = ((Nz- 3)-1)*Nxy+ij
       ijkNm4  = ((Nz- 4)-1)*Nxy+ij
       ijkNm5  = ((Nz- 5)-1)*Nxy+ij
       ijkNm6  = ((Nz- 6)-1)*Nxy+ij
       ijkNm7  = ((Nz- 7)-1)*Nxy+ij
       ijkNm8  = ((Nz- 8)-1)*Nxy+ij
       ijkNm9  = ((Nz- 9)-1)*Nxy+ij
       ijkNm10 = ((Nz-10)-1)*Nxy+ij             
!
!      for k = Ny
       fp(ijkN) = FCD10th_BN(f(ijkNm10),f(ijkNm9),f(ijkNm8),f(ijkNm7),f(ijkNm6), &
                             f(ijkNm5),f(ijkNm4),f(ijkNm3),f(ijkNm2),f(ijkNm1),f(ijkN),dz)
!
!      for k = Ny-1
       fp(ijkNm1) = FCD10th_BN1(f(ijkNm10),f(ijkNm9),f(ijkNm8),f(ijkNm7),f(ijkNm6), &
                                f(ijkNm5),f(ijkNm4),f(ijkNm3),f(ijkNm2),f(ijkNm1),f(ijkN),dz)
!
!      for k = Ny-2
       fp(ijkNm2) = FCD10th_BN2(f(ijkNm10),f(ijkNm9),f(ijkNm8),f(ijkNm7),f(ijkNm6), &
                                f(ijkNm5),f(ijkNm4),f(ijkNm3),f(ijkNm2),f(ijkNm1),f(ijkN),dz)
!
!      for k = Ny-3
       fp(ijkNm3) = FCD10th_BN3(f(ijkNm10),f(ijkNm9),f(ijkNm8),f(ijkNm7),f(ijkNm6), &
                                f(ijkNm5),f(ijkNm4),f(ijkNm3),f(ijkNm2),f(ijkNm1),f(ijkN),dz)
!
!      for k = Ny-4
       fp(ijkNm4) = FCD10th_BN4(f(ijkNm10),f(ijkNm9),f(ijkNm8),f(ijkNm7),f(ijkNm6), &
                                f(ijkNm5),f(ijkNm4),f(ijkNm3),f(ijkNm2),f(ijkNm1),f(ijkN),dz)
 30  continue     
!-------------------------------------------------------------------------------
   end subroutine FCDz10th3D
!===============================================================================
!===============================================================================
   function FCD10th(fm5,fm4,fm3,fm2,fm1,fp1,fp2,fp3,fp4,fp5,dx)
     implicit double precision (A-H,O-Z)
     a1 = 5.d0/6.d0    ! 2100/2520
     a2 = 5.d0/21.d0   !  600/2520
     a3 = 5.d0/84.d0   !  150/2520
     a4 = 5.d0/504.d0  !   25/2520
     a5 = 1.d0/1260.d0 !    2/2520
!
     sum1 = sumABCDE(-a5*fm5, a4*fm4,-a3*fm3, a2*fm2,-a1*fm1)
     sum2 = sumABCDE( a1*fp1,-a2*fp2, a3*fp3,-a4*fp4, a5*fp5)
     FCD10th = sumAB(sum1,sum2)/dx
   end function FCD10th
!===============================================================================
!===============================================================================
   function FCD10th_B1(f,fp1,fp2,fp3,fp4,fp5,fp6,fp7,fp8,fp9,fp10,dx)
     implicit double precision (A-H,O-Z)
     a0 = -7381.d0
     a1 = 53200.d0
     a2 = -56700.d0
     a3 = 100800.d0
     a4 = -132300.d0
     a5 = 127008.d0
     a6 = -88200.d0
     a7 = 43200.d0
     a8 = -14175.d0
     a9 = 2800.d0
     a10 = -252.d0
!
     sum1 = sumABCDE(a0*f,a1*fp1,a2*fp2,a3*fp3,a4*fp4)
     sum2 = sumABC(a5*fp5,a6*fp6,a7*fp7)
     sum3 = sumABC(a8*fp8,a9*fp9,a10*fp10)
     FCD10th_B1 = sumABC(sum1,sum2,sum3)/dx/2520.d0
   end function FCD10th_B1
!===============================================================================
!===============================================================================
   function FCD10th_B2(fm1,f,fp1,fp2,fp3,fp4,fp5,fp6,fp7,fp8,fp9,dx)
     implicit double precision (A-H,O-Z)
     a0 = -252.d0
     a1 = -4609.d0
     a2 = 11340.d0
     a3 = -15120.d0
     a4 = 17640.d0
     a5 = -15876.d0
     a6 = 10584.d0
     a7 = -5040.d0
     a8 = 1620.d0
     a9 = -315.d0
     a10 = 28.d0
!
     sum1 = sumABCDE(a0*fm1,a1*f,a2*fp1,a3*fp2,a4*fp3)
     sum2 = sumABC(a5*fp4,a6*fp5,a7*fp6)
     sum3 = sumABC(a8*fp7,a9*fp8,a10*fp9)
     FCD10th_B2 = sumABC(sum1,sum2,sum3)/dx/2520.d0
   end function FCD10th_B2
!===============================================================================
!===============================================================================
   function FCD10th_B3(fm2,fm1,f,fp1,fp2,fp3,fp4,fp5,fp6,fp7,fp8,dx)
     implicit double precision (A-H,O-Z)
     a0 = 28.d0
     a1 = -560.d0
     a2 = -3609.d0
     a3 = 6720.d0
     a4 = -5880.d0
     a5 = 4704.d0
     a6 = -2940.d0
     a7 = 1344.d0
     a8 = -420.d0
     a9 = 80.d0
     a10 = -7.d0
!
     sum1 = sumABCDE(a0*fm2,a1*fm1,a2*f,a3*fp1,a4*fp2)
     sum2 = sumABC(a5*fp3,a6*fp4,a7*fp5)
     sum3 = sumABC(a8*fp6,a9*fp7,a10*fp8)
     FCD10th_B3 = sumABC(sum1,sum2,sum3)/dx/2520.d0
   end function FCD10th_B3
!===============================================================================
!===============================================================================
   function FCD10th_B4(fm3,fm2,fm1,f,fp1,fp2,fp3,fp4,fp5,fp6,fp7,dx)
     implicit double precision (A-H,O-Z)
     a0 = -7.d0
     a1 = 105.d0
     a2 = -945.d0
     a3 = -1914.d0
     a4 = 4410.d0
     a5 = -2646.d0
     a6 = 1470.d0
     a7 = -630.d0
     a8 = 189.d0
     a9 = -35.d0
     a10 = 3.d0
!
     sum1 = sumABCDE(a0*fm3,a1*fm2,a2*fm1,a3*f,a4*fp1)
     sum2 = sumABC(a5*fp2,a6*fp3,a7*fp4)
     sum3 = sumABC(a8*fp5,a9*fp6,a10*fp7)
     FCD10th_B4 = sumABC(sum1,sum2,sum3)/dx/2520.d0
   end function FCD10th_B4
!===============================================================================
!===============================================================================
   function FCD10th_B5(fm4,fm3,fm2,fm1,f,fp1,fp2,fp3,fp4,fp5,fp6,dx)
     implicit double precision (A-H,O-Z)
     a0 = 3.d0
     a1 = -40.d0
     a2 = 270.d0
     a3 = -1440.d0
     a4 = -924.d0
     a5 = 3024.d0
     a6 = -1260.d0
     a7 = 480.d0
     a8 = -135.d0
     a9 = 24.d0
     a10 = -2.d0
!
     sum1 = sumABCDE(a0*fm4,a1*fm3,a2*fm2,a3*fm1,a4*f)
     sum2 = sumABC(a5*fp1,a6*fp2,a7*fp3)
     sum3 = sumABC(a8*fp4,a9*fp5,a10*fp6)
     FCD10th_B5 = sumABC(sum1,sum2,sum3)/dx/2520.d0
   end function FCD10th_B5
!===============================================================================
!===============================================================================
   function FCD10th_BN(fm10,fm9,fm8,fm7,fm6,fm5,fm4,fm3,fm2,fm1,f,dx)
     implicit double precision (A-H,O-Z)
     a0 = 252.d0
     a1 = -2800.d0
     a2 = 14175.d0
     a3 = -43200.d0
     a4 = 88200.d0
     a5 = -127008.d0
     a6 = 132300.d0
     a7 = -100800.d0
     a8 = 56700.d0
     a9 = -25200.d0
     a10 = 7381.d0
!
     sum1 = sumABCDE(a0*fm10,a1*fm9,a2*fm8,a3*fm7,a4*fm6)
     sum2 = sumABC(a5*fm5,a6*fm4,a7*fm3)
     sum3 = sumABC(a8*fm2,a9*fm1,a10*f)
     FCD10th_BN = sumABC(sum1,sum2,sum3)/dx/2520.d0
   end function FCD10th_BN
!===============================================================================
!===============================================================================
   function FCD10th_BN1(fm9,fm8,fm7,fm6,fm5,fm4,fm3,fm2,fm1,f,fp1,dx)
     implicit double precision (A-H,O-Z)
     a0 = -28.d0
     a1 = 315.d0
     a2 = -1620.d0
     a3 = 5040.d0
     a4 = -10584.d0
     a5 = 15876.d0
     a6 = -17640.d0
     a7 = 15120.d0
     a8 = -11340.d0
     a9 = 4609.d0
     a10 = 252.d0
!
     sum1 = sumABCDE(a0*fm9,a1*fm8,a2*fm7,a3*fm6,a4*fm5)
     sum2 = sumABC(a5*fm4,a6*fm3,a7*fm2)
     sum3 = sumABC(a8*fm1,a9*f,a10*fp1)
     FCD10th_BN1 = sumABC(sum1,sum2,sum3)/dx/2520.d0
   end function FCD10th_BN1
!===============================================================================
!===============================================================================
   function FCD10th_BN2(fm8,fm7,fm6,fm5,fm4,fm3,fm2,fm1,f,fp1,fp2,dx)
     implicit double precision (A-H,O-Z)
     a0 = 7.d0
     a1 = -80.d0
     a2 = 420.d0
     a3 = -1344.d0
     a4 = 2940.d0
     a5 = -4704.d0
     a6 = 5880.d0
     a7 = -6720.d0
     a8 = 3069.d0
     a9 = 560.d0
     a10 = -28.d0
!
     sum1 = sumABCDE(a0*fm8,a1*fm7,a2*fm6,a3*fm5,a4*fm4)
     sum2 = sumABC(a5*fm3,a6*fm2,a7*fm1)
     sum3 = sumABC(a8*f,a9*fp1,a10*fp2)
     FCD10th_BN2 = sumABC(sum1,sum2,sum3)/dx/2520.d0
   end function FCD10th_BN2
!===============================================================================
!===============================================================================
   function FCD10th_BN3(fm7,fm6,fm5,fm4,fm3,fm2,fm1,f,fp1,fp2,fp3,dx)
     implicit double precision (A-H,O-Z)
     a0 = -3.d0
     a1 = 35.d0
     a2 = -189.d0
     a3 = 630.d0
     a4 = -1470.d0
     a5 = 2646.d0
     a6 = -4410.d0
     a7 = 1914.d0
     a8 = 945.d0
     a9 = -105.d0
     a10 = 7.d0
!
     sum1 = sumABCDE(a0*fm7,a1*fm6,a2*fm5,a3*fm4,a4*fm3)
     sum2 = sumABC(a5*fm2,a6*fm1,a7*f)
     sum3 = sumABC(a8*fp1,a9*fp2,a10*fp3)
     FCD10th_BN3 = sumABC(sum1,sum2,sum3)/dx/2520.d0
   end function FCD10th_BN3
!===============================================================================
!===============================================================================
   function FCD10th_BN4(fm6,fm5,fm4,fm3,fm2,fm1,f,fp1,fp2,fp3,fp4,dx)
     implicit double precision (A-H,O-Z)
     a0 = 2.d0
     a1 = -24.d0
     a2 = 135.d0
     a3 = -480.d0
     a4 = 1260.d0
     a5 = -3024.d0
     a6 = 924.d0
     a7 = 1440.d0
     a8 = -270.d0
     a9 = 40.d0
     a10 = -3.d0
!
     sum1 = sumABCDE(a0*fm6,a1*fm5,a2*fm4,a3*fm3,a4*fm2)
     sum2 = sumABC(a5*fm1,a6*f,a7*fp1)
     sum3 = sumABC(a8*fp2,a9*fp3,a10*fp4)
     FCD10th_BN4 = sumABC(sum1,sum2,sum3)/dx/2520.d0
   end function FCD10th_BN4
!===============================================================================
!===============================================================================
!===============================================================================
!  This is the step of fifth-order interpolation method
!    for free boundary condition.
!-------------------------------------------------------------------------------
!  input arrays  : f(Nxyz)
!  output arrays : fINTP(Nxyz)
!-------------------------------------------------------------------------------
   subroutine FINTPxm5th3D_MPI(f,fINTP,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fINTP(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy = Nx*Ny
     Nyz = Ny*Nz
     Nyz2 = 2*Nyz
     Nyz3 = 3*Nyz
     allocate (fm3(Nyz),fm2(Nyz),fm1(Nyz),fp1(Nyz),fp2(Nyz),fp3(Nyz))
     allocate (fsend(Nyz3),frecvp(Nyz3),frecvm(Nyz3))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCx3D(f(1),fsend(     1),1,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz +1),2,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nyz3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fp1(1:Nyz) = frecvp(1     :Nyz )
     fp2(1:Nyz) = frecvp(Nyz +1:Nyz2)
     fp3(1:Nyz) = frecvp(Nyz2+1:Nyz3)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCx3D(f(1),fsend(     1),Nx,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz +1),Nx-1,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz2+1),Nx-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nyz3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fm1(1:Nyz) = frecvm(1     :Nyz )
     fm2(1:Nyz) = frecvm(Nyz +1:Nyz2)
     fm3(1:Nyz) = frecvm(Nyz2+1:Nyz3)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
     if (iD == iDstart) then
       call getBCx3D(f(1),fm1(1),1,Nx,Ny,Nz)      
       call getBCx3D(f(1),fm2(1),1,Nx,Ny,Nz)
       call getBCx3D(f(1),fm3(1),1,Nx,Ny,Nz)      
     endif
!
     if (iD == iDend) then           
       call getBCx3D(f(1),fp1(1),Nx,Nx,Ny,Nz)      
       call getBCx3D(f(1),fp2(1),Nx,Nx,Ny,Nz)
       call getBCx3D(f(1),fp3(1),Nx,Nx,Ny,Nz)      
     endif                                         
!-------------------------------------------------------------------------------
     ii = 0
!$OMP parallel do private(ii,jk,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2) 
     do 30 k = 1, Nz
     do 30 j = 1, Ny
       ii = j+(k-1)*Ny
       jk = (j-1)*Nx+(k-1)*Nxy
!      
       do i = 4, Nx-2
!         ijkp3 = (i+3)+jk
         ijkp2 = (i+2)+jk
         ijkp1 = (i+1)+jk
         ijk   = (i  )+jk
         ijkm1 = (i-1)+jk
         ijkm2 = (i-2)+jk
         ijkm3 = (i-3)+jk
         fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
       enddo      
!-------------------------------------------------------------------------------
!     for i = 1
!       ijkp3 = 4+jk
       ijkp2 = 3+jk
       ijkp1 = 2+jk
       ijk   = 1+jk
!
       fINTP(ijk) = smooth5th(fm3(ii),fm2(ii),fm1(ii),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for i = 2
!       ijkp3 = 5+jk
       ijkp2 = 4+jk
       ijkp1 = 3+jk
       ijk   = 2+jk
       ijkm1 = 1+jk     
!
       fINTP(ijk) = smooth5th(fm2(ii),fm1(ii),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for i = 3
!       ijkp3 = 6+jk
       ijkp2 = 5+jk
       ijkp1 = 4+jk
       ijk   = 3+jk
       ijkm1 = 2+jk
       ijkm2 = 1+jk     
!
       fINTP(ijk) = smooth5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for i = Nx-2
!       ijkp2 = (Nx  )+jk
!       ijkp1 = (Nx-1)+jk
!       ijk   = (Nx-2)+jk
!       ijkm1 = (Nx-3)+jk
!       ijkm2 = (Nx-4)+jk
!       ijkm3 = (Nx-5)+jk              
!
!       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for i = Nx-1
       ijkp1 = (Nx  )+jk
       ijk   = (Nx-1)+jk
       ijkm1 = (Nx-2)+jk
       ijkm2 = (Nx-3)+jk
       ijkm3 = (Nx-4)+jk              
!
       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),fp1(ii))
!-------------------------------------------------------------------------------
!      for i = Nx
       ijk   = (Nx  )+jk
       ijkm1 = (Nx-1)+jk
       ijkm2 = (Nx-2)+jk
       ijkm3 = (Nx-3)+jk              
!
       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),fp1(ii),fp2(ii))
 30  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm)  
   end subroutine FINTPxm5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fifth-order interpolation method
!    for periodic boundary condition.
!-------------------------------------------------------------------------------
   subroutine PINTPxm5th3D_MPI(f,fINTP,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input arrays  : f(Nxyz)
!  output arrays : fINTP(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fINTP(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy  = Nx*Ny
     Nyz  = Ny*Nz
     Nyz2 = 2*Nyz
     Nyz3 = 3*Nyz
     allocate (fm3(Nyz),fm2(Nyz),fm1(Nyz),fp1(Nyz),fp2(Nyz),fp3(Nyz))
     allocate (fsend(Nyz3),frecvp(Nyz3),frecvm(Nyz3))
!-------------------------------------------------------------------------------
!  if P = 1, it doesn't transfer data in MPI
!
     if (iDend == iDstart) then       
       call getBCx3D(f(1),fp3(1),4,Nx,Ny,Nz)
       call getBCx3D(f(1),fp2(1),3,Nx,Ny,Nz)
       call getBCx3D(f(1),fp1(1),2,Nx,Ny,Nz)
!       
       call getBCx3D(f(1),fm1(1),Nx-1,Nx,Ny,Nz)
       call getBCx3D(f(1),fm2(1),Nx-2,Nx,Ny,Nz)
       call getBCx3D(f(1),fm3(1),Nx-3,Nx,Ny,Nz)
       goto 123
     endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCx3D(f(1),fsend(1     ),1,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz +1),2,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz2+1),3,Nx,Ny,Nz)
!
     call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nyz3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR) 
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCx3D(f(1),fsend(1     ),Nx  ,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz +1),Nx-1,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz2+1),Nx-2,Nx,Ny,Nz)
!
     call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nyz3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR) 
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
     if (iD == iDstart) then
!
       call getBCx3D(f(1),fsend(1     ),2,Nx,Ny,Nz)
       call getBCx3D(f(1),fsend(Nyz +1),3,Nx,Ny,Nz)
       call getBCx3D(f(1),fsend(Nyz2+1),4,Nx,Ny,Nz)
!
       call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDend,120, &
                         frecvm(1),Nyz3,MPI_REAL8,iDend,120, &
                         MPI_WORLD,ISTATUS,IERR)      
     endif
!
     if (iD == iDend) then      
!
       call getBCx3D(f(1),fsend(1     ),Nx-1,Nx,Ny,Nz)
       call getBCx3D(f(1),fsend(Nyz +1),Nx-2,Nx,Ny,Nz)
       call getBCx3D(f(1),fsend(Nyz2+1),Nx-3,Nx,Ny,Nz)
!
       call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDstart,120, &
                         frecvp(1),Nyz3,MPI_REAL8,iDstart,120, &
                         MPI_WORLD,ISTATUS,IERR)
     endif
!
     fp1(1:Nyz) = frecvp(1     :Nyz )
     fp2(1:Nyz) = frecvp(Nyz +1:Nyz2)
     fp3(1:Nyz) = frecvp(Nyz2+1:Nyz3)
     fm1(1:Nyz) = frecvm(1     :Nyz )
     fm2(1:Nyz) = frecvm(Nyz +1:Nyz2)
     fm3(1:Nyz) = frecvm(Nyz2+1:Nyz3)
!    
123  continue                                          
!-------------------------------------------------------------------------------
     ii = 0
!$OMP parallel do private(ii,jk,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 k = 1, Nz
     do 30 j = 1, Ny
       ii = j+(k-1)*Ny
       jk = (j-1)*Nx+(k-1)*Nxy
!      
       do i = 4, Nx-2
!         ijkp3 = (i+3)+jk
         ijkp2 = (i+2)+jk
         ijkp1 = (i+1)+jk
         ijk   = (i  )+jk
         ijkm1 = (i-1)+jk
         ijkm2 = (i-2)+jk
         ijkm3 = (i-3)+jk
         fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
       enddo      
!-------------------------------------------------------------------------------
!     for i = 1
!       ijkp3 = 4+jk
       ijkp2 = 3+jk
       ijkp1 = 2+jk
       ijk   = 1+jk
!
       fINTP(ijk) = smooth5th(fm3(ii),fm2(ii),fm1(ii),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for i = 2
!       ijkp3 = 5+jk
       ijkp2 = 4+jk
       ijkp1 = 3+jk
       ijk   = 2+jk
       ijkm1 = 1+jk     
!
       fINTP(ijk) = smooth5th(fm2(ii),fm1(ii),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for i = 3
!       ijkp3 = 6+jk
       ijkp2 = 5+jk
       ijkp1 = 4+jk
       ijk   = 3+jk
       ijkm1 = 2+jk
       ijkm2 = 1+jk     
!
       fINTP(ijk) = smooth5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for i = Nx-2
!       ijkp2 = (Nx  )+jk
!       ijkp1 = (Nx-1)+jk
!       ijk   = (Nx-2)+jk
!       ijkm1 = (Nx-3)+jk
!       ijkm2 = (Nx-4)+jk
!       ijkm3 = (Nx-5)+jk              
!
!       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for i = Nx-1
       ijkp1 = (Nx  )+jk
       ijk   = (Nx-1)+jk
       ijkm1 = (Nx-2)+jk
       ijkm2 = (Nx-3)+jk
       ijkm3 = (Nx-4)+jk              
!
       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),fp1(ii))
!-------------------------------------------------------------------------------
!      for i = Nx
       ijk   = (Nx  )+jk
       ijkm1 = (Nx-1)+jk
       ijkm2 = (Nx-2)+jk
       ijkm3 = (Nx-3)+jk              
!
       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),fp1(ii),fp2(ii))
 30  continue
!$OMP end parallel do 
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm) 
   end subroutine PINTPxm5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fifth-order interpolation method
!    for free boundary condition.
!-------------------------------------------------------------------------------
!  input arrays  : f(Nxyz)
!  output arrays : fINTP(Nxyz)
!-------------------------------------------------------------------------------
   subroutine FINTPxp5th3D_MPI(f,fINTP,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fINTP(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy = Nx*Ny
     Nyz = Ny*Nz
     Nyz2 = 2*Nyz
     Nyz3 = 3*Nyz
     allocate (fm3(Nyz),fm2(Nyz),fm1(Nyz),fp1(Nyz),fp2(Nyz),fp3(Nyz))
     allocate (fsend(Nyz3),frecvp(Nyz3),frecvm(Nyz3))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCx3D(f(1),fsend(     1),1,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz +1),2,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nyz3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fp1(1:Nyz) = frecvp(1     :Nyz )
     fp2(1:Nyz) = frecvp(Nyz +1:Nyz2)
     fp3(1:Nyz) = frecvp(Nyz2+1:Nyz3)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCx3D(f(1),fsend(     1),Nx,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz +1),Nx-1,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz2+1),Nx-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nyz3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fm1(1:Nyz) = frecvm(1     :Nyz )
     fm2(1:Nyz) = frecvm(Nyz +1:Nyz2)
     fm3(1:Nyz) = frecvm(Nyz2+1:Nyz3)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
     if (iD == iDstart) then
       call getBCx3D(f(1),fm1(1),1,Nx,Ny,Nz)      
       call getBCx3D(f(1),fm2(1),1,Nx,Ny,Nz)
       call getBCx3D(f(1),fm3(1),1,Nx,Ny,Nz)      
     endif
!
     if (iD == iDend) then           
       call getBCx3D(f(1),fp1(1),Nx,Nx,Ny,Nz)      
       call getBCx3D(f(1),fp2(1),Nx,Nx,Ny,Nz)
       call getBCx3D(f(1),fp3(1),Nx,Nx,Ny,Nz)      
     endif                                         
!-------------------------------------------------------------------------------
     ii = 0
!$OMP parallel do private(ii,jk,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2) 
     do 30 k = 1, Nz
     do 30 j = 1, Ny
       ii = j+(k-1)*Ny
       jk = (j-1)*Nx+(k-1)*Nxy
!      
       do i = 3, Nx-3
         ijkp3 = (i+3)+jk
         ijkp2 = (i+2)+jk
         ijkp1 = (i+1)+jk
         ijk   = (i  )+jk
         ijkm1 = (i-1)+jk
         ijkm2 = (i-2)+jk
!         ijkm3 = (i-3)+jk
         fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
       enddo      
!-------------------------------------------------------------------------------
!     for i = 1
       ijkp3 = 4+jk
       ijkp2 = 3+jk
       ijkp1 = 2+jk
       ijk   = 1+jk
!
       fINTP(ijk) = smooth5th(fm2(ii),fm1(ii),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for i = 2
       ijkp3 = 5+jk
       ijkp2 = 4+jk
       ijkp1 = 3+jk
       ijk   = 2+jk
       ijkm1 = 1+jk     
!
       fINTP(ijk) = smooth5th(fm1(ii),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for i = 3
!       ijkp3 = 6+jk
!       ijkp2 = 5+jk
!       ijkp1 = 4+jk
!       ijk   = 3+jk
!       ijkm1 = 2+jk
!       ijkm2 = 1+jk     
!
!       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for i = Nx-2
       ijkp2 = (Nx  )+jk
       ijkp1 = (Nx-1)+jk
       ijk   = (Nx-2)+jk
       ijkm1 = (Nx-3)+jk
       ijkm2 = (Nx-4)+jk
!       ijkm3 = (Nx-5)+jk              
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),fp1(ii))
!-------------------------------------------------------------------------------
!      for i = Nx-1
       ijkp1 = (Nx  )+jk
       ijk   = (Nx-1)+jk
       ijkm1 = (Nx-2)+jk
       ijkm2 = (Nx-3)+jk
!       ijkm3 = (Nx-4)+jk              
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),fp1(ii),fp2(ii))
!-------------------------------------------------------------------------------
!      for i = Nx
       ijk   = (Nx  )+jk
       ijkm1 = (Nx-1)+jk
       ijkm2 = (Nx-2)+jk
!       ijkm3 = (Nx-3)+jk              
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),fp1(ii),fp2(ii),fp3(ii))
 30  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm)  
   end subroutine FINTPxp5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fifth-order interpolation method
!    for periodic boundary condition.
!-------------------------------------------------------------------------------
   subroutine PINTPxp5th3D_MPI(f,fINTP,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input arrays  : f(Nxyz)
!  output arrays : fINTP(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fINTP(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy  = Nx*Ny
     Nyz  = Ny*Nz
     Nyz2 = 2*Nyz
     Nyz3 = 3*Nyz
     allocate (fm3(Nyz),fm2(Nyz),fm1(Nyz),fp1(Nyz),fp2(Nyz),fp3(Nyz))
     allocate (fsend(Nyz3),frecvp(Nyz3),frecvm(Nyz3))
!-------------------------------------------------------------------------------
!  if P = 1, it doesn't transfer data in MPI
!
     if (iDend == iDstart) then       
       call getBCx3D(f(1),fp3(1),4,Nx,Ny,Nz)
       call getBCx3D(f(1),fp2(1),3,Nx,Ny,Nz)
       call getBCx3D(f(1),fp1(1),2,Nx,Ny,Nz)
!       
       call getBCx3D(f(1),fm1(1),Nx-1,Nx,Ny,Nz)
       call getBCx3D(f(1),fm2(1),Nx-2,Nx,Ny,Nz)
       call getBCx3D(f(1),fm3(1),Nx-3,Nx,Ny,Nz)
       goto 123
     endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCx3D(f(1),fsend(1     ),1,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz +1),2,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz2+1),3,Nx,Ny,Nz)
!
     call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nyz3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR) 
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCx3D(f(1),fsend(1     ),Nx  ,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz +1),Nx-1,Nx,Ny,Nz)
     call getBCx3D(f(1),fsend(Nyz2+1),Nx-2,Nx,Ny,Nz)
!
     call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nyz3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR) 
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
     if (iD == iDstart) then
!
       call getBCx3D(f(1),fsend(1     ),2,Nx,Ny,Nz)
       call getBCx3D(f(1),fsend(Nyz +1),3,Nx,Ny,Nz)
       call getBCx3D(f(1),fsend(Nyz2+1),4,Nx,Ny,Nz)
!
       call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDend,120, &
                         frecvm(1),Nyz3,MPI_REAL8,iDend,120, &
                         MPI_WORLD,ISTATUS,IERR)      
     endif
!
     if (iD == iDend) then      
!
       call getBCx3D(f(1),fsend(1     ),Nx-1,Nx,Ny,Nz)
       call getBCx3D(f(1),fsend(Nyz +1),Nx-2,Nx,Ny,Nz)
       call getBCx3D(f(1),fsend(Nyz2+1),Nx-3,Nx,Ny,Nz)
!
       call MPI_SENDRECV(fsend (1),Nyz3,MPI_REAL8,iDstart,120, &
                         frecvp(1),Nyz3,MPI_REAL8,iDstart,120, &
                         MPI_WORLD,ISTATUS,IERR)
     endif
!
     fp1(1:Nyz) = frecvp(1     :Nyz )
     fp2(1:Nyz) = frecvp(Nyz +1:Nyz2)
     fp3(1:Nyz) = frecvp(Nyz2+1:Nyz3)
     fm1(1:Nyz) = frecvm(1     :Nyz )
     fm2(1:Nyz) = frecvm(Nyz +1:Nyz2)
     fm3(1:Nyz) = frecvm(Nyz2+1:Nyz3)
!    
123  continue                                          
!-------------------------------------------------------------------------------
     ii = 0
!$OMP parallel do private(ii,jk,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 k = 1, Nz
     do 30 j = 1, Ny
       ii = j+(k-1)*Ny
       jk = (j-1)*Nx+(k-1)*Nxy
!      
       do i = 3, Nx-3
         ijkp3 = (i+3)+jk
         ijkp2 = (i+2)+jk
         ijkp1 = (i+1)+jk
         ijk   = (i  )+jk
         ijkm1 = (i-1)+jk
         ijkm2 = (i-2)+jk
!         ijkm3 = (i-3)+jk
         fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
       enddo      
!-------------------------------------------------------------------------------
!     for i = 1
       ijkp3 = 4+jk
       ijkp2 = 3+jk
       ijkp1 = 2+jk
       ijk   = 1+jk
!
       fINTP(ijk) = smooth5th(fm2(ii),fm1(ii),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for i = 2
       ijkp3 = 5+jk
       ijkp2 = 4+jk
       ijkp1 = 3+jk
       ijk   = 2+jk
       ijkm1 = 1+jk     
!
       fINTP(ijk) = smooth5th(fm1(ii),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for i = 3
!       ijkp3 = 6+jk
!       ijkp2 = 5+jk
!       ijkp1 = 4+jk
!       ijk   = 3+jk
!       ijkm1 = 2+jk
!       ijkm2 = 1+jk     
!
!       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for i = Nx-2
       ijkp2 = (Nx  )+jk
       ijkp1 = (Nx-1)+jk
       ijk   = (Nx-2)+jk
       ijkm1 = (Nx-3)+jk
       ijkm2 = (Nx-4)+jk
!       ijkm3 = (Nx-5)+jk              
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),fp1(ii))
!-------------------------------------------------------------------------------
!      for i = Nx-1
       ijkp1 = (Nx  )+jk
       ijk   = (Nx-1)+jk
       ijkm1 = (Nx-2)+jk
       ijkm2 = (Nx-3)+jk
!       ijkm3 = (Nx-4)+jk              
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),fp1(ii),fp2(ii))
!-------------------------------------------------------------------------------
!      for i = Nx
       ijk   = (Nx  )+jk
       ijkm1 = (Nx-1)+jk
       ijkm2 = (Nx-2)+jk
!       ijkm3 = (Nx-3)+jk              
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),fp1(ii),fp2(ii),fp3(ii))
 30  continue
!$OMP end parallel do 
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm) 
   end subroutine PINTPxp5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fifth-order interpolation method
!    for free boundary condition.
!-------------------------------------------------------------------------------
   subroutine FINTPym5th3D_MPI(f,fINTP,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input arrays  : f(Nxyz)
!  output arrays : fINTP(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fINTP(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy  = Nx*Ny
     Nxz  = Nx*Nz
     Nxz2 = 2*Nxz
     Nxz3 = 3*Nxz
     allocate (fm3(Nxz),fm2(Nxz),fm1(Nxz),fp1(Nxz),fp2(Nxz),fp3(Nxz))
     allocate (fsend(Nxz3),frecvp(Nxz3),frecvm(Nxz3))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCy3D(f(1),fsend(     1),1,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz +1),2,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxz3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nxz3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fp1(1:Nxz) = frecvp(1     :Nxz )
     fp2(1:Nxz) = frecvp(Nxz +1:Nxz2)
     fp3(1:Nxz) = frecvp(Nxz2+1:Nxz3)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCy3D(f(1),fsend(     1),Ny  ,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz +1),Ny-1,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz2+1),Ny-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxz3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nxz3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fm1(1:Nxz) = frecvm(1     :Nxz )
     fm2(1:Nxz) = frecvm(Nxz +1:Nxz2)
     fm3(1:Nxz) = frecvm(Nxz2+1:Nxz3)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
     if (iD == iDstart) then
       call getBCy3D(f(1),fm1(1),1,Nx,Ny,Nz)      
       call getBCy3D(f(1),fm2(1),1,Nx,Ny,Nz)
       call getBCy3D(f(1),fm3(1),1,Nx,Ny,Nz)     
     endif
!
     if (iD == iDend) then           
       call getBCy3D(f(1),fp1(1),Ny,Nx,Ny,Nz)      
       call getBCy3D(f(1),fp2(1),Ny,Nx,Ny,Nz)
       call getBCy3D(f(1),fp3(1),Ny,Nx,Ny,Nz)      
     endif                                           
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,ik,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 k = 1, Nz
     do 30 i = 1, Nx
       ii = i+(k-1)*Nx
       ik = i+(k-1)*Nxy
!      
       do j = 4, Ny-2
!         ijkp3 = ((j+3)-1)*Nx+ik
         ijkp2 = ((j+2)-1)*Nx+ik
         ijkp1 = ((j+1)-1)*Nx+ik
         ijk   = ((j  )-1)*Nx+ik
         ijkm1 = ((j-1)-1)*Nx+ik
         ijkm2 = ((j-2)-1)*Nx+ik
         ijkm3 = ((j-3)-1)*Nx+ik       
!
         fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
       enddo      
!-------------------------------------------------------------------------------
!      for j = 1
!       ijkp3 = (4-1)*Nx+ik
       ijkp2 = (3-1)*Nx+ik
       ijkp1 = (2-1)*Nx+ik
       ijk   = (1-1)*Nx+ik  
!
       fINTP(ijk) = smooth5th(fm3(ii),fm2(ii),fm1(ii),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for j = 2
!       ijkp3 = (5-1)*Nx+ik
       ijkp2 = (4-1)*Nx+ik
       ijkp1 = (3-1)*Nx+ik
       ijk   = (2-1)*Nx+ik
       ijkm1 = (1-1)*Nx+ik    
!
       fINTP(ijk) = smooth5th(fm2(ii),fm1(ii),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for j = 3
!       ijkp3 = (6-1)*Nx+ik
       ijkp2 = (5-1)*Nx+ik
       ijkp1 = (4-1)*Nx+ik
       ijk   = (3-1)*Nx+ik
       ijkm1 = (2-1)*Nx+ik
       ijkm2 = (1-1)*Nx+ik    
!
       fINTP(ijk) = smooth5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for j = Ny-2
!       ijkp2 = ((Ny  )-1)*Nx+ik
!       ijkp1 = ((Ny-1)-1)*Nx+ik
!       ijk   = ((Ny-2)-1)*Nx+ik
!       ijkm1 = ((Ny-3)-1)*Nx+ik
!       ijkm2 = ((Ny-4)-1)*Nx+ik
!       ijkm3 = ((Ny-5)-1)*Nx+ik             
!
!       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for j = Ny-1
       ijkp1 = ((Ny  )-1)*Nx+ik
       ijk   = ((Ny-1)-1)*Nx+ik
       ijkm1 = ((Ny-2)-1)*Nx+ik
       ijkm2 = ((Ny-3)-1)*Nx+ik
       ijkm3 = ((Ny-4)-1)*Nx+ik             
!
       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),fp1(ii))
!-------------------------------------------------------------------------------
!      for j = Ny
       ijk   = ((Ny  )-1)*Nx+ik
       ijkm1 = ((Ny-1)-1)*Nx+ik
       ijkm2 = ((Ny-2)-1)*Nx+ik
       ijkm3 = ((Ny-3)-1)*Nx+ik             
!
       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),fp1(ii),fp2(ii))
 30  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm)
   end subroutine FINTPym5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fifth-order interpolation method
!    for periodic boundary condition.
!-------------------------------------------------------------------------------
   subroutine PINTPym5th3D_MPI(f,fINTP,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input arrays  : f(Nxyz)
!  output arrays : fINTP(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fINTP(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy  = Nx*Ny
     Nxz  = Nx*Nz
     Nxz2 = Nxz*2
     Nxz3 = Nxz*3
     allocate (fm3(Nxz),fm2(Nxz),fm1(Nxz),fp1(Nxz),fp2(Nxz),fp3(Nxz))
     allocate (fsend(Nxz3),frecvp(Nxz3),frecvm(Nxz3))
!-------------------------------------------------------------------------------
!  if Q = 1, it doesn't transfer data in MPI
!
     if (iDend == iDstart) then
       call getBCy3D(f(1),fp3(1),4,Nx,Ny,Nz)
       call getBCy3D(f(1),fp2(1),3,Nx,Ny,Nz)
       call getBCy3D(f(1),fp1(1),2,Nx,Ny,Nz)
!
       call getBCy3D(f(1),fm1(1),Ny-1,Nx,Ny,Nz)
       call getBCy3D(f(1),fm2(1),Ny-2,Nx,Ny,Nz)
       call getBCy3D(f(1),fm3(1),Ny-3,Nx,Ny,Nz)
       goto 123
     endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCy3D(f(1),fsend(1     ),1,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz +1),2,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend (1),Nxz3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nxz3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCy3D(f(1),fsend(1     ),Ny  ,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz +1),Ny-1,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz2+1),Ny-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend (1),Nxz3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nxz3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
     if (iD == iDstart) then
!
       call getBCy3D(f(1),fsend(1     ),2,Nx,Ny,Nz)
       call getBCy3D(f(1),fsend(Nxz +1),3,Nx,Ny,Nz)
       call getBCy3D(f(1),fsend(Nxz2+1),4,Nx,Ny,Nz)
       call MPI_SENDRECV(fsend (1),Nxz3,MPI_REAL8,iDend,120, &
                         frecvm(1),Nxz3,MPI_REAL8,iDend,120, &
                         MPI_WORLD,ISTATUS,IERR)
     endif
!
     if (iD == iDend) then 
!
       call getBCy3D(f(1),fsend(1     ),Ny-1,Nx,Ny,Nz)
       call getBCy3D(f(1),fsend(Nxz +1),Ny-2,Nx,Ny,Nz)
       call getBCy3D(f(1),fsend(Nxz2+1),Ny-3,Nx,Ny,Nz)
       call MPI_SENDRECV(fsend (1),Nxz3,MPI_REAL8,iDstart,120, &
                         frecvp(1),Nxz3,MPI_REAL8,iDstart,120, &
                         MPI_WORLD,ISTATUS,IERR) 
     endif
!-------------------------------------------------------------------------------
     fp1(1:Nxz) = frecvp(1     :Nxz )
     fp2(1:Nxz) = frecvp(Nxz +1:Nxz2)
     fp3(1:Nxz) = frecvp(Nxz2+1:Nxz3)
     fm1(1:Nxz) = frecvm(1     :Nxz )
     fm2(1:Nxz) = frecvm(Nxz +1:Nxz2)
     fm3(1:Nxz) = frecvm(Nxz2+1:Nxz3)
!   
123  continue                                          
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,ik,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 k = 1, Nz
     do 30 i = 1, Nx
       ii = i+(k-1)*Nx
       ik = i+(k-1)*Nxy
!      
       do j = 4, Ny-2
!         ijkp3 = ((j+3)-1)*Nx+ik
         ijkp2 = ((j+2)-1)*Nx+ik
         ijkp1 = ((j+1)-1)*Nx+ik
         ijk   = ((j  )-1)*Nx+ik
         ijkm1 = ((j-1)-1)*Nx+ik
         ijkm2 = ((j-2)-1)*Nx+ik
         ijkm3 = ((j-3)-1)*Nx+ik       
!
         fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
       enddo      
!-------------------------------------------------------------------------------
!      for j = 1
!       ijkp3 = (4-1)*Nx+ik
       ijkp2 = (3-1)*Nx+ik
       ijkp1 = (2-1)*Nx+ik
       ijk   = (1-1)*Nx+ik  
!
       fINTP(ijk) = smooth5th(fm3(ii),fm2(ii),fm1(ii),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for j = 2
!       ijkp3 = (5-1)*Nx+ik
       ijkp2 = (4-1)*Nx+ik
       ijkp1 = (3-1)*Nx+ik
       ijk   = (2-1)*Nx+ik
       ijkm1 = (1-1)*Nx+ik    
!
       fINTP(ijk) = smooth5th(fm2(ii),fm1(ii),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for j = 3
!       ijkp3 = (6-1)*Nx+ik
       ijkp2 = (5-1)*Nx+ik
       ijkp1 = (4-1)*Nx+ik
       ijk   = (3-1)*Nx+ik
       ijkm1 = (2-1)*Nx+ik
       ijkm2 = (1-1)*Nx+ik    
!
       fINTP(ijk) = smooth5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for j = Ny-2
!       ijkp2 = ((Ny  )-1)*Nx+ik
!       ijkp1 = ((Ny-1)-1)*Nx+ik
!       ijk   = ((Ny-2)-1)*Nx+ik
!       ijkm1 = ((Ny-3)-1)*Nx+ik
!       ijkm2 = ((Ny-4)-1)*Nx+ik
!       ijkm3 = ((Ny-5)-1)*Nx+ik             
!
!       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for j = Ny-1
       ijkp1 = ((Ny  )-1)*Nx+ik
       ijk   = ((Ny-1)-1)*Nx+ik
       ijkm1 = ((Ny-2)-1)*Nx+ik
       ijkm2 = ((Ny-3)-1)*Nx+ik
       ijkm3 = ((Ny-4)-1)*Nx+ik             
!
       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),fp1(ii))
!-------------------------------------------------------------------------------
!      for j = Ny
       ijk   = ((Ny  )-1)*Nx+ik
       ijkm1 = ((Ny-1)-1)*Nx+ik
       ijkm2 = ((Ny-2)-1)*Nx+ik
       ijkm3 = ((Ny-3)-1)*Nx+ik             
!
       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),fp1(ii),fp2(ii))
 30  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm) 
   end subroutine PINTPym5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fifth-order interpolation method
!    for free boundary condition.
!-------------------------------------------------------------------------------
   subroutine FINTPyp5th3D_MPI(f,fINTP,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input arrays  : f(Nxyz)
!  output arrays : fINTP(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fINTP(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy  = Nx*Ny
     Nxz  = Nx*Nz
     Nxz2 = 2*Nxz
     Nxz3 = 3*Nxz
     allocate (fm3(Nxz),fm2(Nxz),fm1(Nxz),fp1(Nxz),fp2(Nxz),fp3(Nxz))
     allocate (fsend(Nxz3),frecvp(Nxz3),frecvm(Nxz3))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCy3D(f(1),fsend(     1),1,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz +1),2,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxz3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nxz3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fp1(1:Nxz) = frecvp(1     :Nxz )
     fp2(1:Nxz) = frecvp(Nxz +1:Nxz2)
     fp3(1:Nxz) = frecvp(Nxz2+1:Nxz3)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCy3D(f(1),fsend(     1),Ny  ,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz +1),Ny-1,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz2+1),Ny-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxz3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nxz3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fm1(1:Nxz) = frecvm(1     :Nxz )
     fm2(1:Nxz) = frecvm(Nxz +1:Nxz2)
     fm3(1:Nxz) = frecvm(Nxz2+1:Nxz3)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
     if (iD == iDstart) then
       call getBCy3D(f(1),fm1(1),1,Nx,Ny,Nz)      
       call getBCy3D(f(1),fm2(1),1,Nx,Ny,Nz)
       call getBCy3D(f(1),fm3(1),1,Nx,Ny,Nz)     
     endif
!
     if (iD == iDend) then           
       call getBCy3D(f(1),fp1(1),Ny,Nx,Ny,Nz)      
       call getBCy3D(f(1),fp2(1),Ny,Nx,Ny,Nz)
       call getBCy3D(f(1),fp3(1),Ny,Nx,Ny,Nz)      
     endif                                           
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,ik,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 k = 1, Nz
     do 30 i = 1, Nx
       ii = i+(k-1)*Nx
       ik = i+(k-1)*Nxy
!      
       do j = 3, Ny-3
         ijkp3 = ((j+3)-1)*Nx+ik
         ijkp2 = ((j+2)-1)*Nx+ik
         ijkp1 = ((j+1)-1)*Nx+ik
         ijk   = ((j  )-1)*Nx+ik
         ijkm1 = ((j-1)-1)*Nx+ik
         ijkm2 = ((j-2)-1)*Nx+ik
!         ijkm3 = ((j-3)-1)*Nx+ik       
!
         fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
       enddo      
!-------------------------------------------------------------------------------
!      for j = 1
       ijkp3 = (4-1)*Nx+ik
       ijkp2 = (3-1)*Nx+ik
       ijkp1 = (2-1)*Nx+ik
       ijk   = (1-1)*Nx+ik  
!
       fINTP(ijk) = smooth5th(fm2(ii),fm1(ii),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for j = 2
       ijkp3 = (5-1)*Nx+ik
       ijkp2 = (4-1)*Nx+ik
       ijkp1 = (3-1)*Nx+ik
       ijk   = (2-1)*Nx+ik
       ijkm1 = (1-1)*Nx+ik    
!
       fINTP(ijk) = smooth5th(fm1(ii),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for j = 3
!       ijkp3 = (6-1)*Nx+ik
!       ijkp2 = (5-1)*Nx+ik
!       ijkp1 = (4-1)*Nx+ik
!       ijk   = (3-1)*Nx+ik
!       ijkm1 = (2-1)*Nx+ik
!       ijkm2 = (1-1)*Nx+ik    
!
!       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for j = Ny-2
       ijkp2 = ((Ny  )-1)*Nx+ik
       ijkp1 = ((Ny-1)-1)*Nx+ik
       ijk   = ((Ny-2)-1)*Nx+ik
       ijkm1 = ((Ny-3)-1)*Nx+ik
       ijkm2 = ((Ny-4)-1)*Nx+ik
!       ijkm3 = ((Ny-5)-1)*Nx+ik             
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),fp1(ii))
!-------------------------------------------------------------------------------
!      for j = Ny-1
       ijkp1 = ((Ny  )-1)*Nx+ik
       ijk   = ((Ny-1)-1)*Nx+ik
       ijkm1 = ((Ny-2)-1)*Nx+ik
       ijkm2 = ((Ny-3)-1)*Nx+ik
!       ijkm3 = ((Ny-4)-1)*Nx+ik             
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),fp1(ii),fp2(ii))
!-------------------------------------------------------------------------------
!      for j = Ny
       ijk   = ((Ny  )-1)*Nx+ik
       ijkm1 = ((Ny-1)-1)*Nx+ik
       ijkm2 = ((Ny-2)-1)*Nx+ik
!       ijkm3 = ((Ny-3)-1)*Nx+ik             
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),fp1(ii),fp2(ii),fp3(ii))
 30  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm)
   end subroutine FINTPyp5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fifth-order interpolation method
!    for periodic boundary condition.
!-------------------------------------------------------------------------------
   subroutine PINTPyp5th3D_MPI(f,fINTP,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input arrays  : f(Nxyz)
!  output arrays : fINTP(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fINTP(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy  = Nx*Ny
     Nxz  = Nx*Nz
     Nxz2 = Nxz*2
     Nxz3 = Nxz*3
     allocate (fm3(Nxz),fm2(Nxz),fm1(Nxz),fp1(Nxz),fp2(Nxz),fp3(Nxz))
     allocate (fsend(Nxz3),frecvp(Nxz3),frecvm(Nxz3))
!-------------------------------------------------------------------------------
!  if Q = 1, it doesn't transfer data in MPI
!
     if (iDend == iDstart) then
       call getBCy3D(f(1),fp3(1),4,Nx,Ny,Nz)
       call getBCy3D(f(1),fp2(1),3,Nx,Ny,Nz)
       call getBCy3D(f(1),fp1(1),2,Nx,Ny,Nz)
!
       call getBCy3D(f(1),fm1(1),Ny-1,Nx,Ny,Nz)
       call getBCy3D(f(1),fm2(1),Ny-2,Nx,Ny,Nz)
       call getBCy3D(f(1),fm3(1),Ny-3,Nx,Ny,Nz)
       goto 123
     endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCy3D(f(1),fsend(1     ),1,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz +1),2,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend (1),Nxz3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nxz3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCy3D(f(1),fsend(1     ),Ny  ,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz +1),Ny-1,Nx,Ny,Nz)
     call getBCy3D(f(1),fsend(Nxz2+1),Ny-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend (1),Nxz3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nxz3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
     if (iD == iDstart) then
!
       call getBCy3D(f(1),fsend(1     ),2,Nx,Ny,Nz)
       call getBCy3D(f(1),fsend(Nxz +1),3,Nx,Ny,Nz)
       call getBCy3D(f(1),fsend(Nxz2+1),4,Nx,Ny,Nz)
       call MPI_SENDRECV(fsend (1),Nxz3,MPI_REAL8,iDend,120, &
                         frecvm(1),Nxz3,MPI_REAL8,iDend,120, &
                         MPI_WORLD,ISTATUS,IERR)
     endif
!
     if (iD == iDend) then 
!
       call getBCy3D(f(1),fsend(1     ),Ny-1,Nx,Ny,Nz)
       call getBCy3D(f(1),fsend(Nxz +1),Ny-2,Nx,Ny,Nz)
       call getBCy3D(f(1),fsend(Nxz2+1),Ny-3,Nx,Ny,Nz)
       call MPI_SENDRECV(fsend (1),Nxz3,MPI_REAL8,iDstart,120, &
                         frecvp(1),Nxz3,MPI_REAL8,iDstart,120, &
                         MPI_WORLD,ISTATUS,IERR) 
     endif
!-------------------------------------------------------------------------------
     fp1(1:Nxz) = frecvp(1     :Nxz )
     fp2(1:Nxz) = frecvp(Nxz +1:Nxz2)
     fp3(1:Nxz) = frecvp(Nxz2+1:Nxz3)
     fm1(1:Nxz) = frecvm(1     :Nxz )
     fm2(1:Nxz) = frecvm(Nxz +1:Nxz2)
     fm3(1:Nxz) = frecvm(Nxz2+1:Nxz3)
!   
123  continue                                          
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,ik,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 k = 1, Nz
     do 30 i = 1, Nx
       ii = i+(k-1)*Nx
       ik = i+(k-1)*Nxy
!      
       do j = 3, Ny-3
         ijkp3 = ((j+3)-1)*Nx+ik
         ijkp2 = ((j+2)-1)*Nx+ik
         ijkp1 = ((j+1)-1)*Nx+ik
         ijk   = ((j  )-1)*Nx+ik
         ijkm1 = ((j-1)-1)*Nx+ik
         ijkm2 = ((j-2)-1)*Nx+ik
!         ijkm3 = ((j-3)-1)*Nx+ik       
!
         fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
       enddo      
!-------------------------------------------------------------------------------
!      for j = 1
       ijkp3 = (4-1)*Nx+ik
       ijkp2 = (3-1)*Nx+ik
       ijkp1 = (2-1)*Nx+ik
       ijk   = (1-1)*Nx+ik  
!
       fINTP(ijk) = smooth5th(fm2(ii),fm1(ii),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for j = 2
       ijkp3 = (5-1)*Nx+ik
       ijkp2 = (4-1)*Nx+ik
       ijkp1 = (3-1)*Nx+ik
       ijk   = (2-1)*Nx+ik
       ijkm1 = (1-1)*Nx+ik    
!
       fINTP(ijk) = smooth5th(fm1(ii),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for j = 3
!       ijkp3 = (6-1)*Nx+ik
!       ijkp2 = (5-1)*Nx+ik
!       ijkp1 = (4-1)*Nx+ik
!       ijk   = (3-1)*Nx+ik
!       ijkm1 = (2-1)*Nx+ik
!       ijkm2 = (1-1)*Nx+ik    
!
!       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for j = Ny-2
       ijkp2 = ((Ny  )-1)*Nx+ik
       ijkp1 = ((Ny-1)-1)*Nx+ik
       ijk   = ((Ny-2)-1)*Nx+ik
       ijkm1 = ((Ny-3)-1)*Nx+ik
       ijkm2 = ((Ny-4)-1)*Nx+ik
!       ijkm3 = ((Ny-5)-1)*Nx+ik             
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),fp1(ii))
!-------------------------------------------------------------------------------
!      for j = Ny-1
       ijkp1 = ((Ny  )-1)*Nx+ik
       ijk   = ((Ny-1)-1)*Nx+ik
       ijkm1 = ((Ny-2)-1)*Nx+ik
       ijkm2 = ((Ny-3)-1)*Nx+ik
!       ijkm3 = ((Ny-4)-1)*Nx+ik             
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),fp1(ii),fp2(ii))
!-------------------------------------------------------------------------------
!      for j = Ny
       ijk   = ((Ny  )-1)*Nx+ik
       ijkm1 = ((Ny-1)-1)*Nx+ik
       ijkm2 = ((Ny-2)-1)*Nx+ik
!       ijkm3 = ((Ny-3)-1)*Nx+ik             
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),fp1(ii),fp2(ii),fp3(ii))
 30  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm) 
   end subroutine PINTPyp5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of the fifth-order interpolation method
!    for free boundary condition.
!-------------------------------------------------------------------------------
   subroutine FINTPzm5th3D_MPI(f,fINTP,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input arrays  : f(Nxyz)
!  output arrays : fINTP(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fINTP(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy  = Nx*Ny
     Nxy2 = 2*Nxy
     Nxy3 = 3*Nxy
     allocate (fm3(Nxy),fm2(Nxy),fm1(Nxy),fp1(Nxy),fp2(Nxy),fp3(Nxy))
     allocate (fsend(Nxy3),frecvp(Nxy3),frecvm(Nxy3))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCz3D(f(1),fsend(     1),1,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy +1),2,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nxy3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fp1(1:Nxy) = frecvp(1     :Nxy )
     fp2(1:Nxy) = frecvp(Nxy +1:Nxy2)
     fp3(1:Nxy) = frecvp(Nxy2+1:Nxy3)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCz3D(f(1),fsend(     1),Nz  ,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy +1),Nz-1,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy2+1),Nz-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nxy3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fm1(1:Nxy) = frecvm(1     :Nxy )
     fm2(1:Nxy) = frecvm(Nxy +1:Nxy2)
     fm3(1:Nxy) = frecvm(Nxy2+1:Nxy3)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
     if (iD == iDstart) then
       call getBCz3D(f(1),fm1(1),1,Nx,Ny,Nz)      
       call getBCz3D(f(1),fm2(1),1,Nx,Ny,Nz)
       call getBCz3D(f(1),fm3(1),1,Nx,Ny,Nz)      
     endif
!
     if (iD == iDend) then           
       call getBCz3D(f(1),fp1(1),Nz,Nx,Ny,Nz)      
       call getBCz3D(f(1),fp2(1),Nz,Nx,Ny,Nz)
       call getBCz3D(f(1),fp3(1),Nz,Nx,Ny,Nz)      
     endif                                           
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,ij,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 j = 1, Ny
     do 30 i = 1, Nx
       ii = i+(j-1)*Nx
       ij = ii !i+(j-1)*Nx
!      
       do k = 4, Nz-2
!         ijkp3 = ((k+3)-1)*Nxy+ij
         ijkp2 = ((k+2)-1)*Nxy+ij
         ijkp1 = ((k+1)-1)*Nxy+ij
         ijk   = ((k  )-1)*Nxy+ij
         ijkm1 = ((k-1)-1)*Nxy+ij
         ijkm2 = ((k-2)-1)*Nxy+ij
         ijkm3 = ((k-3)-1)*Nxy+ij       
!
         fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
       enddo      
!-------------------------------------------------------------------------------
!      for k = 1
!       ijkp3 = (4-1)*Nxy+ij
       ijkp2 = (3-1)*Nxy+ij
       ijkp1 = (2-1)*Nxy+ij
       ijk   = (1-1)*Nxy+ij  
!
       fINTP(ijk) = smooth5th(fm3(ii),fm2(ii),fm1(ii),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for k = 2
!       ijkp3 = (5-1)*Nxy+ij
       ijkp2 = (4-1)*Nxy+ij
       ijkp1 = (3-1)*Nxy+ij
       ijk   = (2-1)*Nxy+ij
       ijkm1 = (1-1)*Nxy+ij    
!
       fINTP(ijk) = smooth5th(fm2(ii),fm1(ii),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for k = 3
!       ijkp3 = (6-1)*Nxy+ij
       ijkp2 = (5-1)*Nxy+ij
       ijkp1 = (4-1)*Nxy+ij
       ijk   = (3-1)*Nxy+ij
       ijkm1 = (2-1)*Nxy+ij 
       ijkm2 = (1-1)*Nxy+ij   
!
       fINTP(ijk) = smooth5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for k = Nz-2
!       ijkp2 = ((Nz  )-1)*Nxy+ij
!       ijkp1 = ((Nz-1)-1)*Nxy+ij
!       ijk   = ((Nz-2)-1)*Nxy+ij
!       ijkm1 = ((Nz-3)-1)*Nxy+ij
!       ijkm2 = ((Nz-4)-1)*Nxy+ij
!       ijkm3 = ((Nz-5)-1)*Nxy+ij             
!
!       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for k = Nz-1
       ijkp1 = ((Nz  )-1)*Nxy+ij
       ijk   = ((Nz-1)-1)*Nxy+ij
       ijkm1 = ((Nz-2)-1)*Nxy+ij
       ijkm2 = ((Nz-3)-1)*Nxy+ij
       ijkm3 = ((Nz-4)-1)*Nxy+ij             
!
       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),fp1(ii))
!-------------------------------------------------------------------------------
!      for k = Nz
       ijk   = ((Nz  )-1)*Nxy+ij
       ijkm1 = ((Nz-1)-1)*Nxy+ij
       ijkm2 = ((Nz-2)-1)*Nxy+ij
       ijkm3 = ((Nz-3)-1)*Nxy+ij             
!
       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),fp1(ii),fp2(ii)) 
 30  continue
!$OMP end parallel do     
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm) 
   end subroutine FINTPzm5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of the fifth-order interpolation method
!    for periodic boundary condition.
!-------------------------------------------------------------------------------
   subroutine PINTPzm5th3D_MPI(f,fINTP,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input arrays  : f(Nxyz)
!  output arrays : fINTP(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fINTP(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm 
     Nxy  = Nx*Ny
     Nxy2 = 2*Nxy
     Nxy3 = 3*Nxy
     allocate (fm3(Nxy),fm2(Nxy),fm1(Nxy),fp1(Nxy),fp2(Nxy),fp3(Nxy))
     allocate (fsend(Nxy3),frecvp(Nxy3),frecvm(Nxy3))
!-------------------------------------------------------------------------------
!  if R = 1, it doesn't transfer data in MPI
!
     if (iDend == iDstart) then
       call getBCz3D(f(1),fp1(1),2,Nx,Ny,Nz)
       call getBCz3D(f(1),fp2(1),3,Nx,Ny,Nz)
       call getBCz3D(f(1),fp3(1),4,Nx,Ny,Nz)
!       
       call getBCz3D(f(1),fm1(1),Nz-1,Nx,Ny,Nz)
       call getBCz3D(f(1),fm2(1),Nz-2,Nx,Ny,Nz)
       call getBCz3D(f(1),fm3(1),Nz-3,Nx,Ny,Nz)
       goto 123
     endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCz3D(f(1),fsend(     1),1,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy +1),2,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nxy3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCz3D(f(1),fsend(     1),Nz  ,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy +1),Nz-1,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy2+1),Nz-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nxy3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
     if (iD == iDstart) then
!
       call getBCz3D(f(1),fsend(     1),2,Nx,Ny,Nz)
       call getBCz3D(f(1),fsend(Nxy +1),3,Nx,Ny,Nz)
       call getBCz3D(f(1),fsend(Nxy2+1),4,Nx,Ny,Nz)
       call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDend,120, &
                         frecvm(1),Nxy3,MPI_REAL8,iDend,120, &
                         MPI_WORLD,ISTATUS,IERR)
     endif
!
     if (iD == iDend) then      
!
       call getBCz3D(f(1),fsend(     1),Nz-1,Nx,Ny,Nz)
       call getBCz3D(f(1),fsend(Nxy +1),Nz-2,Nx,Ny,Nz)
       call getBCz3D(f(1),fsend(Nxy2+1),Nz-3,Nx,Ny,Nz)
       call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDstart,120, &
                         frecvp(1),Nxy3,MPI_REAL8,iDstart,120, &
                         MPI_WORLD,ISTATUS,IERR) 
     endif
!-------------------------------------------------------------------------------
     fp1(1:Nxy) = frecvp(1     :Nxy )
     fp2(1:Nxy) = frecvp(Nxy +1:Nxy2)
     fp3(1:Nxy) = frecvp(Nxy2+1:Nxy3)
     fm1(1:Nxy) = frecvm(1     :Nxy )
     fm2(1:Nxy) = frecvm(Nxy +1:Nxy2)
     fm3(1:Nxy) = frecvm(Nxy2+1:Nxy3)
!    
123  continue                                          
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,ij,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 j = 1, Ny
     do 30 i = 1, Nx
       ii = i+(j-1)*Nx
       ij = ii !i+(j-1)*Nx
!      
       do k = 4, Nz-2
!         ijkp3 = ((k+3)-1)*Nxy+ij
         ijkp2 = ((k+2)-1)*Nxy+ij
         ijkp1 = ((k+1)-1)*Nxy+ij
         ijk   = ((k  )-1)*Nxy+ij
         ijkm1 = ((k-1)-1)*Nxy+ij
         ijkm2 = ((k-2)-1)*Nxy+ij
         ijkm3 = ((k-3)-1)*Nxy+ij       
!
         fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
       enddo
!-------------------------------------------------------------------------------
!      for k = 1
!       ijkp3 = (4-1)*Nxy+ij
       ijkp2 = (3-1)*Nxy+ij
       ijkp1 = (2-1)*Nxy+ij
       ijk   = (1-1)*Nxy+ij  
!
       fINTP(ijk) = smooth5th(fm3(ii),fm2(ii),fm1(ii),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for k = 2
!       ijkp3 = (5-1)*Nxy+ij
       ijkp2 = (4-1)*Nxy+ij
       ijkp1 = (3-1)*Nxy+ij
       ijk   = (2-1)*Nxy+ij
       ijkm1 = (1-1)*Nxy+ij    
!
       fINTP(ijk) = smooth5th(fm2(ii),fm1(ii),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for k = 3
       ijkp3 = (6-1)*Nxy+ij
       ijkp2 = (5-1)*Nxy+ij
       ijkp1 = (4-1)*Nxy+ij
       ijk   = (3-1)*Nxy+ij
       ijkm1 = (2-1)*Nxy+ij 
       ijkm2 = (1-1)*Nxy+ij   
!
       fINTP(ijk) = smooth5th(fm1(ii),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for k = Nz-2
!       ijkp2 = ((Nz  )-1)*Nxy+ij
!       ijkp1 = ((Nz-1)-1)*Nxy+ij
!       ijk   = ((Nz-2)-1)*Nxy+ij
!       ijkm1 = ((Nz-3)-1)*Nxy+ij
!       ijkm2 = ((Nz-4)-1)*Nxy+ij
!       ijkm3 = ((Nz-5)-1)*Nxy+ij             
!
!       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2))
!-------------------------------------------------------------------------------
!      for k = Nz-1
       ijkp1 = ((Nz  )-1)*Nxy+ij
       ijk   = ((Nz-1)-1)*Nxy+ij
       ijkm1 = ((Nz-2)-1)*Nxy+ij
       ijkm2 = ((Nz-3)-1)*Nxy+ij
       ijkm3 = ((Nz-4)-1)*Nxy+ij             
!
       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),fp1(ii))
!-------------------------------------------------------------------------------
!      for k = Nz
       ijk   = ((Nz  )-1)*Nxy+ij
       ijkm1 = ((Nz-1)-1)*Nxy+ij
       ijkm2 = ((Nz-2)-1)*Nxy+ij
       ijkm3 = ((Nz-3)-1)*Nxy+ij             
!
       fINTP(ijk) = smooth5th(f(ijkm3),f(ijkm2),f(ijkm1),f(ijk),fp1(ii),fp2(ii)) 
 30  continue
!$OMP end parallel do     
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm)  
   end subroutine PINTPzm5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of the fifth-order interpolation method
!    for free boundary condition.
!-------------------------------------------------------------------------------
   subroutine FINTPzp5th3D_MPI(f,fINTP,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input arrays  : f(Nxyz)
!  output arrays : fINTP(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fINTP(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3 
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm
     Nxy  = Nx*Ny
     Nxy2 = 2*Nxy
     Nxy3 = 3*Nxy
     allocate (fm3(Nxy),fm2(Nxy),fm1(Nxy),fp1(Nxy),fp2(Nxy),fp3(Nxy))
     allocate (fsend(Nxy3),frecvp(Nxy3),frecvm(Nxy3))
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCz3D(f(1),fsend(     1),1,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy +1),2,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nxy3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fp1(1:Nxy) = frecvp(1     :Nxy )
     fp2(1:Nxy) = frecvp(Nxy +1:Nxy2)
     fp3(1:Nxy) = frecvp(Nxy2+1:Nxy3)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCz3D(f(1),fsend(     1),Nz  ,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy +1),Nz-1,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy2+1),Nz-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nxy3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!
     fm1(1:Nxy) = frecvm(1     :Nxy )
     fm2(1:Nxy) = frecvm(Nxy +1:Nxy2)
     fm3(1:Nxy) = frecvm(Nxy2+1:Nxy3)
!-------------------------------------------------------------------------------
!  set the data for free boundary condition
!
     if (iD == iDstart) then
       call getBCz3D(f(1),fm1(1),1,Nx,Ny,Nz)      
       call getBCz3D(f(1),fm2(1),1,Nx,Ny,Nz)
       call getBCz3D(f(1),fm3(1),1,Nx,Ny,Nz)      
     endif
!
     if (iD == iDend) then           
       call getBCz3D(f(1),fp1(1),Nz,Nx,Ny,Nz)      
       call getBCz3D(f(1),fp2(1),Nz,Nx,Ny,Nz)
       call getBCz3D(f(1),fp3(1),Nz,Nx,Ny,Nz)      
     endif                                           
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,ij,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 j = 1, Ny
     do 30 i = 1, Nx
       ii = i+(j-1)*Nx
       ij = ii !i+(j-1)*Nx
!      
       do k = 3, Nz-3
         ijkp3 = ((k+3)-1)*Nxy+ij
         ijkp2 = ((k+2)-1)*Nxy+ij
         ijkp1 = ((k+1)-1)*Nxy+ij
         ijk   = ((k  )-1)*Nxy+ij
         ijkm1 = ((k-1)-1)*Nxy+ij
         ijkm2 = ((k-2)-1)*Nxy+ij
!         ijkm3 = ((k-3)-1)*Nxy+ij       
!
         fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
       enddo      
!-------------------------------------------------------------------------------
!      for k = 1
       ijkp3 = (4-1)*Nxy+ij
       ijkp2 = (3-1)*Nxy+ij
       ijkp1 = (2-1)*Nxy+ij
       ijk   = (1-1)*Nxy+ij  
!
       fINTP(ijk) = smooth5th(fm2(ii),fm1(ii),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for k = 2
       ijkp3 = (5-1)*Nxy+ij
       ijkp2 = (4-1)*Nxy+ij
       ijkp1 = (3-1)*Nxy+ij
       ijk   = (2-1)*Nxy+ij
       ijkm1 = (1-1)*Nxy+ij    
!
       fINTP(ijk) = smooth5th(fm1(ii),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for k = 3
!       ijkp3 = (6-1)*Nxy+ij
!       ijkp2 = (5-1)*Nxy+ij
!       ijkp1 = (4-1)*Nxy+ij
!       ijk   = (3-1)*Nxy+ij
!       ijkm1 = (2-1)*Nxy+ij 
!       ijkm2 = (1-1)*Nxy+ij   
!
!       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for k = Nz-2
       ijkp2 = ((Nz  )-1)*Nxy+ij
       ijkp1 = ((Nz-1)-1)*Nxy+ij
       ijk   = ((Nz-2)-1)*Nxy+ij
       ijkm1 = ((Nz-3)-1)*Nxy+ij
       ijkm2 = ((Nz-4)-1)*Nxy+ij
!       ijkm3 = ((Nz-5)-1)*Nxy+ij             
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),fp1(ii))
!-------------------------------------------------------------------------------
!      for k = Nz-1
       ijkp1 = ((Nz  )-1)*Nxy+ij
       ijk   = ((Nz-1)-1)*Nxy+ij
       ijkm1 = ((Nz-2)-1)*Nxy+ij
       ijkm2 = ((Nz-3)-1)*Nxy+ij
!       ijkm3 = ((Nz-4)-1)*Nxy+ij             
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),fp1(ii),fp2(ii))
!-------------------------------------------------------------------------------
!      for k = Nz
       ijk   = ((Nz  )-1)*Nxy+ij
       ijkm1 = ((Nz-1)-1)*Nxy+ij
       ijkm2 = ((Nz-2)-1)*Nxy+ij
!       ijkm3 = ((Nz-3)-1)*Nxy+ij             
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),fp1(ii),fp2(ii),fp3(ii)) 
 30  continue
!$OMP end parallel do     
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm) 
   end subroutine FINTPzp5th3D_MPI
!===============================================================================
!===============================================================================
!  This is the step of the fifth-order interpolation method
!    for periodic boundary condition.
!-------------------------------------------------------------------------------
   subroutine PINTPzp5th3D_MPI(f,fINTP,Nx,Ny,Nz,iDm1,iD,iDp1,iDstart,iDend)
!-------------------------------------------------------------------------------
!  input arrays  : f(Nxyz)
!  output arrays : fINTP(Nxyz)
!-------------------------------------------------------------------------------
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fINTP(1)
     real*8, allocatable, dimension(:) :: fm3, fm2, fm1, fp1, fp2, fp3
     real*8, allocatable, dimension(:) :: fsend, frecvp, frecvm 
     Nxy  = Nx*Ny
     Nxy2 = 2*Nxy
     Nxy3 = 3*Nxy
     allocate (fm3(Nxy),fm2(Nxy),fm1(Nxy),fp1(Nxy),fp2(Nxy),fp3(Nxy))
     allocate (fsend(Nxy3),frecvp(Nxy3),frecvm(Nxy3))
!-------------------------------------------------------------------------------
!  if R = 1, it doesn't transfer data in MPI
!
     if (iDend == iDstart) then
       call getBCz3D(f(1),fp1(1),2,Nx,Ny,Nz)
       call getBCz3D(f(1),fp2(1),3,Nx,Ny,Nz)
       call getBCz3D(f(1),fp3(1),4,Nx,Ny,Nz)
!       
       call getBCz3D(f(1),fm1(1),Nz-1,Nx,Ny,Nz)
       call getBCz3D(f(1),fm2(1),Nz-2,Nx,Ny,Nz)
       call getBCz3D(f(1),fm3(1),Nz-3,Nx,Ny,Nz)
       goto 123
     endif
!-------------------------------------------------------------------------------
!  send the lower boundary to the previous node
!
     call getBCz3D(f(1),fsend(     1),1,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy +1),2,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy2+1),3,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDm1,100, &
                       frecvp(1),Nxy3,MPI_REAL8,iDp1,100, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send the upper boundary to the next node
!
     call getBCz3D(f(1),fsend(     1),Nz  ,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy +1),Nz-1,Nx,Ny,Nz)
     call getBCz3D(f(1),fsend(Nxy2+1),Nz-2,Nx,Ny,Nz)
     call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDp1,110, &
                       frecvm(1),Nxy3,MPI_REAL8,iDm1,110, &
                       MPI_WORLD,ISTATUS,IERR)
!-------------------------------------------------------------------------------
!  send and receive data for periodic boundary condition
!
     if (iD == iDstart) then
!
       call getBCz3D(f(1),fsend(     1),2,Nx,Ny,Nz)
       call getBCz3D(f(1),fsend(Nxy +1),3,Nx,Ny,Nz)
       call getBCz3D(f(1),fsend(Nxy2+1),4,Nx,Ny,Nz)
       call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDend,120, &
                         frecvm(1),Nxy3,MPI_REAL8,iDend,120, &
                         MPI_WORLD,ISTATUS,IERR)
     endif
!
     if (iD == iDend) then      
!
       call getBCz3D(f(1),fsend(     1),Nz-1,Nx,Ny,Nz)
       call getBCz3D(f(1),fsend(Nxy +1),Nz-2,Nx,Ny,Nz)
       call getBCz3D(f(1),fsend(Nxy2+1),Nz-3,Nx,Ny,Nz)
       call MPI_SENDRECV(fsend(1) ,Nxy3,MPI_REAL8,iDstart,120, &
                         frecvp(1),Nxy3,MPI_REAL8,iDstart,120, &
                         MPI_WORLD,ISTATUS,IERR) 
     endif
!-------------------------------------------------------------------------------
     fp1(1:Nxy) = frecvp(1     :Nxy )
     fp2(1:Nxy) = frecvp(Nxy +1:Nxy2)
     fp3(1:Nxy) = frecvp(Nxy2+1:Nxy3)
     fm1(1:Nxy) = frecvm(1     :Nxy )
     fm2(1:Nxy) = frecvm(Nxy +1:Nxy2)
     fm3(1:Nxy) = frecvm(Nxy2+1:Nxy3)
!    
123  continue                                          
!-------------------------------------------------------------------------------
!$OMP parallel do private(ii,ij,ijkp3,ijkp2,ijkp1,ijk,ijkm1,ijkm2,ijkm3) collapse(2)
     do 30 j = 1, Ny
     do 30 i = 1, Nx
       ii = i+(j-1)*Nx
       ij = ii !i+(j-1)*Nx
!      
       do k = 3, Nz-3
         ijkp3 = ((k+3)-1)*Nxy+ij
         ijkp2 = ((k+2)-1)*Nxy+ij
         ijkp1 = ((k+1)-1)*Nxy+ij
         ijk   = ((k  )-1)*Nxy+ij
         ijkm1 = ((k-1)-1)*Nxy+ij
         ijkm2 = ((k-2)-1)*Nxy+ij
!         ijkm3 = ((k-3)-1)*Nxy+ij       
!
         fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
       enddo
!-------------------------------------------------------------------------------
!      for k = 1
       ijkp3 = (4-1)*Nxy+ij
       ijkp2 = (3-1)*Nxy+ij
       ijkp1 = (2-1)*Nxy+ij
       ijk   = (1-1)*Nxy+ij  
!
       fINTP(ijk) = smooth5th(fm2(ii),fm1(ii),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for k = 2
       ijkp3 = (5-1)*Nxy+ij
       ijkp2 = (4-1)*Nxy+ij
       ijkp1 = (3-1)*Nxy+ij
       ijk   = (2-1)*Nxy+ij
       ijkm1 = (1-1)*Nxy+ij    
!
       fINTP(ijk) = smooth5th(fm1(ii),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for k = 3
!       ijkp3 = (6-1)*Nxy+ij
!       ijkp2 = (5-1)*Nxy+ij
!       ijkp1 = (4-1)*Nxy+ij
!       ijk   = (3-1)*Nxy+ij
!       ijkm1 = (2-1)*Nxy+ij 
!       ijkm2 = (1-1)*Nxy+ij   
!
!       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),f(ijkp3))
!-------------------------------------------------------------------------------
!      for k = Nz-2
       ijkp2 = ((Nz  )-1)*Nxy+ij
       ijkp1 = ((Nz-1)-1)*Nxy+ij
       ijk   = ((Nz-2)-1)*Nxy+ij
       ijkm1 = ((Nz-3)-1)*Nxy+ij
       ijkm2 = ((Nz-4)-1)*Nxy+ij
!       ijkm3 = ((Nz-5)-1)*Nxy+ij             
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),f(ijkp2),fp1(ii))
!-------------------------------------------------------------------------------
!      for k = Nz-1
       ijkp1 = ((Nz  )-1)*Nxy+ij
       ijk   = ((Nz-1)-1)*Nxy+ij
       ijkm1 = ((Nz-2)-1)*Nxy+ij
       ijkm2 = ((Nz-3)-1)*Nxy+ij
!       ijkm3 = ((Nz-4)-1)*Nxy+ij             
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),f(ijkp1),fp1(ii),fp2(ii))
!-------------------------------------------------------------------------------
!      for k = Nz
       ijk   = ((Nz  )-1)*Nxy+ij
       ijkm1 = ((Nz-1)-1)*Nxy+ij
       ijkm2 = ((Nz-2)-1)*Nxy+ij
!       ijkm3 = ((Nz-3)-1)*Nxy+ij             
!
       fINTP(ijk) = smooth5th(f(ijkm2),f(ijkm1),f(ijk),fp1(ii),fp2(ii),fp3(ii)) 
 30  continue
!$OMP end parallel do     
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm)  
   end subroutine PINTPzp5th3D_MPI
!===============================================================================
!===============================================================================
   function smooth5th(fm3,fm2,fm1,f,fp1,fp2)
     implicit double precision (A-H,O-Z)
     a150 = 150.d0
     a25  =  25.d0
     a3   =   3.d0
     a256 = 256.d0
!
     sumA = sumABC(a150*f  ,-a25*fp1,a3*fp2)
     sumB = sumABC(a150*fm1,-a25*fm2,a3*fm3)
     smooth5th = sumAB(sumA,sumB)/a256
   end function smooth5th
!===============================================================================
!===============================================================================
end module CenDif3D_MPI
!*******************************************************************************