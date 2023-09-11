module WENO3D_MPI
      use MPI3D
      use coef
      use omp_lib
    contains
!===============================================================================
!  This is the step of fifth-order WENO Scheme
!    to calculate the first derivative fp and second derivative fpp
!    along x-direction for free boundary condition.
!-------------------------------------------------------------------------------
    subroutine WENOx5th6D_MPI(f,fp,dx,Nx,Ny,Nz,Nux,Nuy,Nuz,iDm1,iD,iDp1,iDstart,iDend)
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
      allocate (fm3(NN),fm2(NN),fm1(NN),fp1(NN),fp2(NN),fp3(NN))
      allocate (fsend(NN3),frecvp(NN3),frecvm(NN3))
!-------------------------------------------------------------------------------
!  if P = 1, it doesnt transfer data in MPI
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
!$OMP parallel do private(iim3,iim2,iim1,ii,iip1,iip2,iip3,ijku,jk1,jk2,ii1,v) collapse(3) 
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
         v=iu-(Nux+1)/2.d0
!        
         do i = 4, Nx-3
           iip3 = ijku+((i+3)+jk2-1)*Nuxyz
           iip2 = ijku+((i+2)+jk2-1)*Nuxyz
           iip1 = ijku+((i+1)+jk2-1)*Nuxyz
           ii   = ijku+((i  )+jk2-1)*Nuxyz
           iim1 = ijku+((i-1)+jk2-1)*Nuxyz
           iim2 = ijku+((i-2)+jk2-1)*Nuxyz
           iim3 = ijku+((i-3)+jk2-1)*Nuxyz
           fp(ii)  = WENO5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx,v)
         enddo      
!-------------------------------------------------------------------------------
!  for i = 1
         iip3 = ijku+(4+jk2-1)*Nuxyz
         iip2 = ijku+(3+jk2-1)*Nuxyz
         iip1 = ijku+(2+jk2-1)*Nuxyz
         ii   = ijku+(1+jk2-1)*Nuxyz
!
         fp(ii)  = WENO5th(fm3(ii1),fm2(ii1),fm1(ii1),f(ii),f(iip1),f(iip2),f(iip3),dx,v)
         !fpp(ii) = FSD5th(fm3(ii1),fm2(ii1),fm1(ii1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = 2
         iip3 = ijku+(5+jk2-1)*Nuxyz
         iip2 = ijku+(4+jk2-1)*Nuxyz
         iip1 = ijku+(3+jk2-1)*Nuxyz
         ii   = ijku+(2+jk2-1)*Nuxyz
         iim1 = ijku+(1+jk2-1)*Nuxyz     
!
         fp(ii)  = WENO5th(fm2(ii1),fm1(ii1),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx,v)
         !fpp(ii) = FSD5th(fm2(ii1),fm1(ii1),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = 3
         iip3 = ijku+(6+jk2-1)*Nuxyz
         iip2 = ijku+(5+jk2-1)*Nuxyz
         iip1 = ijku+(4+jk2-1)*Nuxyz
         ii   = ijku+(3+jk2-1)*Nuxyz
         iim1 = ijku+(2+jk2-1)*Nuxyz
         iim2 = ijku+(1+jk2-1)*Nuxyz     
!
         fp(ii)  = WENO5th(fm1(ii1),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx,v)
         !fpp(ii) = FSD5th(fm1(ii1),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx-2
         iip2 = ijku+((Nx  )+jk2-1)*Nuxyz
         iip1 = ijku+((Nx-1)+jk2-1)*Nuxyz
         ii   = ijku+((Nx-2)+jk2-1)*Nuxyz
         iim1 = ijku+((Nx-3)+jk2-1)*Nuxyz
         iim2 = ijku+((Nx-4)+jk2-1)*Nuxyz
         iim3 = ijku+((Nx-5)+jk2-1)*Nuxyz              
!
         fp(ii)  =WENO5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),fp1(ii1),dx,v)
         !fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),fp1(ii1),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx-1
         iip1 = ijku+((Nx  )+jk2-1)*Nuxyz
         ii   = ijku+((Nx-1)+jk2-1)*Nuxyz
         iim1 = ijku+((Nx-2)+jk2-1)*Nuxyz
         iim2 = ijku+((Nx-3)+jk2-1)*Nuxyz
         iim3 = ijku+((Nx-4)+jk2-1)*Nuxyz              
!
         fp(ii)  = WENO5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),fp1(ii1),fp2(ii1),dx,v)
         !fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),fp1(ii1),fp2(ii1),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx
         ii   = ijku+((Nx  )+jk2-1)*Nuxyz
         iim1 = ijku+((Nx-1)+jk2-1)*Nuxyz
         iim2 = ijku+((Nx-2)+jk2-1)*Nuxyz
         iim3 = ijku+((Nx-3)+jk2-1)*Nuxyz              
!
         fp(ii)  = WENO5th(f(iim3),f(iim2),f(iim1),f(ii),fp1(ii1),fp2(ii1),fp3(ii1),dx,v)
         !fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),fp1(ii1),fp2(ii1),fp3(ii1),dx2)
 30    continue
!
 35  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm)  
   end subroutine WENOx5th6D_MPI
!===============================================================================
!===============================================================================
!  This is the step of fourth-order C-WENO Scheme
!    to calculate the first derivative fp and second derivative fpp
!    along x-direction for free boundary condition.
!-------------------------------------------------------------------------------
    subroutine CWENOx6th6D_MPI(f,fp,dx,Nx,Ny,Nz,Nux,Nuy,Nuz,iDm1,iD,iDp1,iDstart,iDend)
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
      allocate (fm3(NN),fm2(NN),fm1(NN),fp1(NN),fp2(NN),fp3(NN))
      allocate (fsend(NN3),frecvp(NN3),frecvm(NN3))
!-------------------------------------------------------------------------------
!  if P = 1, it doesnt transfer data in MPI
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
!$OMP parallel do private(iim3,iim2,iim1,ii,iip1,iip2,iip3,ijku,jk1,jk2,ii1,v) collapse(3) 
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
         v=iu-(Nux+1)/2.d0
!        
         do i = 4, Nx-3
           iip3 = ijku+((i+3)+jk2-1)*Nuxyz
           iip2 = ijku+((i+2)+jk2-1)*Nuxyz
           iip1 = ijku+((i+1)+jk2-1)*Nuxyz
           ii   = ijku+((i  )+jk2-1)*Nuxyz
           iim1 = ijku+((i-1)+jk2-1)*Nuxyz
           iim2 = ijku+((i-2)+jk2-1)*Nuxyz
           iim3 = ijku+((i-3)+jk2-1)*Nuxyz
           fp(ii)  = CWENO6th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx)
         enddo      
!-------------------------------------------------------------------------------
!  for i = 1
         iip3 = ijku+(4+jk2-1)*Nuxyz
         iip2 = ijku+(3+jk2-1)*Nuxyz
         iip1 = ijku+(2+jk2-1)*Nuxyz
         ii   = ijku+(1+jk2-1)*Nuxyz
!
         fp(ii)  = CWENO6th(fm3(ii1),fm2(ii1),fm1(ii1),f(ii),f(iip1),f(iip2),f(iip3),dx)
         !fpp(ii) = FSD5th(fm3(ii1),fm2(ii1),fm1(ii1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = 2
         iip3 = ijku+(5+jk2-1)*Nuxyz
         iip2 = ijku+(4+jk2-1)*Nuxyz
         iip1 = ijku+(3+jk2-1)*Nuxyz
         ii   = ijku+(2+jk2-1)*Nuxyz
         iim1 = ijku+(1+jk2-1)*Nuxyz     
!
         fp(ii)  = CWENO6th(fm2(ii1),fm1(ii1),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx)
         !fpp(ii) = FSD5th(fm2(ii1),fm1(ii1),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = 3
         iip3 = ijku+(6+jk2-1)*Nuxyz
         iip2 = ijku+(5+jk2-1)*Nuxyz
         iip1 = ijku+(4+jk2-1)*Nuxyz
         ii   = ijku+(3+jk2-1)*Nuxyz
         iim1 = ijku+(2+jk2-1)*Nuxyz
         iim2 = ijku+(1+jk2-1)*Nuxyz     
!
         fp(ii)  = CWENO6th(fm1(ii1),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx)
         !fpp(ii) = FSD5th(fm1(ii1),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx-2
         iip2 = ijku+((Nx  )+jk2-1)*Nuxyz
         iip1 = ijku+((Nx-1)+jk2-1)*Nuxyz
         ii   = ijku+((Nx-2)+jk2-1)*Nuxyz
         iim1 = ijku+((Nx-3)+jk2-1)*Nuxyz
         iim2 = ijku+((Nx-4)+jk2-1)*Nuxyz
         iim3 = ijku+((Nx-5)+jk2-1)*Nuxyz              
!
         fp(ii)  = CWENO6th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),fp1(ii1),dx)
         !fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),fp1(ii1),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx-1
         iip1 = ijku+((Nx  )+jk2-1)*Nuxyz
         ii   = ijku+((Nx-1)+jk2-1)*Nuxyz
         iim1 = ijku+((Nx-2)+jk2-1)*Nuxyz
         iim2 = ijku+((Nx-3)+jk2-1)*Nuxyz
         iim3 = ijku+((Nx-4)+jk2-1)*Nuxyz              
!
         fp(ii)  = CWENO6th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),fp1(ii1),fp2(ii1),dx)
         !fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),fp1(ii1),fp2(ii1),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx
         ii   = ijku+((Nx  )+jk2-1)*Nuxyz
         iim1 = ijku+((Nx-1)+jk2-1)*Nuxyz
         iim2 = ijku+((Nx-2)+jk2-1)*Nuxyz
         iim3 = ijku+((Nx-3)+jk2-1)*Nuxyz              
!
         fp(ii)  = CWENO6th(f(iim3),f(iim2),f(iim1),f(ii),fp1(ii1),fp2(ii1),fp3(ii1),dx)
         !fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),fp1(ii1),fp2(ii1),fp3(ii1),dx2)
 30    continue
!
 35  continue
!$OMP end parallel do
!-------------------------------------------------------------------------------
     deallocate (fm3,fm2,fm1,fp1,fp2,fp3,fsend,frecvp,frecvm)  
   end subroutine CWENOx6th6D_MPI
!===============================================================================
    subroutine CWENOux6th(f,fp,du,Nux,Nuy,Nuz)
      implicit double precision (A-H,O-Z)
      dimension :: f(1), fp(1)
      Nuxy  = Nux*Nuy
      Nuxyz = Nux*Nuy*Nuz    
!------------------------------------------------------------------------------- 
         do iu = 4, Nux-3
           iip3 = iu+3
           iip2 = iu+2
           iip1 = iu+1
           ii   = iu+0
           iim1 = iu-1
           iim2 = iu-2
           iim3 = iu-3
           fp(ii)  = CWENO6th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),du)
         enddo      
!-------------------------------------------------------------------------------
!  for i = 1
         iip3 = 4
         iip2 = 3
         iip1 = 2
         ii   = 1
!
         fp(ii)  = CWENO6th(f(ii),f(ii),f(ii),f(ii),f(iip1),f(iip2),f(iip3),du)
         !fpp(ii) = FSD5th(fm3(ii1),fm2(ii1),fm1(ii1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = 2
         iip3 = 5
         iip2 = 4
         iip1 = 3
         ii   = 2
         iim1 = 1     
!
         fp(ii)  = CWENO6th(f(iim1),f(iim1),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),du)
         !fpp(ii) = FSD5th(fm2(ii1),fm1(ii1),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = 3
         iip3 = 6
         iip2 = 5
         iip1 = 4
         ii   = 3
         iim1 = 2
         iim2 = 1     
!
         fp(ii)  = CWENO6th(f(iim2),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),du)
         !fpp(ii) = FSD5th(fm1(ii1),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx-2
         iip2 = Nux
         iip1 = Nux-1
         ii   = Nux-2
         iim1 = Nux-3
         iim2 = Nux-4
         iim3 = Nux-5              
!
         fp(ii)  = CWENO6th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip2),du)
         !fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),fp1(ii1),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx-1
         iip1 = Nux
         ii   = Nux-1
         iim1 = Nux-2
         iim2 = Nux-3
         iim3 = Nux-4           
!
         fp(ii)  = CWENO6th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip1),f(iip1),du)
         !fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),fp1(ii1),fp2(ii1),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx
         ii   = Nux
         iim1 = Nux-1
         iim2 = Nux-2
         iim3 = Nux-3             
!
         fp(ii)  = CWENO6th(f(iim3),f(iim2),f(iim1),f(ii),f(ii),f(ii),f(ii),du)
         !fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),fp1(ii1),fp2(ii1),fp3(ii1),dx2)
   end subroutine CWENOux6th
!===============================================================================

!===============================================================================
    subroutine WENOux5th(f,fp,du,Nux,Nuy,Nuz,E)
      implicit double precision (A-H,O-Z)
      dimension :: f(1), fp(1)
      Nuxy  = Nux*Nuy
      Nuxyz = Nux*Nuy*Nuz    
!------------------------------------------------------------------------------- 
         do iu = 4, Nux-3
           iip3 = iu+3
           iip2 = iu+2
           iip1 = iu+1
           ii   = iu+0
           iim1 = iu-1
           iim2 = iu-2
           iim3 = iu-3
           fp(ii)  = WENO5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),du,E)
         enddo      
!-------------------------------------------------------------------------------
!  for i = 1
         iip3 = 4
         iip2 = 3
         iip1 = 2
         ii   = 1
!
         fp(ii)  = WENO5th(f(ii),f(ii),f(ii),f(ii),f(iip1),f(iip2),f(iip3),du,E)
         !fpp(ii) = FSD5th(fm3(ii1),fm2(ii1),fm1(ii1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = 2
         iip3 = 5
         iip2 = 4
         iip1 = 3
         ii   = 2
         iim1 = 1     
!
         fp(ii)  = WENO5th(f(iim1),f(iim1),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),du,E)
         !fpp(ii) = FSD5th(fm2(ii1),fm1(ii1),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = 3
         iip3 = 6
         iip2 = 5
         iip1 = 4
         ii   = 3
         iim1 = 2
         iim2 = 1     
!
         fp(ii)  = WENO5th(f(iim2),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),du,E)
         !fpp(ii) = FSD5th(fm1(ii1),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip3),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx-2
         iip2 = Nux
         iip1 = Nux-1
         ii   = Nux-2
         iim1 = Nux-3
         iim2 = Nux-4
         iim3 = Nux-5              
!
         fp(ii)  = WENO5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),f(iip2),du,E)
         !fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip2),fp1(ii1),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx-1
         iip1 = Nux
         ii   = Nux-1
         iim1 = Nux-2
         iim2 = Nux-3
         iim3 = Nux-4           
!
         fp(ii)  = WENO5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),f(iip1),f(iip1),du,E)
         !fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),f(iip1),fp1(ii1),fp2(ii1),dx2)
!-------------------------------------------------------------------------------
!  for i = Nx
         ii   = Nux
         iim1 = Nux-1
         iim2 = Nux-2
         iim3 = Nux-3             
!
         fp(ii)  = WENO5th(f(iim3),f(iim2),f(iim1),f(ii),f(ii),f(ii),f(ii),du,E)
         !fpp(ii) = FSD5th(f(iim3),f(iim2),f(iim1),f(ii),fp1(ii1),fp2(ii1),fp3(ii1),dx2)
   end subroutine WENOux5th
!===============================================================================
!===============================================================================
    function WENOiph(fm2,fm1,f,fp1,fp2,fp3,v)
      implicit double precision (a-h,o-z)
      ah = 0.5d0
!
      temp1 = sumABCD(-fm1,7*f,7*fp1,-fp2)/12.d0
!
      fm32 = DIFAB(fm1,fm2)
      fm12 = DIFAB(f  ,fm1)
      fp12 = DIFAB(fp1,f  )
      fp32 = DIFAB(fp2,fp1)
      fp52 = DIFAB(fp3,fp2)
      if(v .GE. 0) then
        temp2 = phiN(fm32,fm12,fp12,fp32)
      else 
        temp2 = -1.d0*phiN(fp52,fp32,fp12,fm12)
      end if
      WENOiph = DIFAB(temp1,temp2)
    end function WENOiph
!===============================================================================
!===============================================================================
    function WENO5th(fm3,fm2,fm1,f,fp1,fp2,fp3,dx,v)
      implicit double precision (a-h,o-z)
!      WENO5th=(WENOiph(fm2,fm1,f,fp1,fp2)-WENOiph(fm3,fm2,fm1,f,fp1))/dx
      WENO5th=DIFAB(WENOiph(fm2,fm1,f,fp1,fp2,fp3,v),WENOiph(fm3,fm2,fm1,f,fp1,fp2,v))/dx
    end function WENO5th
!===============================================================================
!===============================================================================
    function CWENOiph(fm2,fm1,f,fp1,fp2,fp3)
      implicit double precision (a-h,o-z)
      eps = 1.d-6 
      d_ll = 0.05d0
      d_l  = 0.45d0
      d_r  = 0.45d0
      d_rr = 0.05d0
      tau = (-1.d0*fm2+5.d0*fm1-10.d0*f+10.d0*fp1-5.d0*fp2+fp3)**2.d0
      coef1 = 13.d0/12.d0
      coef2 = 0.25d0
      b_ll = coef1*(fm2-2.d0*fm1+f)**2+coef2*(fm2-4.d0*fm1+3.d0*f)**2
      b_l  = coef1*(fm1-2.d0*f+fp1)**2+coef2*(fm1-fp1)**2
      b_r  = coef1*(f-2.d0*fp1+fp2)**2+coef2*(3.d0*f-4.d0*fp1+fp2)**2
      b_rr = coef1*(fp1-2.d0*fp2+fp3)**2+coef2*(-5.d0*fp1+8.d0*fp2-3.d0*fp3)**2
      a_ll = d_ll*(1+(tau)/(eps+b_ll))
      a_l  = d_l*(1+(tau)/(eps+b_l))
      a_r  = d_r*(1+(tau)/(eps+b_r))
      a_rr = d_rr*(1+(tau)/(eps+b_rr))
      total = sumABCD(a_ll,a_l,a_r,a_rr)
      w_ll = a_ll/total
      w_l  = a_l/total
      w_r  = a_r/total
      w_rr = a_rr/total
      f_ll = (1.d0/6.d0)*(2.d0*fm2-7.d0*fm1+11.d0*f)
      f_l  = (1.d0/6.d0)*(-1.d0*fm1+5.d0*f+2.d0*fp1)
      f_r  = (1.d0/6.d0)*(2.d0*f+5.d0*fp1-1.d0*fp2)
      f_rr = (1.d0/6.d0)*(11.d0*fp1-7.d0*fp2+2.d0*fp3)
      CWENOiph = sumABCD(w_ll*f_ll,w_l*f_l,w_r*f_r,w_rr*f_rr)
    end function CWENOiph
!===============================================================================
!===============================================================================
    function CWENO6th(fm3,fm2,fm1,f,fp1,fp2,fp3,dx)
      implicit double precision (a-h,o-z)
!      WENO5th=(WENOiph(fm2,fm1,f,fp1,fp2)-WENOiph(fm3,fm2,fm1,f,fp1))/dx
      CWENO6th=DIFAB(CWENOiph(fm2,fm1,f,fp1,fp2,fp3),CWENOiph(fm3,fm2,fm1,f,fp1,fp2))/dx
    end function CWENO6th
!===============================================================================
!===============================================================================
    function phiN(a,b,c,d)
      implicit double precision (a-h,o-z)
      e = 1.d-6 
      aIS0 = 13.d0*(DIFAB(a,b))**2+3.d0*(DIFAB(a,3*b))**2
      aIS1 = 13.d0*(DIFAB(b,c))**2+3.d0*(DIFAB(b,-c))**2
      aIS2 = 13.d0*(DIFAB(c,d))**2+3.d0*(DIFAB(3*c,d))**2
!
      alpha0 = 1.d0/(e+aIS0)**2
      alpha1 = 6.d0/(e+aIS1)**2
      alpha2 = 3.d0/(e+aIS2)**2
      alpha_sum = alpha0+alpha1+alpha2
!
      omega0 = alpha0/alpha_sum
      omega2 = alpha2/alpha_sum
!
      temp1 = omega0*sumABC(a,-2*b,c)/3.d0
      temp2 = DIFAB(omega2,0.5d0)*sumABC(b,-2*c,d)/6.d0
      phiN = sumAB(temp1,temp2)
    end function phiN
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
end module WENO3D_MPI
