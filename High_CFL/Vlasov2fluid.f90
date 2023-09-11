module Vlasov2fluid
   use omp_lib
!   use GRID6D
   use CubicSpline
   contains
!===============================================================================
!===============================================================================
!--SUB. dfdv in 3D velocity space
   subroutine gradu(f,fpx,fpy,fpz,Bx,Cx,hx,By,Cy,hy,Bz,Cz,hz,nux,nuy,nuz)
!
!  calculate the first derivative of the distribution function
!    in 3D velocity space
!
!  input array  : f(nux,nuy,nuz)
!  input array  : Bx(nux), Cx(nux), hx(nux)
!  input array  : By(nuy), Cy(nuy), hy(nuy)
!  input array  : Bz(nuz), Cz(nuz), hz(nuz)
!  output array : fpx(nux,nuy,nuz), fpy(nux,nuy,nuz), fpz(nux,nuy,nuz)
!
     implicit double precision (a-h,o-z)
     dimension :: f(1), fpx(1), fpy(1), fpz(1)
     dimension :: Bx(1), Cx(1), hx(1), By(1), Cy(1), hy(1), Bz(1), Cz(1), hz(1)
!
     call gradux(f,fpx,Bx,Cx,hx,nux,nuy,nuz)
!
     call graduy(f,fpy,By,Cy,hy,nux,nuy,nuz)
!
     call graduz(f,fpz,Bz,Cz,hz,nux,nuy,nuz)
!
   end subroutine gradu
!===============================================================================
!===============================================================================
!--SUB. dfdv in 3D velocity space
   subroutine FP_ux_CS(f,fpxx,Bx,Cx,hx,nux,nuy,nuz,vx)
!
!  calculate the first derivative of the distribution function
!    in 3D velocity space
!
!  input array  : f(nux,nuy,nuz)
!  input array  : Bx(nux), Cx(nux), hx(nux)
!  input array  : By(nuy), Cy(nuy), hy(nuy)
!  input array  : Bz(nuz), Cz(nuz), hz(nuz)
!  output array : fpx(nux,nuy,nuz), fpy(nux,nuy,nuz), fpz(nux,nuy,nuz)
!
     implicit double precision (a-h,o-z)
     dimension :: f(1), fpxx(1),vx(1)
     dimension :: Bx(1), Cx(1), hx(1)
     real(8), allocatable, dimension(:) :: FPx, Rx,FP0x
     nuxyz = nux*nuy*nuz
     allocate (FPx(nuxyz),Rx(nuxyz),FP0x(nuxyz))
     wkfx = 0.d0
     Rx = 0.d0
     FPx = 0.d0
     Rx = 0.d0
     FP0x = 0.d0
     !call gradux(f,FPx,Bx,Cx,hx,nux,nuy,nuz)
     !fpxx(1) = 0.d0
     !fpxx(2) = 0.d0
     !fpxx(nux) = 0.d0
     !fpxx(nux-1)=0.d0
     !du = vx(2)-vx(1)
     !do iu = 2,nux-1
     !   fpxx(iu) = (FPx(iu+1)-FPx(iu-1))/(2*du)
     !   fpxx(iu) = fpxx(iu) + f(iu) + vx(iu)*FPx(iu)
     !   if (vx(iu) > 6.d0 ) then
     !      fpxx(iu) = 0.d0
     !   endif
     !enddo
     call getFP0(f,hx,FP0x,nux)
     call FCSYP1(FP0x,Bx,Cx,FPx,nux,Rx)
     call YP2YPP(FPx,FP0x,hx,fpxx,nux)
     do iu =1,nux
        fpxx(iu) = f(iu)+vx(iu)*FPx(iu)+1.d0*fpxx(iu)
     enddo
     deallocate ( FPx, Rx,FP0x) 
!
!
   end subroutine FP_ux_CS
!===============================================================================
!===============================================================================
!--SUB. dfdvx in 3D velocity space
   subroutine gradux(f,fp,Bx,Cx,hx,nux,nuy,nuz)
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
     dimension :: Bx(1), Cx(1), hx(1)
     real(8), allocatable, dimension(:) :: Rx, FPx, FP0x !, wkfx
     nuxy = nux*nuy
!-------------------------------------------------------------------------------
     allocate (Rx(nux), FPx(nux), FP0x(nux)) !, wkfx(nux))     
!
     do k = 1, nuz
     do j = 1, nuy
       jk = nux*(j-1)+nuxy*(k-1)
       ijks =   1+jk
       ijke = nux+jk
       FPx(1)   = 0.d0
       FPx(nux) = 0.d0
!       wkfx(1:nux) = f(ijks:ijke)
!       
       call getFP0(f(ijks:ijke),hx,FP0x,nux)
       call FCSYP1(FP0x,Bx,Cx,FPx,nux,Rx)
!
       fp(ijks:ijke) = FPx(1:nux)
     enddo
     enddo
!
     deallocate (Rx, FPx, FP0x) !, wkfx)
!-------------------------------------------------------------------------------  
   end subroutine gradux
!===============================================================================
!--SUB. dfdvy in 3D velocity space
   subroutine graduy(f,fp,By,Cy,hy,nux,nuy,nuz)
!
!  calculate the first derivative of the distribution function
!    along y-axis in 3D velocity space
!
!  input array   : f(nux,nuy,nuz)
!  input array   : By(nuy), Cy(nuy), hy(nuy)
!  working array : Ry(nuy), FPy(nuy), FP0y(nuy), wkfy(nuy)
!  output array  : fp(nux,nuy,nuz)
!
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1)
     dimension :: By(1), Cy(1), hy(1)
     real(8), allocatable, dimension(:) :: Ry, FPy, FP0y, wkfy
     nuxy = nux*nuy
     nuxyz = nux*nuy*nuz
!-------------------------------------------------------------------------------
     if (nuy == 1) then
       fp(1:nuxyz) = 0.d0
       return
     endif
!------------------------------------------------------------------------------- 
     allocate (Ry(nuy), FPy(nuy), FP0y(nuy), wkfy(nuy))     
!
     do k = 1, nuz
     do i = 1, nux
       ik = i+(k-1)*nuxy
       do j = 1, nuy
         ijk = (j-1)*nux+ik
         wkfy(j) = f(ijk)
       enddo
!       
       FPy(1)   = 0.d0
       FPy(nuy) = 0.d0
       call getFP0(wkfy,hy,FP0y,nuy)
       call FCSYP1(FP0y,By,Cy,FPy,nuy,Ry)
!
       do j = 1, nuy
         ijk = (j-1)*nux+ik
         fp(ijk) = FPy(j)
       enddo
     enddo
     enddo
!
     deallocate (Ry, FPy, FP0y, wkfy)
!-------------------------------------------------------------------------------  
   end subroutine graduy
!===============================================================================
!--SUB. dfdvy in 3D velocity space
   subroutine graduz(f,fp,Bz,Cz,hz,nux,nuy,nuz)
!
!  calculate the first derivative of the distribution function
!    along z-axis in 3D velocity space
!
!  input array   : f(nux,nuy,nuz)
!  input array   : Bz(nuz), Cz(nuz), hz(nuz)
!  working array : Rz(nuz), FPz(nuz), FP0z(nuz), wkfz(nuz)
!  output array  : fp(nux,nuy,nuz)
!
     implicit double precision (A-H,O-Z)
     dimension :: f(1), fp(1)
     dimension :: Bz(1), Cz(1), hz(1)
     real(8), allocatable, dimension(:) :: Rz, FPz, FP0z, wkfz
     nuxy = nux*nuy
     nuxyz = nux*nuy*nuz
!-------------------------------------------------------------------------------
     if (nuz == 1) then
       fp(1:nuxyz) = 0.d0
       return
     endif
!-------------------------------------------------------------------------------
     allocate (Rz(nuz), FPz(nuz), FP0z(nuz), wkfz(nuz))     
!
     do j = 1, nuy
     do i = 1, nux
       ij = i+(j-1)*nux
       do k = 1, nuz
         ijk = (k-1)*nuxy+ij
         wkfz(k) = f(ijk)
       enddo
!       
       FPz(1)   = 0.d0
       FPz(nuz) = 0.d0
       call getFP0(wkfz,hz,FP0z,nuz)
       call FCSYP1(FP0z,Bz,Cz,FPz,nuz,Rz)
!
       do k = 1, nuz
         ijk = (k-1)*nuxy+ij
         fp(ijk) = FPz(k)
       enddo
     enddo
     enddo
!
     deallocate (Rz, FPz, FP0z, wkfz)
!-------------------------------------------------------------------------------  
   end subroutine graduz
!===============================================================================
!===============================================================================
   subroutine get_distribution_1D(nx,ny,nz,nux,nuy,nuz,f,Bx,Cx,hx,By,Cy,hy,Bz,Cz,hz, &
                                  fux,fuy,fuz)
!
!  calculate the 1D distribution fux(nux), fuy(nuy), fuz(nuz) 
!    from 3D distribution function fu(nux,nuy,nuz)
!
!  input array  : f(nux,nuy,nuz,nx,ny,nz)
!  input array  : Bx(nux), Cx(nux), hx(nux)
!  input array  : By(nuy), Cy(nuy), hy(nuy)
!  input array  : Bz(nuz), Cz(nuz), hz(nuz)
!  output array : fux(nux,nx,ny,nz), fuy(nuy,nx,ny,nz), fuz(nuz,nx,ny,nz)
!
     implicit double precision (a-h,o-z)
     real(8) :: f(1), fux(1), fuy(1), fuz(1)
     real(8) :: Bx(1), Cx(1), hx(1), By(1), Cy(1), hy(1), Bz(1), Cz(1), hz(1)
     real(8), allocatable, dimension(:) :: wkfux, wkfuy, wkfuz
!
     nxy   = nx*ny
     nuxyz = nux*nuy*nuz
!$OMP parallel private(ijk,iis,iie,iixs,iixe,iiys,iiye,iizs,iize,wkfux,wkfuy,wkfuz)
     allocate (wkfux(nux), wkfuy(nuy), wkfuz(nuz))
     wkfux = 0.d0
     wkfuy = 0.d0
     wkfuz = 0.d0
!
!$OMP do collapse(3)
     do 131 k = 1, nz
     do 131 j = 1, ny
     do 131 i = 1, nx
       ijk = i+(j-1)*nx+(k-1)*nxy
       iis = 1+nuxyz*(ijk-1)
       iie =   nuxyz* ijk
!
       call integralyz(nux,nuy,nuz,f(iis:iie),wkfux,By,Cy,hy,Bz,Cz,hz)
       call integralxz(nux,nuy,nuz,f(iis:iie),wkfuy,Bx,Cx,hx,Bz,Cz,hz)
       call integralxy(nux,nuy,nuz,f(iis:iie),wkfuz,Bx,Cx,hx,By,Cy,hy)
!
       iixs = 1+nux*(ijk-1)
       iixe =   nux*(ijk) 
       iiys = 1+nuy*(ijk-1)
       iiye =   nuy*(ijk)
       iizs = 1+nuz*(ijk-1)
       iize =   nuz*(ijk)
!
       fux(iixs:iixe) = wkfux(1:nux)
       fuy(iiys:iiye) = wkfuy(1:nuy)
       fuz(iizs:iize) = wkfuz(1:nuz)
 131 continue
!$OMP end do
     deallocate (wkfux, wkfuy, wkfuz)
!$OMP end parallel
!     
   end subroutine get_distribution_1D
!===============================================================================
   subroutine currentX(nx,ny,nz,nux,nuy,nuz,f,vx,Bx,Cx,hx,By,Cy,hy,Bz,Cz,hz,Jx)
!
!  input array  : f(nux,nuy,nuz,nx,ny,nz)
!  input array  : vx(nux,nuy,nuz)
!  input array  : Bx(nux), Cx(nux), hx(nux)
!  input array  : By(nuy), Cy(nuy), hy(nuy)
!  input array  : Bz(nuz), Cz(nuz), hz(nuz)
!  output array : Jx(nx,ny,nz)
!
     implicit double precision (a-h,o-z)
     real(8) :: f(1), vx(1), Jx(1)
     real(8) :: Bx(1), Cx(1), hx(1), By(1), Cy(1), hy(1), Bz(1), Cz(1), hz(1)
     real(8), allocatable, dimension(:) :: wkfx
!
     nxy   = nx*ny
     nuxy  = nux*nuy
     nuxyz = nux*nuy*nuz
!
!$OMP parallel private(ii,ijk,ijku,wkfx)
     allocate (wkfx(nuxyz))
     wkfx = 0.d0
!
!$OMP do collapse(3)
     do 131 k = 1, nz
     do 131 j = 1, ny
     do 131 i = 1, nx
!
       ijk = i+(j-1)*nx+(k-1)*nxy
!      
       do 132 ku = 1, nuz
       do 132 ju = 1, nuy
       do 132 iu = 1, nux
         ijku = iu+(ju-1)*nux+(ku-1)*nuxy
         ii = ijku+nuxyz*(ijk-1)
!
         wkfx(ijku) = f(ii)*vx(ijku)
 132   continue
!
       call integral3D(nux,nuy,nuz,wkfx,Jx(ijk),Bx,Cx,hx,By,Cy,hy,Bz,Cz,hz)
 131 continue
!$OMP end do 
     deallocate (wkfx)
!$OMP end parallel
!
   end subroutine currentX
!===============================================================================
   subroutine current(nx,ny,nz,nux,nuy,nuz,f,vx,vy,vz, &
                      Bx,Cx,hx,By,Cy,hy,Bz,Cz,hz,Jx,Jy,Jz)
!
!  input array  : f(nux,nuy,nuz,nx,ny,nz)
!  input array  : vx(nux,nuy,nuz), vy(nux,nuy,nuz), vz(nux,nuy,nuz)
!  input array  : Bx(nux), Cx(nux), hx(nux)
!  input array  : By(nuy), Cy(nuy), hy(nuy)
!  input array  : Bz(nuz), Cz(nuz), hz(nuz)
!  output array : Jx(nx,ny,nz), Jy(nx,ny,nz), Jz(nx,ny,nz)
!
     implicit double precision (a-h,o-z)
     real(8) :: f(1), vx(1), vy(1), vz(1), Jx(1), Jy(1), Jz(1)
     real(8) :: Bx(1), Cx(1), hx(1), By(1), Cy(1), hy(1), Bz(1), Cz(1), hz(1)
     real(8), allocatable, dimension(:) :: wkfx, wkfy, wkfz
!
     nxy   = nx*ny
     nuxy  = nux*nuy
     nuxyz = nux*nuy*nuz
!
!$OMP parallel private(ii,ijk,ijku,wkfx,wkfy,wkfz)
     allocate (wkfx(nuxyz), wkfy(nuxyz), wkfz(nuxyz))
     wkfx = 0.d0
     wkfy = 0.d0
     wkfz = 0.d0
!
!$OMP do collapse(3)
     do 131 k = 1, nz
     do 131 j = 1, ny
     do 131 i = 1, nx
!
       ijk = i+(j-1)*nx+(k-1)*nxy
!      
       do 132 ku = 1, nuz
       do 132 ju = 1, nuy
       do 132 iu = 1, nux
         ijku = iu+(ju-1)*nux+(ku-1)*nuxy
         ii = ijku+nuxyz*(ijk-1)
!
         wkfx(ijku) = f(ii)*vx(ijku)
         wkfy(ijku) = f(ii)*vy(ijku)
         wkfz(ijku) = f(ii)*vz(ijku)      
 132   continue
!
       call integral3D(nux,nuy,nuz,wkfx,Jx(ijk),Bx,Cx,hx,By,Cy,hy,Bz,Cz,hz)
       call integral3D(nux,nuy,nuz,wkfy,Jy(ijk),Bx,Cx,hx,By,Cy,hy,Bz,Cz,hz)
       call integral3D(nux,nuy,nuz,wkfz,Jz(ijk),Bx,Cx,hx,By,Cy,hy,Bz,Cz,hz)
 131 continue
!$OMP end do 
     deallocate (wkfx, wkfy, wkfz)
!$OMP end parallel
!
   end subroutine current
!===============================================================================
   subroutine density(nx,ny,nz,nux,nuy,nuz,f,Bx,Cx,hx,By,Cy,hy,Bz,Cz,hz,Rho)
!
!  input array  : f(nux,nuy,nuz,nx,ny,nz)
!  input array  : Bx(nux), Cx(nux), hx(nux)
!  input array  : By(nuy), Cy(nuy), hy(nuy)
!  input array  : Bz(nuz), Cz(nuz), hz(nuz)
!  output array : Rho(nx,ny,nz)
!
     implicit double precision (a-h,o-z)
     real(8) :: f(1), Rho(1)
     real(8) :: Bx(1), Cx(1), hx(1), By(1), Cy(1), hy(1), Bz(1), Cz(1), hz(1)
     nxy   = nx*ny
!     nuxy  = nux*nuy
     nuxyz = nux*nuy*nuz
!
!$OMP parallel do private(ijk,iis,iie) collapse(3)
     do 131 k = 1, nz
     do 131 j = 1, ny
     do 131 i = 1, nx
       ijk = i+(j-1)*nx+(k-1)*nxy
       iis = 1+nuxyz*(ijk-1)
       iie =   nuxyz* ijk
       call integral3D(nux,nuy,nuz,f(iis:iie),Rho(ijk),Bx,Cx,hx,By,Cy,hy,Bz,Cz,hz)
 131 continue
!$OMP end parallel do
   end subroutine density
!===============================================================================
   subroutine pressure(nx,ny,nz,nux,nuy,nuz,f,am,Bx,Cx,hx,By,Cy,hy,Bz,Cz,hz, &
                       ux,uy,uz,avgux,avguy,avguz,vx,vy,vz,avgvx,avgvy,avgvz, &
                       Px,Py,Pz)
!
!  input array  : f(nux,nuy,nuz,nx,ny,nz)
!  input array  : ux(nux), uy(nuy), uz(nuz)
!  input array  : vx(nux,nuy,nuz), vy(nux,nuy,nuz), vz(nux,nuy,nuz)
!  input array  : avgux(nx,ny,nz), avguy(nx,ny,nz), avguz(nx,ny,nz)
!  input array  : avgvx(nx,ny,nz), avgvy(nx,ny,nz), avgvz(nx,ny,nz)
!  input array  : Bx(nux), Cx(nux), hx(nux)
!  input array  : By(nuy), Cy(nuy), hy(nuy)
!  input array  : Bz(nuz), Cz(nuz), hz(nuz)
!  input data   : am(mass ratio)
!  output array : Px(nx,ny,nz), Py(nx,ny,nz), Pz(nx,ny,nz)
!
     implicit double precision (a-h,o-z)
     real(8), allocatable, dimension(:) :: wkfx, wkfy, wkfz
     real(8) :: f(1), Px(1), Py(1), Pz(1)
     real(8) :: ux(1), uy(1), uz(1), avgux(1), avguy(1), avguz(1)
     real(8) :: vx(1), vy(1), vz(1), avgvx(1), avgvy(1), avgvz(1)
     real(8) :: Bx(1), Cx(1), hx(1), By(1), Cy(1), hy(1), Bz(1), Cz(1), hz(1)
     nxy   = nx*ny
     nxyz  = nx*ny*nz
     nuxy  = nux*nuy
     nuxyz = nux*nuy*nuz
!
!$OMP parallel private(ii,ijk,ijku,wkfx,wkfy,wkfz)
     allocate (wkfx(nuxyz), wkfy(nuxyz), wkfz(nuxyz)) 
     wkfx = 0.d0
     wkfy = 0.d0
     wkfz = 0.d0  
!
!$OMP do collapse(3)
     do 15 k = 1, nz
     do 15 j = 1, ny
     do 15 i = 1, nx
!
       ijk = i+(j-1)*nx+(k-1)*nxy
!     
       do 16 ku = 1, nuz
       do 16 ju = 1, nuy
       do 16 iu = 1, nux
         ijku  = iu+(ju-1)*nux+(ku-1)*nuxy
         ii    = ijku+nuxyz*(ijk-1)
         wkfx(ijku) = (ux(iu)-avgux(ijk))*(vx(ijku)-avgvx(ijk))*f(ii)
         wkfy(ijku) = (uy(ju)-avguy(ijk))*(vy(ijku)-avgvy(ijk))*f(ii)
         wkfz(ijku) = (uz(ku)-avguz(ijk))*(vz(ijku)-avgvz(ijk))*f(ii)     
 16    continue
       call integral3D(nux,nuy,nuz,wkfx,Px(ijk),Bx,Cx,hx,By,Cy,hy,Bz,Cz,hz)
       call integral3D(nux,nuy,nuz,wkfy,Py(ijk),Bx,Cx,hx,By,Cy,hy,Bz,Cz,hz)
       call integral3D(nux,nuy,nuz,wkfz,Pz(ijk),Bx,Cx,hx,By,Cy,hy,Bz,Cz,hz)
 15  continue
!$OMP end do
     deallocate (wkfx, wkfy, wkfz)
!$OMP end parallel
!
     Px(1:nxyz) = am*Px(1:nxyz)
     Py(1:nxyz) = am*Py(1:nxyz)
     Pz(1:nxyz) = am*Pz(1:nxyz)
   end subroutine pressure
!===============================================================================
   subroutine pressureX(nx,ny,nz,nux,nuy,nuz,f,am,Bx,Cx,hx,By,Cy,hy,Bz,Cz,hz, &
                        ux,avgux,vx,avgvx,Px)
!
!  input array  : f(nux,nuy,nuz,nx,ny,nz)
!  input array  : ux(nux), vx(nux,nuy,nuz)
!  input array  : avgux(nx,ny,nz), avgvx(nx,ny,nz)
!  input array  : Bx(nux), Cx(nux), hx(nux)
!  input array  : By(nuy), Cy(nuy), hy(nuy)
!  input array  : Bz(nuz), Cz(nuz), hz(nuz)
!  input data   : am(mass ratio)
!  output array : Px(nx,ny,nz)
!
     implicit double precision (a-h,o-z)
     real(8), allocatable, dimension(:) :: wkfx
     real(8) :: f(1), Px(1)
     real(8) :: ux(1), avgux(1), vx(1), avgvx(1)
     real(8) :: Bx(1), Cx(1), hx(1), By(1), Cy(1), hy(1), Bz(1), Cz(1), hz(1)
     nxy   = nx*ny
     nxyz  = nx*ny*nz
     nuxy  = nux*nuy
     nuxyz = nux*nuy*nuz
!
!$OMP parallel private(ii,ijk,ijku,wkfx)
     allocate (wkfx(nuxyz)) 
     wkfx = 0.d0  
!
!$OMP do collapse(3)
     do 15 k = 1, nz
     do 15 j = 1, ny
     do 15 i = 1, nx
!
       ijk = i+(j-1)*nx+(k-1)*nxy
!     
       do 16 ku = 1, nuz
       do 16 ju = 1, nuy
       do 16 iu = 1, nux
         ijku  = iu+(ju-1)*nux+(ku-1)*nuxy
         ii    = ijku+nuxyz*(ijk-1)
         wkfx(ijku) = (ux(iu)-avgux(ijk))*(vx(ijku)-avgvx(ijk))*f(ii)    
 16    continue
       call integral3D(nux,nuy,nuz,wkfx,Px(ijk),Bx,Cx,hx,By,Cy,hy,Bz,Cz,hz)
 15  continue
!$OMP end do
     deallocate (wkfx)
!$OMP end parallel
!
     Px(1:nxyz) = am*Px(1:nxyz)
   end subroutine pressureX
!===============================================================================
!===============================================================================
   subroutine integralyz(nx,ny,nz,f,fx,By,Cy,hy,Bz,Cz,hz)
!
!  integral along y and z direction using cubic spline method
!
!  input data    : nx, ny, nz
!  input array   : f(nx,ny,nz)
!  input array   : By(ny), Cy(ny), hy(ny), Bz(nz), Cz(nz), hz(nz)
!  working array : Ry(ny), FPy(ny), FP0y(ny), alphay(ny), Betay(ny), Areay(ny)
!  working array : Rz(nz), FPz(nz), FP0z(nz), alphaz(nz), Betaz(nz), Areaz(nz)
!  output data   : fx(nx)
!
     implicit double precision (a-h,o-z)
     dimension :: f(1), fx(1)
     dimension :: By(1), Cy(1), hy(1)
     dimension :: Bz(1), Cz(1), hz(1) 
     real*8, allocatable, dimension(:) :: Ry, FPy, FP0y, alphay, Betay, Areay
     real*8, allocatable, dimension(:) :: Rz, FPz, FP0z, alphaz, Betaz, Areaz
     real*8, allocatable, dimension(:) :: wkfy, wkfz, wkfxy
     nxy = nx*ny
!
     if (ny == 1 .and. nz == 1) then
       fx(1:nx) = f(1:nx)
       return
     endif
!
     allocate (wkfxy(nxy))
!
     if (nz == 1) then
       wkfxy(1:nxy) = f(1:nxy)
       goto 123
     endif   
!-------------------------------------------------------------------------------
!  integral along z-axis
!
     allocate (wkfz(nz))
     allocate (Rz(nz), FPz(nz), FP0z(nz), alphaz(nz), Betaz(nz), Areaz(nz))
!
     do j = 1, ny
     do i = 1, nx
!
       ij = i+nx*(j-1)
       do k = 1, nz
         ijk = ij+nxy*(k-1)
         wkfz(k) = f(ijk)     
       enddo
!
       FPz(1)  = 0.d0
       FPz(nz) = 0.d0
       call getFP0(wkfz,hz,FP0z,nz)
       call FCSYP1(FP0z,Bz,Cz,FPz,nz,Rz)
       call getAB(FPz,FP0z,hz,nz,Alphaz,Betaz)
       call getArea(wkfz,hz,Alphaz,Betaz,Areaz,AreaS,nz)
       wkfxy(ij) = AreaS
     enddo
     enddo
!
     deallocate (wkfz)
     deallocate (Rz, FPz, FP0z, alphaz, Betaz, Areaz)
!-------------------------------------------------------------------------------
!  integral along y-axis
!     
 123 continue
     allocate (wkfy(ny))
     allocate (Ry(ny), FPy(ny), FP0y(ny), alphay(ny), Betay(ny), Areay(ny))
     do i = 1, nx
       do j = 1, ny
         ij = i+nx*(j-1)
         wkfy(j) = wkfxy(ij)
       enddo
!
       FPy(1)  = 0.d0
       FPy(ny) = 0.d0
       call getFP0(wkfy,hy,FP0y,ny)
       call FCSYP1(FP0y,By,Cy,FPy,ny,Ry)
       call getAB(FPy,FP0y,hy,ny,Alphay,Betay)
       call getArea(wkfy,hy,Alphay,Betay,Areay,AreaS,ny)
       fx(i) = AreaS
     enddo
!
     deallocate (Ry, FPy, FP0y, alphay, Betay, Areay)
!-------------------------------------------------------------------------------    
     deallocate (wkfy,wkfxy)
   end subroutine integralyz
!===============================================================================
   subroutine integralxz(nx,ny,nz,f,fy,Bx,Cx,hx,Bz,Cz,hz)
!
!  integral along x and z direction using cubic spline method
!
!  input data    : nx, ny, nz
!  input array   : f(nx,ny,nz)
!  input array   : Bx(nx), Cx(nx), hx(nx), Bz(nz), Cz(nz), hz(nz)
!  working array : Rx(nx), FPx(nx), FP0x(nx), alphax(nx), Betax(nx), Areax(nx)
!  working array : Rz(nz), FPz(nz), FP0z(nz), alphaz(nz), Betaz(nz), Areaz(nz)
!  output data : fy(ny)
!
     implicit double precision (a-h,o-z)
     dimension :: f(1), fy(1)
     dimension :: Bx(1), Cx(1), hx(1)
     dimension :: Bz(1), Cz(1), hz(1)
     real*8, allocatable, dimension(:) :: Rx, FPx, FP0x, alphax, Betax, Areax
     real*8, allocatable, dimension(:) :: Rz, FPz, FP0z, alphaz, Betaz, Areaz
     real*8, allocatable, dimension(:) :: wkfyz, wkfx, wkfz
!
     if (ny == 1 .and. nz == 1) return
!
     nyz = ny*nz
     nxy = nx*ny
     allocate (wkfyz(nyz), wkfx(nx))    
!-------------------------------------------------------------------------------
!  integral along x-axis
!     
     allocate (Rx(nx), FPx(nx), FP0x(nx), alphax(nx), Betax(nx), Areax(nx))
!
     do k = 1, nz
     do j = 1, ny
       jk1 = nx*(j-1)+nxy*(k-1)
       jk2 = j+ny*(k-1)
       ijks =  1+jk1
       ijke = nx+jk1
       wkfx(1:nx) = f(ijks:ijke)
!       
       FPx(1)  = 0.d0
       FPx(nx) = 0.d0
       call getFP0(wkfx,hx,FP0x,nx)
       call FCSYP1(FP0x,Bx,Cx,FPx,nx,Rx)
       call getAB(FPx,FP0x,hx,nx,Alphax,Betax)
       call getArea(wkfx,hx,Alphax,Betax,Areax,AreaS,nx)
!
       wkfyz(jk2) = AreaS 
     enddo
     enddo
!
     deallocate (Rx, FPx, FP0x, alphax, Betax, Areax)
!-------------------------------------------------------------------------------
!  integral along z-axis
!
     if (nz == 1) then
       fy(1:ny) = wkfyz(1:nyz)
       return
     endif
!     
     allocate (wkfz(nz))
     allocate (Rz(nz), FPz(nz), FP0z(nz), alphaz(nz), Betaz(nz), Areaz(nz))
!
     do j = 1, ny
       do k = 1, nz
         jk = j+ny*(k-1)
         wkfz(k) = wkfyz(jk)
       enddo
!
       FPz(1)  = 0.d0
       FPz(nz) = 0.d0
       call getFP0(wkfz,hz,FP0z,nz)
       call FCSYP1(FP0z,Bz,Cz,FPz,nz,Rz)
       call getAB(FPz,FP0z,hz,nz,Alphaz,Betaz)
       call getArea(wkfz,hz,Alphaz,Betaz,Areaz,AreaS,nz)
!
       fy(j) = AreaS
     enddo
!
     deallocate (Rz, FPz, FP0z, alphaz, Betaz, Areaz)
!-------------------------------------------------------------------------------    
     deallocate (wkfyz,wkfz,wkfx)
   end subroutine integralxz
!===============================================================================
!===============================================================================
   subroutine integralxy(nx,ny,nz,f,fz,Bx,Cx,hx,By,Cy,hy)
!
!  integral along x and y direction using cubic spline method
!
!  input data    : nx, ny, nz
!  input array   : f(nx,ny,nz)
!  input array   : Bx(nx), Cx(nx), hx(nx), By(ny), Cy(ny), hy(ny)
!  working array : Rx(nx), FPx(nx), FP0x(nx), alphax(nx), Betax(nx), Areax(nx)
!  working array : Ry(ny), FPy(ny), FP0y(ny), alphay(ny), Betay(ny), Areay(ny)
!  output data   : fz(nz)
!
     implicit double precision (a-h,o-z)
     dimension :: f(1), fz(1)
     dimension :: Bx(1), Cx(1), hx(1)
     dimension :: By(1), Cy(1), hy(1)
     real*8, allocatable, dimension(:) :: Rx, FPx, FP0x, alphax, Betax, Areax
     real*8, allocatable, dimension(:) :: Ry, FPy, FP0y, alphay, Betay, Areay
     real*8, allocatable, dimension(:) ::  wkfyz, wkfx, wkfy
!
     if (ny == 1 .or. nz == 1) return
!
     nyz = ny*nz
     nxy = nx*ny
     allocate (wkfyz(nyz), wkfx(nx), wkfy(ny))    
!-------------------------------------------------------------------------------
!  integral along x-axis
!     
     allocate (Rx(nx), FPx(nx), FP0x(nx), alphax(nx), Betax(nx), Areax(nx))
!
     do k = 1, nz
     do j = 1, ny
       jk1 = nx*(j-1)+nxy*(k-1)
       jk2 = j+ny*(k-1)
       ijks =  1+jk1
       ijke = nx+jk1
       wkfx(1:nx) = f(ijks:ijke)
!       
       FPx(1)  = 0.d0
       FPx(nx) = 0.d0
       call getFP0(wkfx,hx,FP0x,nx)
       call FCSYP1(FP0x,Bx,Cx,FPx,nx,Rx)
       call getAB(FPx,FP0x,hx,nx,Alphax,Betax)
       call getArea(wkfx,hx,Alphax,Betax,Areax,AreaS,nx)
!
       wkfyz(jk2) = AreaS
     enddo
     enddo
!
     deallocate (Rx, FPx, FP0x, alphax, Betax, Areax)
!-------------------------------------------------------------------------------
!  integral along y-axis
!
     allocate (Ry(ny), FPy(ny), FP0y(ny), alphay(ny), Betay(ny), Areay(ny))
!
     do k = 1, nz
       jks =  1+ny*(k-1)
       jke =    ny* k
       wkfy(1:ny) = wkfyz(jks:jke)
!
       FPy(1)  = 0.d0
       FPy(ny) = 0.d0
       call getFP0(wkfy,hy,FP0y,ny)
       call FCSYP1(FP0y,By,Cy,FPy,ny,Ry)
       call getAB(FPy,FP0y,hy,ny,Alphay,Betay)
       call getArea(wkfy,hy,Alphay,Betay,Areay,AreaS,ny)
!
       fz(k) = AreaS
     enddo
!
     deallocate (Ry, FPy, FP0y, alphay, Betay, Areay)
!-------------------------------------------------------------------------------    
     deallocate (wkfyz,wkfx,wkfy)
   end subroutine integralxy
!===============================================================================
!===============================================================================
   subroutine integral3D(nx,ny,nz,f,csint,Bx,Cx,hx,By,Cy,hy,Bz,Cz,hz)
!
!  calculate intgration using cubic spline method
!
!  input data    : nx, ny, nz, dx, dy, dz
!  input array   : f(nx,ny,nz)
!  input array   : Bx(nx), Cx(nx), hx(nx)
!  input array   : By(ny), Cy(ny), hy(ny)
!  input array   : Bz(nz), Cz(nz), hz(nz)
!  working array : Rx(nx), FPx(nx), FP0x(nx), alphax(nx), Betax(nx), Areax(nx)
!  working array : Ry(ny), FPy(ny), FP0y(ny), alphay(ny), Betay(ny), Areay(ny)
!  working array : Rz(nz), FPz(nz), FP0z(nz), alphaz(nz), Betaz(nz), Areaz(nz)
!  output data : csint
!
     implicit double precision (a-h,o-z)
     dimension :: f(1)
     dimension :: Bx(1), Cx(1), hx(1)
     dimension :: By(1), Cy(1), hy(1)
     dimension :: Bz(1), Cz(1), hz(1)
     real*8, allocatable, dimension(:) :: Rx, FPx, FP0x, alphax, Betax, Areax
     real*8, allocatable, dimension(:) :: Ry, FPy, FP0y, alphay, Betay, Areay
     real*8, allocatable, dimension(:) :: Rz, FPz, FP0z, alphaz, Betaz, Areaz
     real*8, allocatable, dimension(:) :: wkfyz, wkfx, wkfy, wkfz
!
     nyz = ny*nz
     nxy = nx*ny
     allocate (wkfyz(nyz), wkfx(nx))    
!-------------------------------------------------------------------------------
!  integral along x-axis
!
     allocate (Rx(nx), FPx(nx), FP0x(nx), alphax(nx), Betax(nx), Areax(nx))     
!
     do k = 1, nz
     do j = 1, ny
       jk1 = nx*(j-1)+nxy*(k-1)
       jk2 = j+ny*(k-1)
       ijks =  1+jk1
       ijke = nx+jk1
       wkfx(1:nx) = f(ijks:ijke)
!       
       FPx(1)  = 0.d0
       FPx(nx) = 0.d0
       call getFP0(wkfx,hx,FP0x,nx)
       call FCSYP1(FP0x,Bx,Cx,FPx,nx,Rx)
       call getAB(FPx,FP0x,hx,nx,Alphax,Betax)
       call getArea(wkfx,hx,Alphax,Betax,Areax,AreaS,nx)
!
       wkfyz(jk2) = AreaS 
     enddo
     enddo
!
     deallocate (wkfx)
     deallocate (Rx, FPx, FP0x, alphax, Betax, Areax)
!
!  if ny and nz are equal to 1, then return.
!
     if (ny == 1 .and. nz == 1) then
       csint = wkfyz(1)
       deallocate (wkfyz)
       return
     endif
     allocate (wkfy(ny), wkfz(nz))
!-------------------------------------------------------------------------------
!  integral along y-axis
!
     allocate (Ry(ny), FPy(ny), FP0y(ny), alphay(ny), Betay(ny), Areay(ny))
!
     do k = 1, nz
       jks =  1+ny*(k-1)
       jke =    ny* k
       wkfy(1:ny) = wkfyz(jks:jke)
!
       FPy(1)  = 0.d0
       FPy(ny) = 0.d0
       call getFP0(wkfy,hy,FP0y,ny)
       call FCSYP1(FP0y,By,Cy,FPy,ny,Ry)
       call getAB(FPy,FP0y,hy,ny,Alphay,Betay)
       call getArea(wkfy,hy,Alphay,Betay,Areay,AreaS,ny)
!
       wkfz(k) = AreaS
     enddo
!
     deallocate (wkfyz,wkfy)
     deallocate (Ry, FPy, FP0y, alphay, Betay, Areay)
!
!  if nz is equal to 1, then return.
!
     if (nz == 1) then
       csint = wkfz(1)
       deallocate (wkfz)
       return
     endif
!-------------------------------------------------------------------------------
!  integral along z-axis
!
     allocate (Rz(nz), FPz(nz), FP0z(nz), alphaz(nz), Betaz(nz), Areaz(nz))
!
     FPz(1)  = 0.d0
     FPz(nz) = 0.d0
     call getFP0(wkfz,hz,FP0z,nz)
     call FCSYP1(FP0z,Bz,Cz,FPz,nz,Rz)
     call getAB(FPz,FP0z,hz,nz,Alphaz,Betaz)
     call getArea(wkfz,hz,Alphaz,Betaz,Areaz,csint,nz)
!
     deallocate (Rz, FPz, FP0z, alphaz, Betaz, Areaz)
!-------------------------------------------------------------------------------    
     deallocate (wkfz)
   end subroutine integral3D 
!===============================================================================
!===============================================================================
end module Vlasov2fluid
