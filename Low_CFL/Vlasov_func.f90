!===============================================================================
!===============================================================================
   subroutine Vlasov_EM(time,f_total,func_total,Nt)
!     use coef
!     use comm_global
!     use MPI3D
!     use Vlasov
     use GRID6D
     use vector3D
     use Vlasov2fluid
     use omp_lib
     implicit double precision (a-h,o-z)
     dimension :: f_total(Nt), func_total(Nt)
     real(8) :: tempA, tempB, tempC
     func_total(1:NT) = 0.d0
!===============================================================================
!  Vlasov equation for electrons
!-------------------------------------------------------------------------------
!  Obtain dfe/dx, dfe/dy, dfe/dx, d2fe/dx2, d2fe/dy2, d2fe/dz2
!-------------------------------------------------------------------------------
     allocate (wkfepx(nuexyz*ncxyz_mpi), wkfeppx(nuexyz*ncxyz_mpi))
     allocate (wkfepy(nuexyz*ncxyz_mpi), wkfeppy(nuexyz*ncxyz_mpi))
     allocate (wkfepz(nuexyz*ncxyz_mpi), wkfeppz(nuexyz*ncxyz_mpi))
!
     wkfepx  = 0.d0
     wkfeppx = 0.d0
     wkfepy  = 0.d0
     wkfeppy = 0.d0
     wkfepz  = 0.d0
     wkfeppz = 0.d0
!
     call grad6d_5th(ncx_mpi,ncy_mpi,ncz_mpi,nuex,nuey,nuez, &
                     f_total(ncfe),wkfepx,wkfepy,wkfepz, &
                     wkfeppx,wkfeppy,wkfeppz,dx,dy,dz,ixp,iyp,izp)
!
!$OMP parallel do private(ii,ijk,ijku,tempA,tempB,tempC,tempD) 
     do 112 k = 1, ncz_mpi
     do 112 j = 1, ncy_mpi
     do 112 i = 1, ncx_mpi         
     do 112 ku = 1, nuez
     do 112 ju = 1, nuey
     do 112 iu = 1, nuex
       if (iu == 1 .or. iu == nuex) goto 112
       if (nuey /= 1 .and. (ju == 1 .or. ju == nuey)) goto 112
       if (nuez /= 1 .and. (ku == 1 .or. ku == nuez)) goto 112
!
       ijku = iu+(ju-1)*nuex+(ku-1)*nuexy
       ijk  = i+(j-1)*ncx_mpi+(k-1)*ncxy_mpi
       ii   = ijku+(ijk-1)*nuexyz
!
       tempA = -vex(ijku)*wkfepx(ii)
       tempB = -vey(ijku)*wkfepy(ii)
       tempC = -vez(ijku)*wkfepz(ii)
       tempD = +etafe(ijk)*sumABC(wkfeppx(ii),wkfeppy(ii),wkfeppz(ii))
       func_total(ii) = sumABCD(tempA,tempB,tempC,tempD)
 112 continue
!$OMP end parallel do 
     deallocate (wkfepx,wkfeppx,wkfepy,wkfeppy,wkfepz,wkfeppz)
!-------------------------------------------------------------------------------
!  Obtain dfe/duex, dfe/duey, dfe/duez
!-------------------------------------------------------------------------------
!$OMP parallel private(ijk,iis,iie,ii,ijku,iEx,iEy,iEz,iBx,iBy,iBz) &
!$OMP          private(tempvb1,tempvb2,tempA,tempB,tempC,wkfex,wkfey,wkfez)
     allocate (wkfex(nuexyz),wkfey(nuexyz),wkfez(nuexyz))
!$OMP do 
     do k = 1, ncz_mpi
     do j = 1, ncy_mpi
     do i = 1, ncx_mpi
       ijk = i+(j-1)*ncx_mpi+(k-1)*ncxy_mpi
       iEx = ijk-1+ncEx
       iEy = ijk-1+ncEy
       iEz = ijk-1+ncEz
       iBx = ijk-1+ncBx
       iBy = ijk-1+ncBy
       iBz = ijk-1+ncBz
!
       iis = (ijk-1)*nuexyz+1
       iie = (ijk  )*nuexyz
!
       wkfex = 0.d0
       wkfey = 0.d0
       wkfez = 0.d0
!
       call gradu(f_total(iis:iie),wkfex,wkfey,wkfez,Buex,Cuex,huex, &
                  Buey,Cuey,huey,Buez,Cuez,huez,nuex,nuey,nuez)
!
       do 115 ku = 1, nuez
       do 115 ju = 1, nuey
       do 115 iu = 1, nuex
         if (iu == 1 .or. iu == nuex) goto 115
         if (nuey /= 1 .and. (ju == 1 .or. ju == nuey)) goto 115
         if (nuez /= 1 .and. (ku == 1 .or. ku == nuez)) goto 115
!
         ijku = iu+(ju-1)*nuex+(ku-1)*nuexy
         ii   = ijku+(ijk-1)*nuexyz
!
         tempvb1 = vey(ijku)*f_total(iBz)
         tempvb2 = vez(ijku)*f_total(iBy)
         tempA = sumABC(f_total(iEx),tempvb1,-tempvb2)*wkfex(ijku)
!
         tempvb1 = vez(ijku)*f_total(iBx)
         tempvb2 = vex(ijku)*f_total(iBz)
         tempB = sumABC(f_total(iEy),tempvb1,-tempvb2)*wkfey(ijku)
!
         tempvb1 = vex(ijku)*f_total(iBy)
         tempvb2 = vey(ijku)*f_total(iBx)
         tempC = sumABC(f_total(iEz),tempvb1,-tempvb2)*wkfez(ijku)
! 
         func_total(ii) = sumABCD(func_total(ii),tempA,tempB,tempC)
 115   continue
!
     enddo
     enddo
     enddo
!$OMP end do
     deallocate (wkfex,wkfey,wkfez)
!$OMP end parallel
!===============================================================================
!===============================================================================
!  Vlasov equation for ions
!-------------------------------------------------------------------------------
!  Obtain dfi/dx, dfi/dy, dfi/dx, d2fi/dx2, d2fi/dy2, d2fi/dz2
!-------------------------------------------------------------------------------
     allocate (wkfipx(nuixyz*ncxyz_mpi), wkfippx(nuixyz*ncxyz_mpi))
     allocate (wkfipy(nuixyz*ncxyz_mpi), wkfippy(nuixyz*ncxyz_mpi))
     allocate (wkfipz(nuixyz*ncxyz_mpi), wkfippz(nuixyz*ncxyz_mpi))
!
     wkfipx  = 0.d0
     wkfippx = 0.d0
     wkfipy  = 0.d0
     wkfippy = 0.d0
     wkfipz  = 0.d0
     wkfippz = 0.d0
!
     call grad6d_5th(ncx_mpi,ncy_mpi,ncz_mpi,nuix,nuiy,nuiz, &
                     f_total(ncfi),wkfipx,wkfipy,wkfipz, &
                     wkfippx,wkfippy,wkfippz,dx,dy,dz,ixp,iyp,izp)
!
!$OMP parallel do private(ii,ii1,ijk,ijku,tempA,tempB,tempC,tempD) 
     do 122 k = 1, ncz_mpi
     do 122 j = 1, ncy_mpi
     do 122 i = 1, ncx_mpi         
     do 122 ku = 1, nuiz
     do 122 ju = 1, nuiy
     do 122 iu = 1, nuix
       if (iu == 1 .or. iu == nuix) goto 122
       if (nuiy /= 1 .and. (ju == 1 .or. ju == nuiy)) goto 122
       if (nuiz /= 1 .and. (ku == 1 .or. ku == nuiz)) goto 122
!
       ijku = iu+(ju-1)*nuix+(ku-1)*nuixy
       ijk  = i+(j-1)*ncx_mpi+(k-1)*ncxy_mpi
       ii1  = ijku+(ijk-1)*nuixyz
       ii   = ii1+ncfi-1
!
       tempA = -vix(ijku)*wkfipx(ii1)
       tempB = -viy(ijku)*wkfipy(ii1)
       tempC = -viz(ijku)*wkfipz(ii1)
       tempD = +etafi(ijk)*sumABC(wkfippx(ii1),wkfippy(ii1),wkfippz(ii1))
       func_total(ii) = sumABCD(tempA,tempB,tempC,tempD)
 122 continue
!$OMP end parallel do 
     deallocate (wkfipx,wkfippx,wkfipy,wkfippy,wkfipz,wkfippz)
!-------------------------------------------------------------------------------
!  Obtain dfi/duix, dfi/duiy, dfi/duiz
!-------------------------------------------------------------------------------
!$OMP parallel private(ijk,iis,iie,ii,ijku,iEx,iEy,iEz,iBx,iBy,iBz) &
!$OMP          private(tempvb1,tempvb2,tempA,tempB,tempC,wkfix,wkfiy,wkfiz)
     allocate (wkfix(nuixyz),wkfiy(nuixyz),wkfiz(nuixyz))
!$OMP do 
     do k = 1, ncz_mpi
     do j = 1, ncy_mpi
     do i = 1, ncx_mpi
       ijk = i+(j-1)*ncx_mpi+(k-1)*ncxy_mpi
       iEx  = ijk-1+ncEx
       iEy  = ijk-1+ncEy
       iEz  = ijk-1+ncEz
       iBx  = ijk-1+ncBx
       iBy  = ijk-1+ncBy
       iBz  = ijk-1+ncBz
!
       iis = (ijk-1)*nuixyz+ncfi
       iie = (ijk  )*nuixyz+ncfi-1
!
       wkfix = 0.d0
       wkfiy = 0.d0
       wkfiz = 0.d0
!
       call gradu(f_total(iis:iie),wkfix,wkfiy,wkfiz,Buix,Cuix,huix, &
                  Buiy,Cuiy,huiy,Buiz,Cuiz,huiz,nuix,nuiy,nuiz)
!
       do 125 ku = 1, nuiz
       do 125 ju = 1, nuiy
       do 125 iu = 1, nuix
         if (iu == 1 .or. iu == nuix) goto 125
         if (nuiy /= 1 .and. (ju == 1 .or. ju == nuiy)) goto 125
         if (nuiz /= 1 .and. (ku == 1 .or. ku == nuiz)) goto 125
!
         ijku = iu+(ju-1)*nuix+(ku-1)*nuixy
         ii   = ijku+(ijk-1)*nuixyz+ncfi-1
!
         tempvb1 = viy(ijku)*f_total(iBz)
         tempvb2 = viz(ijku)*f_total(iBy)
         tempA = -sumABC(f_total(iEx),tempvb1,-tempvb2)*wkfix(ijku)/ami
!
         tempvb1 = viz(ijku)*f_total(iBx)
         tempvb2 = vix(ijku)*f_total(iBz)
         tempB = -sumABC(f_total(iEy),tempvb1,-tempvb2)*wkfiy(ijku)/ami
!!
         tempvb1 = vix(ijku)*f_total(iBy)
         tempvb2 = viy(ijku)*f_total(iBx)
         tempC = -sumABC(f_total(iEz),tempvb1,-tempvb2)*wkfiz(ijku)/ami
!! 
         func_total(ii) = sumABCD(func_total(ii),tempA,tempB,tempC)
 125   continue
!!
     enddo
     enddo
     enddo
!$OMP end do
     deallocate (wkfix,wkfiy,wkfiz)
!$OMP end parallel
!===============================================================================
!===============================================================================
!  Faraday's law
!-------------------------------------------------------------------------------
!  - curl E
     wkdAx = 0.d0
     wkdAy = 0.d0
     wkdAz = 0.d0
!
     call curl3d_5th(ncx_mpi,ncy_mpi,ncz_mpi, &
                     f_total(ncEx),f_total(ncEy),f_total(ncEz), &
                     wkdAx,wkdAy,wkdAz,dx,dy,dz,ixp,iyp,izp)
!
     func_total(ncBx:ncBx_end) = -wkdAx(1:ncxyz_mpi)
     func_total(ncBy:ncBy_end) = -wkdAy(1:ncxyz_mpi)
     func_total(ncBz:ncBz_end) = -wkdAz(1:ncxyz_mpi) 
!-------------------------------------------------------------------------------
!  etaB * (ddBx/dxx + ddBy/dyy + ddBz/dzz)
     wkdAx = 0.d0
     wkdAy = 0.d0
     wkdAz = 0.d0
!
     call divgradx3d_5th(ncx_mpi,ncy_mpi,ncz_mpi,f_total(ncBx),wkdAx,dx,ixp)
     call divgrady3d_5th(ncx_mpi,ncy_mpi,ncz_mpi,f_total(ncBy),wkdAy,dy,iyp)
     call divgradz3d_5th(ncx_mpi,ncy_mpi,ncz_mpi,f_total(ncBz),wkdAz,dz,izp)
!
!$OMP parallel do private(iBx,iBy,iBz,tempA,tempB,tempC)
     do ii = 1, ncxyz_mpi      
       iBx  = ii-1+ncBx
       iBy  = ii-1+ncBy
       iBz  = ii-1+ncBz
!
       tempA = etaB(ii)*wkdAx(ii)
       tempB = etaB(ii)*wkdAy(ii)
       tempC = etaB(ii)*wkdAz(ii)
       func_total(iBx) = sumAB(func_total(iBx),tempA)
       func_total(iBy) = sumAB(func_total(iBy),tempB)
       func_total(iBz) = sumAB(func_total(iBz),tempC)     
     enddo
!$OMP end parallel do
!===============================================================================
!===============================================================================
!  Ampere's law
!-------------------------------------------------------------------------------
!  curl B
     call curl3d_5th(ncx_mpi,ncy_mpi,ncz_mpi, &
                     f_total(ncBx),f_total(ncBy),f_total(ncBz), &
                     func_total(ncEx),func_total(ncEy),func_total(ncEz), &
                     dx,dy,dz,ixp,iyp,izp)
!
     C2 = C*C 
!$OMP parallel do private(iEx,iEy,iEz)
     do ii = 1, ncxyz_mpi      
       iEx  = ii-1+ncEx
       iEy  = ii-1+ncEy
       iEz  = ii-1+ncEz
!
       func_total(iEx) = C2*func_total(iEx)
       func_total(iEy) = C2*func_total(iEy)
       func_total(iEz) = C2*func_total(iEz)     
     enddo
!$OMP end parallel do
!-------------------------------------------------------------------------------
!  obtain current for electron
!
     call current(ncx_mpi,ncy_mpi,ncz_mpi,nuex,nuey,nuez,f_total(ncfe),vex,vey,vez, &
                  Buex,Cuex,huex,Buey,Cuey,huey,Buez,Cuez,huez,wkJex,wkJey,wkJez)
!     wkJey = 0.d0
!     wkJez = 0.d0
!-------------------------------------------------------------------------------
!  obtain current for ion
!
     call current(ncx_mpi,ncy_mpi,ncz_mpi,nuix,nuiy,nuiz,f_total(ncfi),vix,viy,viz, &
                  Buix,Cuix,huix,Buiy,Cuiy,huiy,Buiz,Cuiz,huiz,wkJix,wkJiy,wkJiz)
!     wkJiy = 0.d0
!     wkJiz = 0.d0
!-------------------------------------------------------------------------------
!$OMP parallel do private(iEx,iEy,iEz)
     do ii = 1, ncxyz_mpi      
       iEx  = ii-1+ncEx
       iEy  = ii-1+ncEy
       iEz  = ii-1+ncEz
!
       func_total(iEx) = sumABCD(func_total(iEx),-wkJix(ii),wkJex(ii),+wkJx0(ii))
       func_total(iEy) = sumABCD(func_total(iEy),-wkJiy(ii),wkJey(ii),+wkJy0(ii))
       func_total(iEz) = sumABCD(func_total(iEz),-wkJiz(ii),wkJez(ii),+wkJz0(ii))     
     enddo
!$OMP end parallel do
!-------------------------------------------------------------------------------
   end subroutine Vlasov_EM
!===============================================================================
!*******************************************************************************
!===============================================================================
   subroutine Vlasov_ES(time,f_total,func_total,Nt)
     use GRID6D
     use vector3D
     use Vlasov2fluid
     use omp_lib
     implicit double precision (a-h,o-z)
     dimension :: f_total(Nt), func_total(Nt)
     real(8) :: tempA, tempB, tempC
     func_total(1:NT) = 0.d0
!===============================================================================
!  Vlasov equation for electrons
!-------------------------------------------------------------------------------
!  Obtain dfe/dx, d2fe/dx2
!-------------------------------------------------------------------------------
     allocate (wkfepx(nuexyz*ncxyz_mpi), wkfeppx(nuexyz*ncxyz_mpi))
     allocate (wkfepy(nuexyz*ncxyz_mpi), wkfeppy(nuexyz*ncxyz_mpi))
     allocate (wkfepz(nuexyz*ncxyz_mpi), wkfeppz(nuexyz*ncxyz_mpi))
!
     wkfepx  = 0.d0
     wkfeppx = 0.d0
     wkfepy  = 0.d0
     wkfeppy = 0.d0
     wkfepz  = 0.d0
     wkfeppz = 0.d0
!
     call grad6d_WENO(ncx_mpi,ncy_mpi,ncz_mpi,nuex,nuey,nuez, &
                     f_total(ncfe),wkfepx,wkfepy,wkfepz, &
                     wkfeppx,wkfeppy,wkfeppz,dx,dy,dz,ixp,iyp,izp)
!
!$OMP parallel do private(ii,ijk,ijku,tempA,tempB,tempC,tempD) 
     do 112 k = 1, ncz_mpi
     do 112 j = 1, ncy_mpi
     do 112 i = 1, ncx_mpi         
     do 112 ku = 1, nuez
     do 112 ju = 1, nuey
     do 112 iu = 1, nuex
       if (iu == 1 .or. iu == nuex) goto 112
       if (nuey /= 1 .and. (ju == 1 .or. ju == nuey)) goto 112
       if (nuez /= 1 .and. (ku == 1 .or. ku == nuez)) goto 112
!
       ijku = iu+(ju-1)*nuex+(ku-1)*nuexy
       ijk  = i+(j-1)*ncx_mpi+(k-1)*ncxy_mpi
       ii   = ijku+(ijk-1)*nuexyz
!
       tempA = -vex(ijku)*wkfepx(ii)
       tempB = -vey(ijku)*wkfepy(ii)
       tempC = -vez(ijku)*wkfepz(ii)
       tempD = +etafe(ijk)*sumABC(wkfeppx(ii),wkfeppy(ii),wkfeppz(ii))
       func_total(ii) = sumABCD(tempA,tempB,tempC,tempD)
 112 continue
!$OMP end parallel do 
     deallocate (wkfepx,wkfeppx,wkfepy,wkfeppy,wkfepz,wkfeppz)
!-------------------------------------------------------------------------------
!  Obtain dfe/duex
!-------------------------------------------------------------------------------
!$OMP parallel private(ijk,iis,iie,ii,ijku,iEx,iEy,iEz) &
!$OMP          private(tempA,tempB,tempC,wkfex,wkfey,wkfez)
     allocate (wkfex(nuexyz),wkfey(nuexyz),wkfez(nuexyz))
!$OMP do 
     do k = 1, ncz_mpi
     do j = 1, ncy_mpi
     do i = 1, ncx_mpi
       ijk = i+(j-1)*ncx_mpi+(k-1)*ncxy_mpi
       iEx = ijk-1+ncEx
       iEy = ijk-1+ncEy
       iEz = ijk-1+ncEz
!
       iis = (ijk-1)*nuexyz+1
       iie = (ijk  )*nuexyz
!
       wkfex = 0.d0
       wkfey = 0.d0
       wkfez = 0.d0
!
       call gradu(f_total(iis:iie),wkfex,wkfey,wkfez,Buex,Cuex,huex, &
                  Buey,Cuey,huey,Buez,Cuez,huez,nuex,nuey,nuez)
!
       do 115 ku = 1, nuez
       do 115 ju = 1, nuey
       do 115 iu = 1, nuex
         if (iu == 1 .or. iu == nuex) goto 115
         if (nuey /= 1 .and. (ju == 1 .or. ju == nuey)) goto 115
         if (nuez /= 1 .and. (ku == 1 .or. ku == nuez)) goto 115
!         
         ijku = iu+(ju-1)*nuex+(ku-1)*nuexy
         ii   = ijku+(ijk-1)*nuexyz
!
         tempA = f_total(iEx)*wkfex(ijku)
         tempB = f_total(iEy)*wkfey(ijku)
         tempC = f_total(iEz)*wkfez(ijku)
! 
         func_total(ii) = sumABCD(func_total(ii),tempA,tempB,tempC)
 115   continue
!
     enddo
     enddo
     enddo
!$OMP end do
     deallocate (wkfex,wkfey,wkfez)
!$OMP end parallel
!===============================================================================
!===============================================================================
!  Vlasov equation for ions
!-------------------------------------------------------------------------------
!  Obtain dfi/dx, d2fi/dx2
!-------------------------------------------------------------------------------
     allocate (wkfipx(nuixyz*ncxyz_mpi), wkfippx(nuixyz*ncxyz_mpi))
     allocate (wkfipy(nuixyz*ncxyz_mpi), wkfippy(nuixyz*ncxyz_mpi))
     allocate (wkfipz(nuixyz*ncxyz_mpi), wkfippz(nuixyz*ncxyz_mpi))
!
     wkfipx  = 0.d0
     wkfippx = 0.d0
     wkfipy  = 0.d0
     wkfippy = 0.d0
     wkfipz  = 0.d0
     wkfippz = 0.d0
!
     call grad6d_WENO(ncx_mpi,ncy_mpi,ncz_mpi,nuix,nuiy,nuiz, &
                     f_total(ncfi),wkfipx,wkfipy,wkfipz, &
                     wkfippx,wkfippy,wkfippz,dx,dy,dz,ixp,iyp,izp)
!
!$OMP parallel do private(ii,ii1,ijk,ijku,tempA,tempB,tempC,tempD) 
     do 122 k = 1, ncz_mpi
     do 122 j = 1, ncy_mpi
     do 122 i = 1, ncx_mpi         
     do 122 ku = 1, nuiz
     do 122 ju = 1, nuiy
     do 122 iu = 1, nuix
       if (iu == 1 .or. iu == nuix) goto 122
       if (nuiy /= 1 .and. (ju == 1 .or. ju == nuiy)) goto 122
       if (nuiz /= 1 .and. (ku == 1 .or. ku == nuiz)) goto 122
!
       ijku = iu+(ju-1)*nuix+(ku-1)*nuixy
       ijk  = i+(j-1)*ncx_mpi+(k-1)*ncxy_mpi
       ii1  = ijku+(ijk-1)*nuixyz
       ii   = ii1+ncfi-1
!
       tempA = -vix(ijku)*wkfipx(ii1)
       tempB = -viy(ijku)*wkfipy(ii1)
       tempC = -viz(ijku)*wkfipz(ii1)
       tempD = +etafi(ijk)*sumABC(wkfippx(ii1),wkfippy(ii1),wkfippz(ii1))
       func_total(ii) = sumABCD(tempA,tempB,tempC,tempD)
 122 continue
!$OMP end parallel do 
     deallocate (wkfipx,wkfippx,wkfipy,wkfippy,wkfipz,wkfippz)
!-------------------------------------------------------------------------------
!  Obtain dfi/duix, dfi/duiy, dfi/duiz
!-------------------------------------------------------------------------------
!$OMP parallel private(ijk,iis,iie,ii,ijku,iEx,iEy,iEz) &
!$OMP          private(tempA,tempB,tempC,wkfix,wkfiy,wkfiz)
     allocate (wkfix(nuixyz),wkfiy(nuixyz),wkfiz(nuixyz))
!$OMP do 
     do k = 1, ncz_mpi
     do j = 1, ncy_mpi
     do i = 1, ncx_mpi
       ijk = i+(j-1)*ncx_mpi+(k-1)*ncxy_mpi
       iEx  = ijk-1+ncEx
       iEy  = ijk-1+ncEy
       iEz  = ijk-1+ncEz
!
       iis = (ijk-1)*nuixyz+ncfi
       iie = (ijk  )*nuixyz+ncfi-1
!
       wkfix = 0.d0
       wkfiy = 0.d0
       wkfiz = 0.d0
!
       call gradu(f_total(iis:iie),wkfix,wkfiy,wkfiz,Buix,Cuix,huix, &
                  Buiy,Cuiy,huiy,Buiz,Cuiz,huiz,nuix,nuiy,nuiz)     
!
       do 125 ku = 1, nuiz
       do 125 ju = 1, nuiy
       do 125 iu = 1, nuix
         if (iu == 1 .or. iu == nuix) goto 125
         if (nuiy /= 1 .and. (ju == 1 .or. ju == nuiy)) goto 125
         if (nuiz /= 1 .and. (ku == 1 .or. ku == nuiz)) goto 125
!
         ijku = iu+(ju-1)*nuix+(ku-1)*nuixy
         ii   = ijku+(ijk-1)*nuixyz+ncfi-1
!
         tempA = -f_total(iEx)*wkfix(ijku)/ami
         tempB = -f_total(iEy)*wkfiy(ijku)/ami
         tempC = -f_total(iEz)*wkfiz(ijku)/ami
! 
         func_total(ii) = sumABCD(func_total(ii),tempA,tempB,tempC)
 125   continue
!
     enddo
     enddo
     enddo
!$OMP end do
     deallocate (wkfix,wkfiy,wkfiz)
!$OMP end parallel
!===============================================================================
!===============================================================================
!  Ampere's law
!-------------------------------------------------------------------------------
!  obtain current for electron
!
     call current(ncx_mpi,ncy_mpi,ncz_mpi,nuex,nuey,nuez,f_total(ncfe),vex,vey,vez, &
                  Buex,Cuex,huex,Buey,Cuey,huey,Buez,Cuez,huez,wkJex,wkJey,wkJez)
!-------------------------------------------------------------------------------
!  obtain current for ion
!
     call current(ncx_mpi,ncy_mpi,ncz_mpi,nuix,nuiy,nuiz,f_total(ncfi),vix,viy,viz, &
                  Buix,Cuix,huix,Buiy,Cuiy,huiy,Buiz,Cuiz,huiz,wkJix,wkJiy,wkJiz)
!-------------------------------------------------------------------------------
!$OMP parallel do private(iEx,iEy,iEz)
     do ii = 1, ncxyz_mpi      
       iEx  = ii-1+ncEx
       iEy  = ii-1+ncEy
       iEz  = ii-1+ncEz
!
       func_total(iEx) = sumABC(-wkJix(ii),wkJex(ii),+wkJx0(ii))
       func_total(iEy) = sumABC(-wkJiy(ii),wkJey(ii),+wkJy0(ii))
       func_total(iEz) = sumABC(-wkJiz(ii),wkJez(ii),+wkJz0(ii))     
     enddo
!$OMP end parallel do
!-------------------------------------------------------------------------------
   end subroutine Vlasov_ES
!===============================================================================
!===============================================================================
!*******************************************************************************
!===============================================================================
   subroutine Vlasov_ESim(time,f_total,func_total,Nt)
!
!  ion : uniform, motionless background
!
     use GRID6D
     use vector3D
     use Vlasov2fluid
     use CubicSpline
     use omp_lib
     implicit double precision (a-h,o-z)
     dimension :: f_total(Nt), func_total(Nt)
     real(8) :: tempA, tempB, tempC
     func_total(1:NT) = 0.d0
!===============================================================================
!  Vlasov equation for electrons
!-------------------------------------------------------------------------------
!  Obtain dfe/dx, d2fe/dx2
!-------------------------------------------------------------------------------
!     call printThreads()
     allocate (wkfepx(nuexyz*ncxyz_mpi), wkfeppx(nuexyz*ncxyz_mpi))
     allocate (wkfepy(nuexyz*ncxyz_mpi), wkfeppy(nuexyz*ncxyz_mpi))
     allocate (wkfepz(nuexyz*ncxyz_mpi), wkfeppz(nuexyz*ncxyz_mpi))
!
     wkfepx  = 0.d0
     wkfeppx = 0.d0
     wkfepy  = 0.d0
     wkfeppy = 0.d0
     wkfepz  = 0.d0
     wkfeppz = 0.d0
!
!     call grad6d_5th(ncx_mpi,ncy_mpi,ncz_mpi,nuex,nuey,nuez, &
!                     f_total(ncfe),wkfepx,wkfepy,wkfepz, &
!                     wkfeppx,wkfeppy,wkfeppz,dx,dy,dz,ixp,iyp,izp)
     call grad6d_WENO(ncx_mpi,ncy_mpi,ncz_mpi,nuex,nuey,nuez, &
                     f_total(ncfe),wkfepx,wkfepy,wkfepz, &
                     wkfeppx,wkfeppy,wkfeppz,dx,dy,dz,ixp,iyp,izp)
!
!$OMP parallel do private(ii,ijk,ijku,tempA,tempB,tempC,tempD) 
     do 112 k = 1, ncz_mpi
     do 112 j = 1, ncy_mpi
     do 112 i = 1, ncx_mpi         
     do 112 ku = 1, nuez
     do 112 ju = 1, nuey
     do 112 iu = 1, nuex
       if (iu == 1 .or. iu == nuex) goto 112
       if (nuey /= 1 .and. (ju == 1 .or. ju == nuey)) goto 112
       if (nuez /= 1 .and. (ku == 1 .or. ku == nuez)) goto 112
!
       ijku = iu+(ju-1)*nuex+(ku-1)*nuexy
       ijk  = i+(j-1)*ncx_mpi+(k-1)*ncxy_mpi
       ii   = ijku+(ijk-1)*nuexyz
!
       tempA = -vex(ijku)*wkfepx(ii)
       tempB = -vey(ijku)*wkfepy(ii)
       tempC = -vez(ijku)*wkfepz(ii)
!       m_k = 0.8d0
!       tempD = dexp(-1.d0*(wkfepx(ii)**2)/m_k**2)*wkfeppx(ii)
!       tempD = +etafe(ijk)*tempD*sumAB(-2.d0*wkfepx(ii)**2/m_k**2,1.d0)
       tempD = 0.d0!+etafe(ijk)*wkfeppx(ii)
       func_total(ii) = sumABCD(tempA,tempB,tempC,tempD)
 112 continue
!$OMP end parallel do
     deallocate (wkfepx,wkfeppx,wkfepy,wkfeppy,wkfepz,wkfeppz)
!-------------------------------------------------------------------------------
!  Obtain dfe/duex
!-------------------------------------------------------------------------------
!$OMP parallel private(ijk,iis,iie,ii,ijku,iEx,iEy,iEz) &
!$OMP          private(tempA,tempB,tempC,tempD,du,wkfex,wkfey,wkfez)
     allocate (wkfex(nuexyz),wkfey(nuexyz),wkfez(nuexyz))
!$OMP do 
     do k = 1, ncz_mpi
     do j = 1, ncy_mpi
     do i = 1, ncx_mpi
       ijk = i+(j-1)*ncx_mpi+(k-1)*ncxy_mpi
       iEx = ijk-1+ncEx
       iEy = ijk-1+ncEy
       iEz = ijk-1+ncEz
!
       iis = (ijk-1)*nuexyz+1
       iie = (ijk  )*nuexyz
       du = vex(2)-vex(1)
!
       wkfex = 0.d0
       wkfey = 0.d0
       wkfez = 0.d0
!
       !call gradu(f_total(iis:iie),wkfex,wkfey,wkfez,Buex,Cuex,huex, &
       !           Buey,Cuey,huey,Buez,Cuez,huez,nuex,nuey,nuez)
       call gradux_WENO5(f_total(iis:iie),wkfex,du,nuex,nuey,nuez,-1.d0*f_total(iEx))
       !call FP_ux_CS(f_total(iis:iie),wkfey,Buex,Cuex,huex,nuex,nuey,nuez,vex)
       !call FP_ux_CWENO(f_total(iis:iie),wkfey,du,nuex,nuey,nuez,vex)
!
       do 115 ku = 1, nuez
       do 115 ju = 1, nuey
       do 115 iu = 1, nuex
         if (iu == 1 .or. iu == nuex) goto 115
         if (nuey /= 1 .and. (ju == 1 .or. ju == nuey)) goto 115
         if (nuez /= 1 .and. (ku == 1 .or. ku == nuez)) goto 115
!         
         ijku = iu+(ju-1)*nuex+(ku-1)*nuexy
         ii   = ijku+(ijk-1)*nuexyz
!
         tempA = f_total(iEx)*wkfex(ijku)
         tempB = 0!f_total(iEy)*wkfey(ijku)
         !tempC = 0!f_total(iEz)*wkfez(ijku)
         tempC = 0!1e-3*wkfey(ijku)
! 
         func_total(ii) = sumABCD(func_total(ii),tempA,tempB,tempC)
 115   continue
!
     enddo
     enddo
     enddo
!$OMP end do
     deallocate (wkfex,wkfey,wkfez)
!$OMP end parallel
!===============================================================================
!===============================================================================
!  Ampere's law
!-------------------------------------------------------------------------------
!  obtain current for electron
!
!     call current(ncx_mpi,ncy_mpi,ncz_mpi,nuex,nuey,nuez,f_total(ncfe),vex,vey,vez, &
!                  Buex,Cuex,huex,Buey,Cuey,huey,Buez,Cuez,huez,wkJex,wkJey,wkJez)
!-------------------------------------------------------------------------------
!     avJx0 = 0.d0
!     do ii = 1, ncxyz_mpi      
!       avJx0 = avJx0 + wkJex(ii)
!     enddo
!     avJx0_global = 0.d0
!     if (MyID .eq. 0) then
!          avJx0 = avJx0 - wkJex(1)
!     endif
!     call MPI_ALLREDUCE(avJx0,avJx0_global,1,MPI_REAL8,MPI_SUM,MPI_WORLD,IERR)
!     avJx0_global = avJx0_global/(ncx-1)
!  obtain current for ion
!
!     call current(ncx_mpi,ncy_mpi,ncz_mpi,nuix,nuiy,nuiz,f_total(ncfi),vix,viy,viz, &
!                  Buix,Cuix,huix,Buiy,Cuiy,huiy,Buiz,Cuiz,huiz,wkJix,wkJiy,wkJiz)
     wkJix = 0.d0
     wkJiy = 0.d0
     wkJiz = 0.d0
!-------------------------------------------------------------------------------
!$OMP parallel do private(iEx,iEy,iEz)
     do ii = 1, ncxyz_mpi      
       iEx  = ii-1+ncEx
       iEy  = ii-1+ncEy
       iEz  = ii-1+ncEz
!
       func_total(iEx) = sumAB(-wkJix(ii),wkJex(ii))
     enddo
!$OMP end parallel do
!-------------------------------------------------------------------------------
   end subroutine Vlasov_ESim
!===============================================================================
!===============================================================================
