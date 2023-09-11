module CubicSpline
      use coef 
   contains
!===============================================================================
!
!  FCS_PCS_YP.f : (geth & getFP0)
!                 (FCSYP & FCSYP0 & FCSYP1)
!                 (PCSYP & PCSYP0 & PCSYP1)
!                 (YP2YPP & getAB & getArea & getFF)
!
!===============================================================================
!--SUB. geth
      subroutine geth(Z,h,N)
!
!  INPUT ARRAY   : Z
!  OUTPUT ARRAY  : h  (to be used in cubic spline)
!  
      implicit double precision (A-H,O-Z)
      dimension Z(1)
      dimension h(1)
!
      if (N .le. 2) then
        print *,'Program geth stop ! N must be greater than 2'
        stop
      endif
!
      do i = 1, N-1
        h(i) = Z(i+1)-Z(i)
      enddo
      h(N) = h(1)  ! prepared for periodic cases
      end subroutine geth
!-------------------------------------------------------------------------------
!--SUB. getFP0
      subroutine getFP0(F,h,FP0,N)
!
!  INPUT ARRAY   : F,h
!  OUTPUT ARRAY  : FP0  (to be used in cubic spline)
!  
      implicit double precision (A-H,O-Z)
      dimension F(1),h(1)
      dimension FP0(1)
!
      if (N .le. 2) then
        print *,'Program getFP0 stop ! N must be greater than 2'
        stop
      endif
!
      do i = 1, N-1
        FP0(i) = DIFAB(F(i+1),F(i))/h(i)
!       FP0(i) = (F(i+1)-F(i))/h(i)
      enddo
      FP0(N) = FP0(1) ! prepared for periodic cases
!
      end subroutine getFP0
!===============================================================================
!--SUB. FCSYP
      subroutine FCSYP(F,Z,FP0,h,FP,N,B,C,R)
!
!  Using cubic spline to solve FP=df/dz
!    with fixed Boundary Condition of F and FP
!
!  INPUT DATA     : FP(1), FP(N)
!  INPUT ARRAYS   : F, Z
!  WORKING ARRAYS : B, C, R
!  OUTPUT ARRAYS  : FP0, h, FP
!  
      implicit double precision (A-H,O-Z)
      dimension Z(1), F(1), FP(1)
      dimension h(1), FP0(1)
      dimension B(1), C(1), R(1)
!
      if (N .le. 2) then
        print *,'Program FCSYP stop ! N must be greater than 2'
        stop
      endif
!
      call geth(Z,h,N)
!
      call FCSYP0(h,B,C,N)
!
      call getFP0(F,h,FP0,N)
!
      call FCSYP1(FP0,B,C,FP,N,R)
!
      end subroutine FCSYP
!-------------------------------------------------------------------------------
!--SUB. FCSYP0
      subroutine FCSYP0(h,B,C,N)
!
!  Using cubic spline to solve FP=df/dz in two steps
!    with fixed boundary condition of F and FP.
!  This is the first step to set the constants' table
!
!  INPUT ARRAY  : h
!  OUTPUT ARRAYS : B, C
!
      implicit double precision (A-H,O-Z)
      dimension B(1), C(1)
      dimension h(1)
      data A1 /1.d0/
!
!  Set the tridiagonal matrix coefficient 
!    for fixed Boundary Condition of FP at z=z1 and z=zn
!
      do i = 2, N-1
        C(i) = h(i-1)/h(i)
        B(i) = 2*DIFAB(A1,-C(i))
!       B(i) = 2.d0+2*C(i)
      enddo
!
      i2 = 2
      call TRID0(B(i2),C(i2),N-2)
!
      end subroutine FCSYP0
!-------------------------------------------------------------------------------
!--SUB. FCSYP1
      subroutine FCSYP1(FP0,B,C,FP,N,R)
!
!  Using cubic spline to solve FP=df/dz in two steps
!    with fixed boundary condition of F and FP.
!  This is the second step to get FP=df/dz
!
!  INPUT DATA    : FP(1), FP(N)
!  INPUT ARRAYS  : FP0, B, C
!  WORKING ARRAY : R
!  OUTPUT ARRAY  : FP
!
      implicit double precision (A-H,O-Z)
      dimension B(1), C(1), R(1)
      dimension FP0(1), FP(1)
!
!  Set the right hand side column
!    for fixed Boundary Condition of FP at z=z1 and z=zn
!
      do i = 2, N-1
        temp = -FP0(i)*C(i)
        R(i) = 3*DIFAB(FP0(i-1),temp)
!        R(i) = 3*FP0(i-1)+3*FP0(i)*C(i)
      enddo
!
      i2 = 2
      R(i2)  = DIFAB(R(i2),FP(1))
      temp   = FP(N)*C(N-1)
      R(N-1) = DIFAB(R(N-1),temp)
!     R(2)   = R(2)-FP(1)
!     R(N-1) = R(N-1)-FP(N)*C(N-1)
!
      call TRID1(B(i2),C(i2),R(i2),FP(i2),N-2)
!
      end subroutine FCSYP1
!===============================================================================
!--SUB. PCSYP
      subroutine PCSYP(F,Z,FP0,h,FP,N,B,C,D,E,R)
!
!  Using cubic spline to solve FP=df/dz
!    with Periodic Boundary Condition : F(1)=F(N),
!    FP(1)=FP(N), FPP(1)=FPP(N)
!
!  INPUT ARRAYS   : F, Z
!  WORKING ARRAYS : B, C, D, E, R
!  OUTPUT ARRAYS  : FP, FP0, h
!  
      implicit double precision (A-H,O-Z)
      dimension Z(1), F(1), FP(1)
      dimension h(1), FP0(1)
      dimension B(1), C(1), D(1), E(1), R(1)
!
      if (N .le. 2) then
        print *,'Program PCSYP stop ! N must be greater than 2'
        stop
      endif
!
      call geth(Z,h,N)
!
      call PCSYP0(h,B,C,D,E,N)
!
      call getFP0(F,h,FP0,N)
!
      call PCSYP1(FP0,B,C,D,E,FP,N,R)
!
      end subroutine PCSYP
!-------------------------------------------------------------------------------
!--SUB. PCSYP0
      subroutine PCSYP0(h,B,C,D,E,N)
!
!  Using cubic spline to solve FP=df/dz
!    with Periodic Boundary Condition : F(1)=F(N),
!    FP(1)=FP(N), FPP(1)=FPP(N)
!
!  INPUT ARRAY     : h
!  OUTPUT ARRAYS   : B, C, D, E
!
      implicit double precision (A-H,O-Z)
      dimension B(1), C(1), D(1), E(1)
      dimension h(1)
      data A1 /1.d0/
!
!  Set the tridiagonal matrix coefficient
!    with Periodic Boundary Condition
!
      C(1) = h(N-1)/h(1)
      B(1) = 2*DIFAB(A1,-C(1))
!     B(1) = 2.d0+2*C(1)

      do i = 2, N-1
        C(i) = h(i-1)/h(i)
        B(i) = 2*DIFAB(A1,-C(i))
!       B(i) = 2.d0+2*C(i)
      enddo
!
      call TRIX0(B,C,D,E,N-1)
!
      end subroutine PCSYP0
!-------------------------------------------------------------------------------
!--SUB. PCSYP1
      subroutine PCSYP1(FP0,B,C,D,E,FP,N,R)
!
!  Using cubic spline to solve FP=df/dz in two steps
!    with periodic boundary condition of F and FP.
!  This is the second step to get FP=df/dz
!
!  INPUT ARRAYS  : FP0, B, C, D, E
!  WORKING ARRAY : R
!  OUTPUT ARRAY  : FP
!
      implicit double precision (A-H,O-Z)
      dimension B(1), C(1), D(1), E(1), R(1)
      dimension FP0(1), FP(1)
!
!  Set the right hand side column
!    for periodic Boundary Condition of FP
!
      temp = -FP0(1)*C(1)
      R(1) = 3*DIFAB(FP0(N-1),temp)
!     R(1) = 3*FP0(N-1)+3*FP0(1)*C(1)

      do i = 2, N-1
        temp = -FP0(i)*C(i)
        R(i) = 3*DIFAB(FP0(i-1),temp)
!       R(i) = 3*FP0(i-1)+3*FP0(i)*C(i)
      enddo
!
      call TRIX1(B,C,D,E,R,FP,N-1)
!
      end subroutine PCSYP1
!===============================================================================
!--SUB. YP2YPP
      subroutine YP2YPP(FP,FP0,h,FPP,N)
!
!  INPUT ARRAYS  : FP, FP0, h
!  OUTPUT ARRAY  : FPP
!
      implicit double precision (A-H,O-Z)
      dimension FP(1), FP0(1), h(1), FPP(1)
!
      do i = 1, N-1
        temp   = 4*DIFAB(FP0(i),FP(i))
        temp1  = 2*DIFAB(FP(i+1),FP0(i))
        FPP(i) = DIFAB(temp,temp1)/h(i)
!       FPP(i) = (6*FP0(i)-4*FP(i)-2*FP(i+1))/h(i)
      enddo
      temp1  = 2*DIFAB(FP(N-1),FP0(N-1))
      temp   = 4*DIFAB(FP0(N-1),FP(N))
      FPP(N) = DIFAB(temp1,temp)/h(N-1)
!     FPP(N) = (2*FP(N-1)+4*FP(N)-6*FP0(N-1))/h(N-1)
!
      end subroutine YP2YPP
!-------------------------------------------------------------------------------
!--SUB. getAB
      subroutine getAB(FP,FP0,h,N,Alpha,Beta)
!
!  INPUT ARRAYS  : FP, FP0, h
!  OUTPUT ARRAYS : Alpha, Beta
!
      implicit double precision (A-H,O-Z)
      dimension FP(1), FP0(1), h(1), Alpha(1), Beta(1)
!
      do i = 1, N-1
        temp     = DIFAB(FP(i+1),FP0(i))
        temp1    = DIFAB(FP0(i),FP(i))
        Alpha(i) = DIFAB(temp, temp1)*h(i)
        Beta(i)  = DIFAB(FP0(i),FP(i))*h(i)
!       Alpha(i) = (FP(i+1)+FP(i)-2*FP0(i))*h(i)
!       Beta(i)  = (FP0(i)-FP(i))*h(i)
      enddo
!
      end subroutine getAB
!-------------------------------------------------------------------------------
!--SUB. getArea
      subroutine getArea(F,h,Alpha,Beta,Area, AreaS,N)
!
!  INPUT ARRAYS : F, h, Alpha, Beta
!  OUTPUT ARRAY : Area
!  OUTPUT DATA  : AreaS
!
      implicit double precision (A-H,O-Z)
      dimension F(1), h(1), Alpha(1), Beta(1), Area(1)
!
      AreaS = 0.d0
!
      do i = 1, N-1
        temp    = 3*DIFAB(F(i),-F(i+1))
        temp    = 2*DIFAB(temp, Beta(i))
        Area(i) = DIFAB(temp, Alpha(i))*h(i)/12.d0
!       Area(i) = ((F(i)+F(i+1))/2.d0-Alpha(i)/12.d0-Beta(i)/6.d0)*h(i)
        AreaS   = AreaS+Area(i)
      enddo
!
      end subroutine getArea
!-------------------------------------------------------------------------------
!--SUB. getFF
      subroutine getFF(F, Z, Alpha, Beta, ZZ, FF, N, M)
!
!  INPUT ARRAYS : F, Z, Alpha, Beta, ZZ
!  OUTPUT ARRAY : FF
!
      implicit double precision (A-H,O-Z)
      dimension F(1), Z(1), Alpha(1), Beta(1)
      dimension FF(M), ZZ(M)
!
      do j=1,M
      do i = 2, N
      IF(ZZ(j).LT.Z(i)) then
      ir=i
      il=i-1
      h=Z(ir)-Z(il)
      FF(j)=(F(ir)*(ZZ(j)-Z(il))-F(il)*(ZZ(j)-Z(ir)))/h  &
            +(Alpha(il)*(ZZ(j)-Z(il))/h+Beta(il))*(ZZ(j)-Z(il))*(ZZ(j)-Z(ir))/h**2
      go to 1
      endif
      enddo
    1 continue
      IF(ZZ(j).GE.Z(N)) FF(j)=F(N)
      enddo
!
      end subroutine getFF
!===============================================================================
!
! TRID_TRIX.f : ( TRID0 & TRID1 )
!               ( TRIX0 & TRIX1 )
!               ( GETRERR & DIFAB )
!
!===============================================================================
!--SUB. TRID0
      subroutine TRID0(B,C,N)
!
!  Solving MX=R, where M is a tridiagonal matrix M=(ABC)
!    with A=1 for all i>1.  M is a N*N matrix
!  This is the step to transform the tridiagonal matrix to an
!    upper triangular matrix.
!
!  INPUT ARRAYS  : B, C for tridiagonal matrix
!  WORKING VARIABLES : temp, temp1,...
!  OUTPUT ARRAY : B(new), for upper triangular matrix
!
      implicit double precision (A-H,O-Z)
      dimension B(1), C(1)
!
        if (B(1) .eq. 0.d0) then
          print *,'Program TRID0 stop! B(1) is equal to zero.'
          stop
        endif
!
      do i = 2, N
        temp = C(i-1)/B(i-1)
        temp1 = DIFAB(B(i),temp)
!        temp1 = B(i)-temp
!
! IN TRID, B(i), CANNOT BE ZERO  ?
!
        if (temp1 .eq. 0.d0) then
          print *,'Program TRID0 stop! Denominator is too small.'
          print *,' at I=',i
          print *,' B(I), B(I-1), C(I-1) =', B(i), B(i-1), C(i-1)
          stop
        endif
!
        B(i) = temp1
      enddo
!
      end subroutine TRID0
!-------------------------------------------------------------------------------
!--SUB. TRID1
      subroutine TRID1(B,C,R,X,N)
!
!  Solving MX=R, where M is a tridiagonal matrix M=(ABC)
!    with A=1 for all i>1.  M is a N*N matrix
!  This is the second step to solve X
!
!  INPUT ARRAYS  : B(new), C, R
!  OUTPUT ARRAY : X
!
      implicit double precision (A-H,O-Z)
      dimension B(1), C(1), R(1)
      dimension X(1)
!
!  When the tridiagonal matrix is transformed to an upper 
!    triangular matrix,
!    the right hand side column need to change its value.
!
      do i = 2, N
        temp = R(i-1)/B(i-1)
        temp1 = DIFAB(R(i),temp)
        R(i) = temp1
!        R(i) = R(i)-R(i-1)/B(i-1)
      enddo
!
!  Solve the upper triangular martix to get X
!
      X(N) = R(N)/B(N)
!
      do i = N-1, 1, -1
        temp = C(i)*X(i+1)
        temp1 = DIFAB(R(i),temp)
        X(i) = temp1/B(i)
!        X(i) = (R(i)-C(i)*X(i+1))/B(i)
      enddo
!
      end subroutine TRID1
!===============================================================================
!--SUB. TRIX0
      subroutine TRIX0(B,C,D,E,N)
!
!  Solving MX=R, where M is a tridiagonal matrix M=(ABC)
!   plus two corners, upper-right corner: D(1) & lower-left corner: C(N)  
!   where A=1 for all i>1 and D(1)=1. Matrix M is a N*N matrix
!  This is the step to transform the tridiagonal matrix to an
!    upper triangular matrix in two steps, A and B.
!
!  INPUT ARRAYS     : B, C, D
!  WORKING VARIABLES: temp, temp1, temp2
!  OUTPUT ARRAYS    : B(new), C(new), D(new), E
!
      implicit double precision (A-H,O-Z)
      dimension B(1), C(1), D(1), E(1)
!  
!  Step A
!
      D(1) = 1.d0
!
        if (B(1) .eq. 0.d0) then
          print *,'Program TRIX0 stop! B(1) is equal to zero.'
          stop
        endif
!
      do i = 2, N-1
        temp = C(i-1)/B(i-1)
        temp1 = DIFAB(B(i),temp)
!        temp1 = B(i)-temp
!
! IN TRID, B(i), CANNOT BE ZERO  ?
!
        if (temp1 .eq. 0.d0) then
          print *,'Program TRIX0 stop! Denominator is too small.'
          print *,' at I=',i
          print *,' B(I), B(I-1), C(I-1)=',B(i), B(i-1), C(i-1)
          stop
        endif
!
        B(i) = temp1
        D(i) = -D(i-1)/B(i-1)
      enddo
!
      temp = DIFAB(C(N-1),-D(N-1))
      C(N-1) = temp
!      C(N-1) = C(N-1)+D(N-1)
      temp = C(N-1)/B(N-1)
      temp1 = DIFAB(B(N),temp)
      B(N) = temp1
!      B(N) = B(N)-C(N-1)/B(N-1)
!
!  Step B
!
      temp2 = C(N)
!: the last column Nth row = Nth row - E(i)*ith row
      do i = 1, N-2
        E(i) = temp2/B(i)
        temp = E(i)*D(i)
        temp1 = DIFAB(B(N),temp)
        B(N) = temp1
!        B(N) = B(N)-E(i)*D(i)
        temp2 = -E(i)*C(i)
      enddo

      E(N-1) = temp2/B(N-1)
!
      temp1 = DIFAB(B(N),temp)
!      temp1 = B(N)-temp
!
! IN TRID, B(i), CANNOT BE ZERO  ?
!
      if (dabs(temp1) .eq. 0.d0) then
        print *,'Program TRIX0 stop! Denominator is too small.'
        print *,' at I=', N
        print *,' B(I), E(I-1), C(I-1)=',B(N-1), E(N-2), C(N-2)
        stop
      endif
!
      B(N) = temp1
!
      end subroutine TRIX0
!-------------------------------------------------------------------------------
!--SUB. TRIX1
      subroutine TRIX1(B,C,D,E,R,X,N)
!
!  Solving MX=R, where M is a tridiagonal matrix M=(ABC)
!    plus two corners D(1) & C(N-1)  (ps. A=1 for all i, D(1)=1)
!  This is the second step to solve X
!
!  INPUT ARRAYS  : B, C, D, E, R
!  OUTPUT ARRAY  : X
!
      implicit double precision (A-H,O-Z)
      dimension B(1), C(1), D(1), E(1), R(1), X(1)
!
!  When the tridiagonal matrix is transformed to an  
!    upper triangular matrix, 
!    the right hand side column need to change its value.
!
!  Step A
!
      do i = 2, N
        temp = R(i-1)/B(i-1)
        temp1 = DIFAB(R(i),temp)
        R(i) = temp1
!        R(i) = R(i)-R(i-1)/B(i-1)
      enddo
!
!  Step B
!
      do i = 1, N-1
        temp = E(i)*R(i)
        temp1 = DIFAB(R(N),temp)
        R(N) = temp1
!        R(N) = R(N)-E(i)*R(i)
      enddo
!
!  Solve the upper triangular martix to get X
!
      X(N) = R(N)/B(N)
      temp = C(N-1)*X(N)
      temp1 = DIFAB(R(N-1),temp)
      X(N-1) = temp1/B(N-1)
!      X(N-1) = (R(N-1)-C(N-1)*X(N))/B(N-1)

      do i = N-2, 1, -1
        temp = C(i)*X(i+1)
        temp1 = DIFAB(R(i),temp)
        temp = D(i)*X(N)
        temp2 = DIFAB(temp1,temp)
        X(i) = temp2/B(i)
!        X(i) = (R(i)-C(i)*X(i+1)-D(i)*X(N))/B(i)
      enddo

      X(N+1) = X(1)
!
      end subroutine TRIX1
!===============================================================================
end module CubicSpline
