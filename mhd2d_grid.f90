module mhd2d_grid
! Two-dimensional MHD code using FDTD algorithm
! C. L. Waters and M. D. Sciffer
! Centre for Space Physics
! University of Newcastle
! New South Wales 
! Australia
!
! March 2012
! Last reviewed: feb 2013
!
! Contains the following subroutines
! new_r      : determines r for dipole field equation
! do_merge   : used by sort routine
! merge_sort : sort routine
! do_basis   : calc basis functions for atmosphere solution
! gen_grid   : generate solution grid and g factors
!
  use mhd2d_constants

  implicit none

  contains

  subroutine New_r(r0,u1,u3,Th_0r,ans)
! Solve for R using Newton's method
! C.L. Waters
!
    real(DBL), intent(in) :: Th_0r, u1, u3
    real(DBL), intent(inout) :: r0
    real(DBL), intent(out) :: ans

    integer :: MaxN, i
    real(DBL) :: Tol, err, fr, df_dr

    MaxN = 200
    Tol = 1.0d-5
    err = 1.0d0

    do i=1,MaxN
      fr=u3*u3*r0**4.*cos(Th_0r)*cos(Th_0r)/RI_s**3 - u1*r0 - RI_s
      df_dr=4.0*u3*u3*r0**3*cos(Th_0r)*cos(Th_0r)/RI_s**3 - u1
      ans=r0-fr/df_dr
      err=abs(ans-r0)
      r0=ans
      if (err < Tol) exit
    enddo
    return
  end subroutine New_r
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  subroutine do_merge(A,na,B,nb,C,nc,Ai,Bi,Ci)
! routine used by sort
! CLW - March 2012
    implicit none

    integer,intent(in) :: na, nb, nc
    real(DBL), intent(inout) :: A(na), C(nc)
    real(DBL), intent(in) :: B(nb)
    integer, intent(in) :: Bi(nb)
    integer, intent(inout) :: Ai(na), Ci(nc)
    integer :: i, j, k, ii

    i=1; j=1; k=1;
    do while(i <= na .and. j <= nb)
      if (A(i) <= B(j)) then
        C(k) = A(i)
        Ci(k) = Ai(i)
        i=i+1
      else
       C(k) = B(j)
       Ci(k) = Bi(j)
       j=j+1
      end if
      k=k+1
    enddo
    do while (i <= na)
      C(k) = A(i)
      Ci(k) = Ai(i)
      k=k+1
      i=i+1
    enddo

  end subroutine do_merge
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  recursive subroutine merge_sort(A,N,T,Ai,Ti)
! merge sort with recursive call
! C.L. Waters - from wiki page on sort
! Arrays A and Ai are replaced by sorted results
! A : input array of real to be sorted
! N : number of elements
! T : temp array for sort
! Ai : array of indexes of A
! Ti : array of indexes for T
!
    implicit none

    integer, intent(in) :: N
    real(DBL), dimension(N), intent(inout) :: A
    integer, dimension(N), intent(inout) :: Ai
    real(DBL), dimension((N+1)/2), intent(out) :: T
    integer, dimension((N+1)/2), intent(out) :: Ti

    integer :: ii, na, nb, val_i
    real(DBL) :: val

    if (N < 2) return    ! only 1 value so no sort

    if (N==2) then       ! 2 values to sort
      if (A(1) > A(2)) then
        val=A(1)
        val_i=Ai(1)
        A(1)=A(2)
        Ai(1)=Ai(2)
        A(2)=val
        Ai(2)=val_i
      endif
      return
    endif
    na=(N+1)/2
    nb=N-na

    call merge_sort(A,na,T,Ai,Ti)
    call merge_sort(A(na+1),nb,T,Ai(na+1),Ti)
    if (A(na) > A(na+1)) then
      T(1:na)=A(1:na)
      Ti(1:na)=Ai(1:na)
      call do_merge(T,na,A(na+1),nb,A,N,Ti,Ai(na+1),Ai)
    endif
    return
  end subroutine merge_sort
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  subroutine do_basis(m_num, k, Num_u1, L_Min, L_Max, dTh, Ynm, Ynm_s, zeros)
! Construct basis set of pV/pr for atmosphere solution
! M. D. Sciffer and C. L. Waters
! Centre for Space Physics
! University of Newcastle
! Australia
!
! March 2012
!

    use mhd2d_constants

! Input variables
    integer, intent(in) :: k, Num_u1
    real(DBL), intent(in) :: m_num, L_Min, L_Max
    real(DBL), intent(inout) :: dTh
! Outputs:
    real(DBL), dimension(k,Num_u1), intent(out) :: Ynm, Ynm_s
    real(DBL), dimension(k), intent(out) :: zeros
! Vars for rest of code
    integer :: ii, jj, ipts, LowBC, HighBC
    real(DBL) :: dTh2, dTh_sq, End1, End2, res, p0, Pn
    real(DBL) :: CLat_Min, CLat_Max, CosT, SinT

    real(DBL), dimension(Num_u1) :: x, Theta
    real(DBL), dimension(K,K) :: ortho_n
    real(DBL), dimension(K) :: RTerm_dr, RTerm

    integer, allocatable :: Iperm(:)
    real(DBL), allocatable :: Array(:,:)
    real(DBL), allocatable :: Lp(:)
! variable for sort routine
    real(DBL), dimension(:), allocatable :: T
    integer, dimension(:), allocatable :: Ti
! Matrix eigen_solver variables
    character :: jobVL='N', jobVR='V'
    integer :: LWork, LDVL, info
    real(DBL), allocatable :: VL(:,:), VR(:,:), &
                              wR(:), wI(:), Work(:)

    ipts = Num_u1-2
    LDVL = 1
    LWORK = 4*ipts
! Index array for sort routine
    allocate(Iperm(ipts))
    do ii=1,ipts
      Iperm(ii)=ii
    enddo

    CLat_Min=asin(sqrt(rI_s/L_Max))*RadToDeg       ! Min Co_Lat in degrees, specifies outer boundary
    CLat_Max=asin(sqrt(rI_s/L_Min))*RadToDeg       ! Max Co_Lat in degrees for Nth Hemisphere (Min Latitude), specifies inner boundary

    dTh=(CLat_Max-CLat_Min)/(Num_u1-1.d0)     ! del_Theta in degrees

! Low here is low theta or poleward end
! High here is high theta or equatorward point
! If Thi=0 then B2(Azimuthal) is zero else B1 (theta) is zero
    LowBC = 0          ! '1' -> Thi = 0 '0' -> Dthi/Dtheta = 0, PoleWard condition : If Thi=0 then this also sets B2=0
    HighBC = 0          ! '1' -> Thi = 0  '0' -> Dthi/Dtheta = 0, Equatorward condition has B1=0 in the inner boundary condition

! Get cos(theta) points for Legendres
    do ii = 1,Num_u1
      Theta(ii) = (CLat_Min+((ii-1.d0)*dth))*DegToRad    ! Northern Hemisphere
      X(ii) = cos(Theta(ii))
    enddo
    dTh = dth*DegToRad     ! delta_theta in radians
    dTh2 = 2.d0*dTh
    dTh_sq = dTh*dTh       ! del_theta squared

! Numerical differentiation approx of Laplace Eqn in Spherical coords
! Interior Points
    allocate ( Array(ipts,ipts) )
    do ii = 2, ipts-1
      SinT  = sin(Theta(ii+1))
      CosT  = cos(Theta(ii+1))
      Array(ii,ii-1)=  1.d0/dTh_sq - CosT/SinT/dTh2
      Array(ii,ii)  = -2.d0/dTh_sq - (M_num**2)/(SinT**2)
      Array(ii,ii+1)=  1.d0/dTh_sq + CosT/SinT/dTh2
    enddo

! Low Latitude End - Forward differences ( dp/dT = 0  -> P(0) =  4/3 P(1) - 1/3 P(2) )
    SinT= Sin(Theta(2))
    CosT= Cos(Theta(2))
    p0 = 0.d0
    If (LowBC == 0) p0=1.d0/dTh_sq - CosT/SinT/dTh2
    Array(1,1) = -2.d0/dTh_sq - (M_num**2)/(SinT**2) + 4.d0/3.d0*p0
    Array(1,2) =  1.d0/dTh_sq + CosT/SinT/dTh2 - 1.d0/3.d0*p0

! High Latitude End - Backwards differences ( dp/dT = 0  -> P(n) =  4/3 P(n-1) - 1/3 P(n-2) )
    SinT= sin(Theta(Num_u1-1))
    CosT= cos(Theta(Num_u1-1))
    p0 = 0.d0
    If (HighBC == 0) pn=1.d0/dTh_sq + CosT/SinT/dTh2              ! Derivative = 0
    Array(ipts,ipts)  = -2.d0/dTh_sq - (M_num**2)/(sinT**2) + 4.d0/3.d0*pn
    Array(ipts,ipts-1)=  1.d0/dTh_sq - cosT/sinT/dTh2 - 1.d0/3.d0*pn

! Solve Ax=lmbda x Eqn for lmbda and x, lbmda=-l(l+1)
    allocate (VR(ipts,ipts))
    allocate (wR(ipts))
    allocate (WORK(4*ipts))
    allocate (VL(1,ipts))
    allocate (wI(ipts))

    call dgeev(jobVL, jobVR, ipts, Array, ipts, wR, wI, VL, LDVL, &
               VR, ipts, WORK, LWORK, info)

    deallocate (wI)
    deallocate (VL)
    deallocate (WORK)

! Calculate eigenvalues
    allocate (Lp(ipts))
    Lp = 0.d0            ! Initialise array to 0.0
    Do ii = 1, ipts
      Lp(ii) = (-1.d0+sqrt(1.d0-4.d0*(wR(ii))))/2.d0  ! Solve quadratic lmbda=-l**2-l
    enddo
    deallocate (wR)

! sort eigenvalues - Lp
    allocate (T((ipts+1)/2))
    allocate (Ti((ipts+1)/2))
    do ii=1,ipts
      Iperm(ii)=ii   ! get idx
    enddo
! trap for constant term (exclude it)
    where (Lp <= 0.1) Lp = 9999.0
    call merge_sort(Lp, ipts, T, Iperm, Ti)
    deallocate(Ti)
    deallocate(T)
! eigenvalues are now sorted, get K of them
    do ii=1,K
      Zeros(ii) = Lp(ii)    ! Store sorted solns to eigenvalue quadratic
    enddo

! Now get the sorted eigenvectors - these form the potential basis function set
    do ii = 1,K 
      do jj = 2, Num_u1-1 
        Ynm(ii,jj) = VR(jj-1,Iperm(ii))
      enddo
    enddo
    deallocate (Lp)
    deallocate (VR)
    deallocate (Array)

! Reconstruct end points of eigenvectors
    If (LowBC == 0) then
      do ii = 1,K 
        Ynm(ii,1)=4.d0/3.d0*Ynm(ii,2)-1.d0/3.d0*Ynm(ii,3)   ! Low  Lat Boundary
      enddo
    else 
      Ynm(:,1) = 0.d0                ! Low  Lat Boundary, Thi=0
    endif

    If (HighBC == 0) then
      do ii = 1,K                    ! Derivative = 0
        Ynm(ii,Num_u1)=4.d0/3.d0*Ynm(ii,Num_u1-1)-1.d0/3.d0*Ynm(ii,Num_u1-2)     ! High Lat Boundary
      enddo
    else 
      Ynm(:,Num_u1) = 0.d0           ! High Lat Boundary
    endif

!   Normalise Basis Functions Numerically (using trapezoidal rule)
    do ii = 1,K
      do jj = 1 ,K
        End1 =(Ynm(ii,1)*sqrt(sin(theta(1))))*(Ynm(jj,1)*sqrt(sin(theta(1))))
        End2 =(Ynm(ii,Num_u1)*sqrt(sin(theta(Num_u1))))*(Ynm(jj,Num_u1)*sqrt(sin(theta(Num_u1))))
        res=(2.0*Dot_Product(Ynm(ii,:)*sqrt(sin(theta(:))),Ynm(jj,:)*sqrt(sin(theta(:)))) - End1 - End2)*(dTh/2.0)
        ortho_n(ii,jj) = res
      enddo
    enddo

    do ii = 1,K
      Ynm(ii,:) = Ynm(ii,:)/sqrt(ortho_n(ii,ii))      ! Basis functions are now an orthonormal basis set
    enddo

    do ii = 1,K
      RTerm_dr(ii) = zeros(ii)*RI_s**(zeros(ii)-1.0)*(1.0-(rE_s/RI_s)**(2.0*zeros(ii)+1.0))          ! Radial Derivative
      RTerm(ii) = RI_s**zeros(ii) + (zeros(ii)/(zeros(ii)+1.0))*rE_s**(2.0*zeros(ii)+1.0)* &
                  RI_s**(-(zeros(ii)+1.0))            ! Radial Term
      Ynm_S(ii,:) = Ynm(ii,:)*RTerm_dr(ii)
    enddo

  end subroutine do_basis
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  subroutine gen_grid ( Num_u1, Num_u1_2, Num_u3, Num_u3_2, &
                       LMin_u1, LMax_u1, LMin_u3, LMax_u3, &
                       del_Th_0_u1, del_Th, &
             u1Arr, u3Arr, RArr, ThArr, LVArr, XArr, YArr, &
             g11_con, g13_con, g22_con, g33_con, &
             g11_cov, g13_cov, g22_cov, g33_cov, &
             Jacb, Colat_N, Colat_S, h_nu, h_ph, h_mu)

! Set up spatial grid and scale factors for the Non-Orthogonal Coords
! C.L. Waters
! June, 2004
!
    implicit none

    integer, intent(in) :: Num_u1, Num_u1_2, Num_u3, Num_u3_2
    real(DBL), intent(in) :: LMin_u1, LMax_u1, LMin_u3, LMax_u3
    real(DBL), intent(inout) :: del_Th_0_u1, del_Th

    real(DBL), intent(inout), dimension(Num_u1) :: LVArr, CoLat_N, CoLat_S
    real(DBL), intent(inout), dimension(Num_u1,Num_u3) :: u1Arr, u3Arr, RArr, ThArr
    real(DBL), intent(inout), dimension(Num_u1,Num_u3) :: xArr, yArr
    real(DBL), intent(inout), dimension(Num_u1,Num_u3) :: g11_con, g13_con, g22_con, g33_con
    real(DBL), intent(inout), dimension(Num_u1,Num_u3) :: g11_cov, g13_cov, g22_cov, g33_cov
    real(DBL), intent(inout), dimension(Num_u1,Num_u3) :: Jacb
    real(DBL), intent(inout), dimension(Num_u1,Num_u3) :: h_nu, h_ph, h_mu

    integer :: ii, jj
    real(DBL) :: CLat_Min_u1, CLat_Max_u1, CLat_Min_u3, CLat_Max_u3
    real(DBL) :: del_Th_0_u3
    real(DBL) :: CLat_0_u1, CLat_0r_u1, CLat_0_u3, CLat_0r_u3
    real(DBL) :: u1, u3, rg, ans, cos_th, sin_th, theta, lat, one_p_3cos_sq

! Determine del_Theta for u1 calcs...
    CLat_Min_u1=asin(sqrt(RI_S/LMax_u1))*RadToDeg     ! Min Co_Lat in degrees for Nth Hemisphere (Max Latitude), specifies outer boundary
    CLat_Max_u1=asin(sqrt(RI_s/LMin_u1))*RadToDeg     ! Max Co_Lat in degrees for Nth Hemisphere (Min Latitude), specifies inner boundary
    del_Th_0_u1=(CLat_Max_u1-CLat_Min_u1)/(Num_u1-1.d0)   ! del_Theta for u1 points
    del_Th = del_Th_0_u1*DegToRad                    ! del_Thet in radians

    Print*,'Non-Orthogonal Grid Generation:'
    Print*,'L_Value_Min, L_Value_Max, CLat_Min_u1, CLat_Max_u1 (Deg) : '
    Print*, LMin_u1, LMax_u1, CLat_Min_u1, CLat_Max_u1
    Print*,'Number of u1 Lines, Del_Theta_u1 : ',Num_u1, del_Th_0_u1

    CLat_Min_u3=asin(sqrt(RI_s/LMax_u3))*RadToDeg
    CLat_Max_u3=asin(sqrt(RI_s/LMin_u3))*RadToDeg
    del_Th_0_u3=(CLat_Max_u3-CLat_Min_u3)/(Num_u3_2-1.d0)
    Print*,'Number of u3 Lines, Del_Theta_u3 : ',Num_u3, del_Th_0_u3

!  u1 (L shell) Loop, Northern Hemisphere
    do jj=1,Num_u1                                 ! Starts at inner-most field line (i.e. smallest L)
      CLat_0_u1=CLat_Max_u1-(jj-1.)*del_Th_0_u1    ! Co_Lat at Nth hemisphere Ionosphere (in degrees)
      CLat_0r_u1=CLat_0_u1*DegToRad                ! Convert CLat_0 to radians
      u1=-sin(CLat_0r_u1)*sin(CLat_0r_u1)          ! Determine u1 (which is constant along an L shell) for this L value
      Colat_N(jj) = CLat_0r_u1
      u1Arr(jj,:)=u1
      LVArr(jj)=1.0/(cos(pi/2.0-CLat_0r_u1)**2)
      do ii=1, Num_u3_2
        CLat_0_u3=CLat_Max_u3-(CLat_Max_u3-CLat_Min_u3)* ((ii-1.)/(Num_u3_2-1.d0))**0.5    ! Co_Lat at Nth hemisphere Ionosphere (in degrees)
        CLat_0r_u3=CLat_0_u3*DegToRad                ! Convert CLat_0 to radians
        u3=sin(CLat_0r_u3)**4                        ! Determine u3 for Nth (determined at CLat=0)
        u3Arr(jj,ii)=u3
        rg=RI_s
        call New_r(rg,u1,u3,CLat_0r_u1,ans)
        RArr(jj,ii)=ans                              ! r value for the intersection
        cos_Th=u3*cos(CLat_0r_u1)*ans*ans/(RI_s*RI_s)! see eqn (11) Lysak, 2004
        one_p_3cos_sq=1.0+3.0*cos_Th*cos_Th
        Theta=dacos(cos_Th)
        sin_Th=sin(Theta)
        ThArr(jj,ii)=Theta
        Lat=pi/2.0-Theta
        xArr(jj,ii)=ans*cos(Lat)
        yArr(jj,ii)=ans*sin(Lat)
        g11_con(jj,ii)= (RI_s**2/ans**4)*sin_Th**2*one_p_3cos_sq
        g13_con(jj,ii)=-(RI_s**3/ans**5)*sin(CLat_0r_u1)**2*cos_Th*one_p_3cos_sq/(2.0*cos(CLat_0r_u1)**3)
        g22_con(jj,ii)=1.0/(ans**2*sin_Th**2)
        g33_con(jj,ii)=RI_s**4/(ans**6*cos(CLat_0r_u1)**6)* &
          (0.25*cos_Th**2*(1.0+3.0*cos(CLat_0r_u1)**2)**2+sin_Th**2*(1.0-RI_s/ans)**2)
        g11_cov(jj,ii)=ans**4/(RI_s**2*cos(CLat_0r_u1)**4*one_p_3cos_sq**2)* &
          ((1.0-RI_s/ans)**2+0.25*(cos_Th/sin_Th)**2*(1.0+3.0*cos(CLat_0r_u1)**2)**2)
        g13_cov(jj,ii)=ans**4*cos_Th/(2.0*RI_s**2*cos(CLat_0r_u1)*one_p_3cos_sq)
        g22_cov(jj,ii)=ans*ans*sin_Th*sin_Th
        g33_cov(jj,ii)=ans**6*cos(CLat_0r_u1)**2/(RI_s**4*one_p_3cos_sq)
        Jacb(jj,ii)=ans**6*cos(CLat_0r_u1)/(RI_s**3*one_p_3cos_sq)

        h_nu(jj,ii)= 1.0d0/sqrt(g11_con(jj,ii))   ! At the Equator this is in the R direction
        h_ph(jj,ii)= 1.0d0/sqrt(g22_con(jj,ii))   ! E_W, Azimuthal
        h_mu(jj,ii)= sqrt(g33_cov(jj,ii))         ! At the Equator this is in the Theta direction (FA)  enddo

! Southern Hemisphere
! Mirror hemisphere by changing Theta
        u3Arr(jj,Num_u3-ii+1)=-u3                 ! Keep same u3 value
        RArr(jj,Num_u3-ii+1)=ans                  ! Keep same r value
        Theta=pi-ThArr(jj,ii)                     ! Change Theta here for Sth hemisphere
        ThArr(jj,Num_u3-ii+1)=Theta
        Lat=(pi/2.0-Theta)
        xArr(jj,Num_u3-ii+1)=ans*cos(Lat)
        yArr(jj,Num_u3-ii+1)=ans*sin(Lat)
        cos_Th=cos(Theta)
        one_p_3cos_sq=1.0+3.0*cos_Th*cos_Th
        sin_Th=sin(Theta)
        g11_con(jj,Num_u3-ii+1)= (RI_s**2/ans**4)*sin_Th**2*one_p_3cos_sq
        g13_con(jj,Num_u3-ii+1)=-(RI_s**3/ans**5)*sin(CLat_0r_u1)**2* &
          cos_Th*one_p_3cos_sq/(2.0*cos(CLat_0r_u1)**3)
        g22_con(jj,Num_u3-ii+1)=1.0/(ans**2*sin_Th**2)
        g33_con(jj,Num_u3-ii+1)=RI_s**4/(ans**6*cos(CLat_0r_u1)**6)* &
          (0.25*cos_Th**2*(1.0+3.0*cos(CLat_0r_u1)**2)**2 + sin_Th**2*(1.0-RI_s/ans)**2)
        g11_cov(jj,Num_u3-ii+1)=ans**4/(RI_s**2*cos(CLat_0r_u1)**4*one_p_3cos_sq**2)* &
          ((1.0-RI_s/ans)**2 + 0.25*(cos_Th/sin_Th)**2*(1.0+3.0*cos(CLat_0r_u1)**2)**2)
        g13_cov(jj,Num_u3-ii+1)=ans**4*cos_Th/(2.0*RI_s**2*cos(CLat_0r_u1)*one_p_3cos_sq)
        g22_cov(jj,Num_u3-ii+1)=ans*ans*sin_Th*sin_Th
        g33_cov(jj,Num_u3-ii+1)=ans**6*cos(CLat_0r_u1)**2/(RI_s**4*one_p_3cos_sq)
        Jacb(jj,Num_u3-ii+1)=ans**6*cos(CLat_0r_u1)/(RI_s**3*one_p_3cos_sq)

        h_nu(jj,Num_u3-ii+1)= 1.0d0/sqrt(g11_con(jj,Num_u3-ii+1))   ! At the Equator this is in the R direction
        h_ph(jj,Num_u3-ii+1)= 1.0d0/sqrt(g22_con(jj,Num_u3-ii+1))   ! E_W, Azimuthal
        h_mu(jj,Num_u3-ii+1)= sqrt(g33_cov(jj,Num_u3-ii+1))         ! At the Equator this is in the Theta direction (FA)
      enddo
    enddo
!
!  Do Equator, located at Num_u3_2 position
    do jj=1,Num_u1
      CLat_0_u1=CLat_Max_u1-(jj-1.0)*del_Th_0_u1            ! Co_Lat at Nth hemisphere Ionosphere (in degrees)
      CLat_0r_u1=CLat_0_u1*DegToRad                         ! Convert CLat_0 to radians
      ThArr(jj,Num_u3_2+1)=pi/2.0
      xArr(jj,Num_u3_2+1)=RI_s/sin(CLat_0r_u1)**2
      yArr(jj,Num_u3_2+1)=0.0d0
      RArr(jj,Num_u3_2+1)=RI_s/sin(CLat_0r_u1)**2
      u1Arr(jj,Num_u3_2+1)=-sin(CLat_0r_u1)**2
! u3 should all be zero along the equator
      Theta=pi/2.0                      ! Change Theta here for Equator
      ThArr(jj,Num_u3_2+1)=Theta
      Lat = 0.0d0
      cos_Th=cos(Theta)
      one_p_3cos_sq=1.0+3.0*cos_Th*cos_Th
      sin_Th=sin(Theta)
      ans = Rarr(jj,Num_u3_2+1)
      g11_con(jj,Num_u3_2+1)=(RI_s**2/ans**4)*sin_Th**2*one_p_3cos_sq
      g13_con(jj,Num_u3_2+1)=-(RI_S**3/ans**5)*sin(CLat_0r_u1)**2*cos_Th*one_p_3cos_sq/(2.0*cos(CLat_0r_u1)**3)
      g22_con(jj,Num_u3_2+1)=1.0/(ans**2*sin_Th**2)
      g33_con(jj,Num_u3_2+1)=RI_S**4/(ans**6*cos(CLat_0r_u1)**6)* &
        (0.25*cos_Th**2*(1.0+3.0*cos(CLat_0r_u1)**2)**2 + sin_Th**2*(1.0-RI_s/ans)**2)
      g11_cov(jj,Num_u3_2+1)=ans**4/(RI_s**2*cos(CLat_0r_u1)**4*one_p_3cos_sq**2)* &
        ((1.0-RI_s/ans)**2 + 0.25*(cos_Th/sin_Th)**2*(1.0+3.0*cos(CLat_0r_u1)**2)**2)
      g13_cov(jj,Num_u3_2+1)=ans**4*cos_Th/(2.0*RI_S**2*cos(CLat_0r_u1)*one_p_3cos_sq)
      g22_cov(jj,Num_u3_2+1)=ans*ans*sin_Th*sin_Th
      g33_cov(jj,Num_u3_2+1)=ans**6*cos(CLat_0r_u1)**2/(RI_s**4*one_p_3cos_sq)
      Jacb(jj,Num_u3_2+1)=ans**6*cos(CLat_0r_u1)/(RI_s**3*one_p_3cos_sq)

      h_nu(jj,Num_u3_2+1)= 1.0d0/sqrt(g11_con(jj,Num_u3_2))   ! At the Equator this is in the R direction
      h_ph(jj,Num_u3_2+1)= 1.0d0/sqrt(g22_con(jj,Num_u3_2))   ! E_W, Azimuthal
      h_mu(jj,Num_u3_2+1)= sqrt(g33_cov(jj,Num_u3_2))         ! At the Equator this is in the Theta direction (FA)
    enddo

  end Subroutine gen_grid

end module mhd2d_grid
!
! ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
