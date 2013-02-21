!   Simulation of ULF wave propagation through the magnetosphere
!   with Realistic Ionospheric Boundaries
!   Version 2.0     5 Aug 2006
!
!   C.L. Waters, M.D. Sciffer
!   Centre for Space Physics
!   University of Newcastle
!   New South Wales, Australia
!
! Mods:
!
!  In order to restart, output the whole data set every FULL_FREQ steps.
!  Some parts/subroutines are deleted if they're not used here
!  See the original IDL code if necessary.
!
! ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
!
  program mhd2d

    use mhd2d_constants
    use mhd2d_va_cond
    use mhd2d_grid

!    include 'dgeevx.f'

    implicit none

    integer ::  ii, jj
    integer ::  Num_u1, Num_u1_2, Num_u3, Num_u3_2, K
    integer(8) :: nt

    real(DBL) :: LMin_u1, LMax_u1, LMin_u3, LMax_u3, m_num
    real(DBL) :: dt, del_Th_0_u1, del_Th

! Va vars
    integer :: i1, i2
    real(DBL) :: xd, va_sp, lat, bmg, PLaw, nH, nO, nO_Ht
    real(DBL), allocatable, dimension(:,:) :: V2, Va_a

! grid arrays
    real(DBL), allocatable, dimension(:) :: LVArr, CoLat_N, CoLat_S, h_mu_outer
    real(DBL), allocatable, dimension(:,:) :: u1Arr, u3Arr, RArr, ThArr
    real(DBL), allocatable, dimension(:,:) :: XArr, YArr
    real(DBL), allocatable, dimension(:,:) :: g11_con, g13_con, g22_con, g33_con
    real(DBL), allocatable, dimension(:,:) :: g11_cov, g13_cov, g22_cov, g33_cov
    real(DBL), allocatable, dimension(:,:) :: Jacb, rho_a
    real(DBL), allocatable, dimension(:,:) :: h_nu, h_ph, h_mu
    complex(8), allocatable, dimension(:) :: odds, evens

! field ord arrays
    real(DBL), allocatable, dimension(:,:) :: B1x_ord, B2x_ord, B3x_ord, &
                                              B1y_ord, B2y_ord, B3y_ord, &
                                              E1x_ord, E2x_ord, &
                                              E1y_ord, E2y_ord
! Ionosphere conductance arrays and vars
    real(DBL), allocatable, dimension(:) :: Sd_N,Sp_N,Sh_N, &
                                            Sd_S,Sp_S,Sh_S, &
                                            h_ph_n,h_th_n,h_ra_n,&
                                            h_ph_s,h_th_s,h_ra_s
    real(DBL), allocatable, dimension(:,:,:) :: INVSS_N,INVSS_S
    character(LEN=80) :: iriva_file
    real(SGL) ::  yr, doy, ut_hr, geolon,f107,f107a,ap
    real(DBL) :: s0n,s1n,s2n,s0s,s1s,s2s

! Courant condition arrays
    real(DBL), allocatable, dimension(:,:) :: Delta_T_CF, Delta_T_FA, &
                                              Delta_L_CF, Delta_L_FA
    real(DBL), allocatable, dimension(:) :: D_out, Out_DL, &
                                            FLR_Per
    real(DBL) :: Av_Va, Delta_x, Delta_y, Bnd_L, Length, Courant_cond, &
                 Freq, Eperp

! Basis set arrays
    real(DBL), allocatable, dimension(:,:) :: Ynm, Ynm_s
    real(DBL), allocatable, dimension(:) :: zeros


! - - - - - - - - - start of code - - - - - - - - - - -

! Parameters for run - to come from file
    iriva_file='/home/mdw/mhd2d/iri_va_sig.mhd'

    open(unit=10,file=iriva_file,status='old',action='read')
    read(10,*) yr,doy,ut_hr,geolon,Num_u1,Num_u3
    read(10,*) f107,f107a,ap
    read(10,*) LMin_u1, LMax_u1, LMin_u3,LMax_u3

!    Num_u1 = 210
!    Num_u3 = 151
!    LMin_u1 = 1.2d0
!    LMax_u1 = 15.0d0
!    LMin_u3 = RI_s
!    LMax_u3 = 95.0d0

! Number of iterations
    nt = 5000
! dt time step
    dt = 0.00025

    m_num = 2.0

! For additional density near ionosphere
    PLaw = 4.0
    nO = 1.0d7             ! O2 density (kg/Re^3)
    nO_Ht = 250.0d3/rE_m   ! scale height of oxygen
  
! - - - - - end parameters - - - - -

! Check grid size and numbers
    Num_u1_2 = int(Num_u1/2.0)
    Num_u1=int(2*Num_u1_2)

    Num_u3_2=int(Num_u3/2.0)
    Num_u3=int(2*Num_u3_2)+1

! spherical harmonic soln order for atmosphere
    K = int(Num_u1/2 - 2)

! setup solution grid arrays
    allocate (LVArr(Num_u1), CoLat_N(Num_u1), CoLat_S(Num_u1) )
    allocate (u1Arr(Num_u1,Num_u3), u3Arr(Num_u1,Num_u3) )
    allocate (RArr(Num_u1,Num_u3), ThArr(Num_u1,Num_u3) )
    allocate (xArr(Num_u1,Num_u3), yArr(Num_u1,Num_u3) )

    allocate (g11_con(Num_u1,Num_u3), g13_con(Num_u1,Num_u3) )
    allocate (g22_con(Num_u1,Num_u3), g33_con(Num_u1,Num_u3) )

    allocate (g11_cov(Num_u1,Num_u3), g13_cov(Num_u1,Num_u3) )
    allocate (g22_cov(Num_u1,Num_u3), g33_cov(Num_u1,Num_u3) )

    allocate (Jacb(Num_u1,Num_u3) )

    allocate (h_nu(Num_u1,Num_u3) )
    allocate (h_ph(Num_u1,Num_u3) )
    allocate (h_mu(Num_u1,Num_u3) )

    allocate (rho_a(Num_U1,Num_u3) )

! Call grid generation routine
    print*,'Generating grid...'
    call gen_grid(Num_u1, Num_u1_2, Num_u3, Num_u3_2, &
           LMin_u1, LMax_u1, LMin_u3, LMax_u3, &
           del_Th_0_u1, del_Th, &
           u1Arr, u3Arr, RArr, ThArr, LVArr, xArr, yArr, &
           g11_con, g13_con, g22_con, g33_con, &
           g11_cov, g13_cov, g22_cov, g33_cov, &
           Jacb, CoLat_N, CoLat_S, h_nu, h_ph, h_mu)
! # # #
! Test output for grid to fort.1
!    do ii=1,Num_u1
!      do jj=1,Num_u3
!        write(1,*) ii,jj, u1Arr(ii,jj), u3Arr(ii,jj)
!        write(2,*) ii,jj, g11_con(ii,jj), g11_cov(ii,jj)
!        write(3,*) ii,jj, g13_con(ii,jj), g13_cov(ii,jj)
!        write(4,*) ii,jj, g22_con(ii,jj), g22_cov(ii,jj)
!        write(7,*) ii,jj, g33_con(ii,jj), g33_cov(ii,jj)
!        write(8,*) ii,jj, Jacb(ii,jj)
!        write(9,*) ii,jj, h_nu(ii,jj),h_ph(ii,jj),h_mu(ii,jj)
!      enddo
!    enddo
!    stop
! # # #

! Read in Va data
    allocate(Va_a(Num_u1,Num_u3))    
    do ii=1,Num_u1
      do jj=1,Num_u3
        read(10,*) i1, i2, lat, va_sp
        Va_a(i1+1,i2+1) = va_sp
      enddo
    enddo

    allocate (Sd_N(Num_u1))
    allocate (Sp_N(Num_u1))
    allocate (Sh_N(Num_u1))
    allocate (Sd_S(Num_u1))
    allocate (Sp_S(Num_u1))
    allocate (Sh_S(Num_u1))

! Read in conductance data
    do ii=1,Num_u1
      read(10,*) lat, s0n,s1n,s2n,s0s,s1s,s2s
      Sd_N(ii)=s0n
      Sp_N(ii)=s1n
      Sh_N(ii)=s2n
      Sd_S(ii)=s0s
      Sp_S(ii)=s1s
      Sh_S(ii)=s2s
    enddo
    close(unit=10)

    print*,'Calculating VSa...'
    do ii=1,Num_u1                   ! chooses which field line
      xd = LVArr(ii)
      lat = 0.0d0                    ! set to equator

!      call clw_va(xd,va_sp)         ! Given xd, calc Va at equator
      va_sp=Va_a(ii,Num_u3_2)        ! use data from input file

      call Get_B0(Lat,xd,bmg)        ! Given xd and Latitude, calc B_0 at equator
      nH   = bmg**2/(mu0*va_sp**2)   ! Calc Hydrogen plasma number density at the equator

      do jj=1,Num_u3                 ! Loop over points along field line
        Lat=pi/2.0-ThArr(ii,jj)
        call Get_B0(Lat,xd,bmg)
!   rho_a(ii,jj)=n0*(xd/RArr(ii,jj))**3

        rho_a(ii,jj) = nH*(xd/RArr(ii,jj))**PLaw + nO/rE_m**3* &
          exp(-(Rarr(ii,jj)-Re_s)/nO_Ht)

        Va_a(ii,jj) = bmg/sqrt(mu0*rho_a(ii,jj))/rE_m  ! modify Va
      enddo                          ! Num_u3 loop

! call eigen_solver to check FLR frequencies
!   call Find_flr(xd,N3,rho_a(ii,:),har)
!   flr_a(ii,0:5)=har(0:5)           ! Fundamental FLR
    enddo                            ! Num_u1 loop

! Get memory for scale factors and conductance matrices
    allocate (h_ph_n(Num_u1))
    allocate (h_th_n(Num_u1))
    allocate (h_ra_n(Num_u1))
    allocate (h_ph_s(Num_u1))
    allocate (h_th_s(Num_u1))
    allocate (h_ra_s(Num_u1))

    allocate (INVSS_N(Num_u1,2,2))
    allocate (INVSS_S(Num_u1,2,2))
!
    print*,'Calculating ionosphere conductance...'
    call Calc_cond_ionos(Num_u1, CoLat_N, CoLat_S, &
                 Sd_N,Sp_N,Sh_N,INVSS_N, &
                 Sd_S,Sp_S,Sh_S,INVSS_S, &
                 h_ph_n,h_th_n,h_ra_n, &
                 h_ph_s,h_th_s,h_ra_s)
! stop

! Discretisation Scheme
!
! The Number of full cells is (N1_2,N3_2)
! In this programing section ii and jj indentifies the field values in the (ii,jj) cell
! and within the cell the location of parameters (such as alfven speed etc ) associated
! with a particular position in space are given by ....
!
!       (2*ii-1, 2*jj-1)      location of E2
!       (2*ii,   2*jj-1)      location of E1
!       (2*ii-1, 2*jj)        location of B1
!       (2*ii,   2*jj)        location of B2
!       (2*ii,   2*jj-1)      location of B3   ( located at the same position in space as E1)
!
!   Setting up Array for plots of fields
!
  allocate (B1x_ord(Num_u1_2,Num_u3_2))
  allocate (B2x_ord(Num_u1_2,Num_u3_2))
  allocate (B1y_ord(Num_u1_2,Num_u3_2))
  allocate (B2y_ord(Num_u1_2,Num_u3_2))

  allocate (B3x_ord(Num_u1_2,Num_u3_2+1))
  allocate (B3y_ord(Num_u1_2,Num_u3_2+1))
  allocate (E1x_ord(Num_u1_2,Num_u3_2+1))
  allocate (E1y_ord(Num_u1_2,Num_u3_2+1))
  allocate (E2x_ord(Num_u1_2,Num_u3_2+1))
  allocate (E2y_ord(Num_u1_2,Num_u3_2+1))

  allocate (odds(Num_u1))
  allocate (evens(Num_u1))
  odds = 0.0
  evens = 0.0
  
  do ii = 1,Num_u1_2
    evens(2*ii) = 1.0
    odds(2*ii-1) = 1.0
    do jj = 1,Num_u3_2                       ! For component on the interior
      B1x_ord(ii,jj) = Xarr(2*ii-1 ,2*jj)    ! X Cordinate for Plot routines of various fields
      B2x_ord(ii,jj) = Xarr(2*ii   ,2*jj)
      B3x_ord(ii,jj) = Xarr(2*ii   ,2*jj-1)
      E1x_ord(ii,jj) = Xarr(2*ii   ,2*jj-1)
      E2x_ord(ii,jj) = Xarr(2*ii-1 ,2*jj-1)

      B1y_ord(ii,jj) = Yarr(2*ii-1 ,2*jj)    ! Y Cordinate for Plot routines of various fields
      B2y_ord(ii,jj) = Yarr(2*ii   ,2*jj)
      B3y_ord(ii,jj) = Yarr(2*ii   ,2*jj-1)
      E1y_ord(ii,jj) = Yarr(2*ii   ,2*jj-1)
      E2y_ord(ii,jj) = Yarr(2*ii-1 ,2*jj-1)
!      write(11,'(2(i4,1X),10(f10.3,1X))') ii,jj, B1x_ord(ii,jj), B2x_ord(ii,jj), B3x_ord(ii,jj), E1x_ord(ii,jj), E2x_ord(ii,jj), B1y_ord(ii,jj), B2y_ord(ii,jj), B3y_ord(ii,jj), E1y_ord(ii,jj), E2y_ord(ii,jj)
    enddo
  enddo 
! CHECKED
!   Addition Half cell in u1 direction
  do ii = 1,Num_u1_2
    E1x_ord(ii,Num_u3_2) = Xarr(2*ii  ,Num_u3)
    E2x_ord(ii,Num_u3_2) = Xarr(2*ii-1,Num_u3)
    B3x_ord(ii,Num_u3_2) = Xarr(2*ii  ,Num_u3)

    E1y_ord(ii,Num_u3_2) = Yarr(2*ii  ,Num_u3)
    E2y_ord(ii,Num_u3_2) = Yarr(2*ii-1,Num_u3)
    B3y_ord(ii,Num_u3_2) = Yarr(2*ii  ,Num_u3)

!    write(11,'(i4,1X,6(f10.3,1X))') ii, E1x_ord(ii,Num_u3_2), E2x_ord(ii,Num_u3_2), B3x_ord(ii,Num_u3_2), E1y_ord(ii,Num_u3_2), E2y_ord(ii,Num_u3_2), B3y_ord(ii,Num_u3_2)
  enddo 

! Calc Va^2 used in time iteration which is time invariant over the full grid
  allocate (V2(Num_u1, Num_u3))
! CHECKED
  do ii=1,Num_u1
    do jj=1,Num_u3
      Eperp  = eps0_s*(1.0 + c_s**2/Va_a(ii,jj)**2)   ! in SI units
      V2(ii,jj) = 1.0/(mu0_s*Eperp)  ! Va^2 in SI
!      write(11,'(2(i4,1X),f10.6,1X,2E14.6)') ii, jj, V2(ii,jj), Eperp, Va_a(ii, jj)
!      Sa(ii,jj)       = 1.0/(u0*Va_a(ii,jj))
    enddo
  enddo

!
! Calculate Info about run
!
  Bnd_L   = 0.0
  allocate (Delta_L_FA(Num_u1,Num_u3))
  allocate (Delta_T_FA(Num_u1,Num_u3))
!  allocate (F_Len(Num_u1))
  allocate(FLR_Per(Num_u1))

  do ii = 1, Num_u1
    do jj = 2, Num_u3-1
      Delta_x = (xArr(ii,jj+1) - xArr(ii,jj-1))           !     Change in X coordinate along field line
      Delta_y = (yArr(ii,jj+1) - yArr(ii,jj-1))           !     Change in Y coordinate along field line
      Delta_L_FA(ii,jj)= sqrt(Delta_x**2 + Delta_y**2)    !     Line segment Length (in Re)
!      F_len(ii) = F_len(ii) + Delta_L_FA(ii,jj)
      Av_Va =max(Va_a(ii,jj+1),Va_a(ii,jj-1))
      Delta_T_FA(ii,jj)= Delta_L_FA(ii,jj)/Av_Va        !    Transit time across Delta_L segment
      FLR_Per(ii) = FLR_Per(ii) + Delta_T_FA(ii,jj)     !    Transit time along field line
    enddo
  enddo

  Delta_L_FA(:,Num_u3) = Maxval(Delta_L_FA)
  Delta_T_FA(:,Num_u3) = Maxval(Delta_T_FA)
  Delta_L_FA(:,1) = Maxval(Delta_L_FA)
  Delta_T_FA(:,1) = Maxval(Delta_T_FA)

!  print*,Maxval(Delta_L_FA),Maxval(Delta_T_FA)

  allocate (Delta_L_CF(Num_u1,Num_u3))
  allocate (Delta_T_CF(Num_u1,Num_u3))

  do ii = 2, Num_u1-1
    do jj = 1, Num_u3
!   Across field line
      Delta_x = (xArr(ii+1,jj) - xArr(ii-1,jj))   !  Change in X coordinate along field line
      Delta_y = (yArr(ii+1,jj) - yArr(ii-1,jj))   !  Change in Y coordinate along field line
      Length = sqrt(Delta_x**2 + Delta_y**2)      !  Line segment Length (in Re)
!      Av_Va = (Va_a(ii+1,jj) + Va_a(ii-1,jj))/2
      Av_Va =Max(Va_a(ii+1,jj),Va_a(ii-1,jj))             !  Average Alfven Speed across Delta_L segment (in m/s)
      Delta_L_CF(ii,jj)= Length                           !  Transit time across Delta_L segment
      Delta_T_CF(ii,jj)= Delta_L_CF(ii,jj)/Av_Va          !  Transit time across Delta_L segment
    enddo
  enddo

  Delta_L_CF(Num_u1,:) = Maxval(Delta_L_CF)
  Delta_T_CF(Num_u1,:) = Maxval(Delta_T_CF)
  Delta_L_CF(1,:) = Maxval(Delta_L_CF)
  Delta_T_CF(1,:) = Maxval(Delta_T_CF)

  print*,Minval(Delta_T_FA),Minval(Delta_T_CF)

  Courant_cond = min(Minval(Delta_T_FA),Minval(Delta_T_CF))/2
  print*,'Courant cond = ',courant_cond

  if (dt.GT.Courant_cond) then
    print*,"WARNING: Courant condition problem. Try reducing 'dt'"
    stop
  endif

  deallocate(Delta_T_CF)
  deallocate(Delta_L_CF)
  deallocate(FLR_Per)
!  deallocate(F_Len)
  deallocate(Delta_T_FA)
  deallocate(Delta_L_FA)

stop

!  allocate (D_out(Num_u1))
!  allocate (Out_DL(Num_u1-1))
!  allocate (h_mu_outer(Num_u3))
!
!  Bnd_L   = 0.0
!  D_out(1)   = 0.0
!
!  do ii = 1, Num_u1-1
!    Delta_x = (xArr(ii+1,Num_u3) - xArr(ii,Num_u3))     ! Change in X coordinate along outer Boundary
!    Delta_y = (yArr(ii+1,Num_u3) - yArr(ii,Num_u3))     ! Change in Y coordinate along outer Boundary
!    Out_DL(ii) = sqrt(Delta_x**2 + Delta_y**2)
!    D_out(ii+1) = D_out(ii) + Out_DL(ii)
!  enddo
!  Bnd_L = Sum(out_DL)                       ! Calculates length of outer boundary in Re
!
!  do ii = 1,Num_u3_2+1
!    h_mu_outer(2*ii-1) = h_mu(Num_u1,2*ii-1)
!  enddo
!
!! Calc basis functions, V and pV/pr
!  allocate (Ynm(K,Num_u1))
!  allocate (Ynm_s(K,Num_u1))
!  allocate (zeros(K))
!  call do_basis(m_num, k, Num_u1, LMin_u1, LMax_u1, del_Th, Ynm, Ynm_s, zeros)
!! Ynm and Ynm_S are (k,N1) arrays where the Legendres are listed from HIGH to LOW Lat
!! CoLat_N goes from LOW to HIGH Latitude i.e. HIGH to LOW CoLatitude
!!
!!  Set up scale factors for the Ionospheric Boundary Conditions for the Spherical Harmonic Expansion
!!   Need to handle the alternating grid pattern on the ionopsheric boundary....
!
!
!!  do kk = 1,K                    ! k = Number of Basis Functions in Expansion (from Sav file)
!! B3 and B2 are on the same grid in u1
!! South and North are the same
!    do ii = 1,Num_u1
!      Ynm_B3(kk,ii) = Ynm(kk,N1+1-ii)            ! Reverse order of Basis Vectors due to U1 direction (0 = low boundary N1-1 = upper boundary)
!      Ynm_B3_dr(kk,ii)=Ynm_S(kk,Num_u1+1-ii)
!      Ynm_g(kk,ii)=Ynm(kk,Num_u1+1-ii)*(2.d0*zeros(kk)+1.d0)/(zeros(kk)+1.d0)
!    enddo
!    do ii = 1,Num_u1_2
!      Ynm_B3_dr2(kk,ii) = Ynm_S(kk,Num_u1+1-(2*ii))  ! Reverse order and use 2nd value of Basis Vectors due to U1 direction (0 = low boundary N1-1 = upper boundary)
!    enddo
!  enddo
!
!  AlphaB3_2=MatMul(Ynm_B3_dr2,Transpose(Ynm_B3_dr2))
!  AlphaB3  =MatMul(Ynm_B3_dr, transpose(Ynm_B3_dr))
!
!  count1 = 0
!  time=0.d0
!
!filename=trim(BEdata)//'coordinate.dat'
!INQUIRE(FILE =filename, EXIST = exists)
!!	print*,	filename, '',exists
!
!if(restart <=0) then     ! to start a new run
!print*,''
!print*,'now we start a new run ....'
!         if(exists) then		 ! Whether the same data is already there....
!print*,' The data exists, Input ''C'' to continue, or quit!'
!      WRITE(*,'(A,$)')'  your command is:> '
!read(*,'(a)') cont
!if(cont/='C' .and. cont/='c')
!     +	 STOP 'Modify DATA_FILE or DIR and try again'
!print*,' So we proceed anyway......'
!   endif
!
!call getID(IDsave)
!      open(IDsave,file=filename,form='binary')
!write(IDsave) Dt*Plot_freq
!write(IDsave) (u1Arr(2*ii-1,1),ii=1,N1_2)
!write(IDsave) (u3Arr(1,2*jj),jj=1,N3_2)
!      close(IDsave)
!
!call getID(IDsave)
!      open(IDsave,file=trim(BEdata)//'fields.dat',form='binary')
!      write(IDsave) B1x_ord,B1y_ord
!      close(IDsave)
!
!call getID(IDsave)	! save the density, Alfven velocity, transit times.
!      open(IDsave,file=trim(BEdata)//'density.dat',form='binary')
!      write(IDsave) N1,N3,xArr,yArr,rho_a,Va_a,Delta_T_CF,Delta_T_FA
!      close(IDsave)
!!	print*,maxval(rho_a),minval(rho_a)
!
!call getID(IDsave)
!      open(IDsave,file=trim(BEdata)//'Parameter.log')
!      write(IDsave,'(i3,2x,i3,2x,a16)') N1_2,N3_2,'N1_2 and N3_2'
!      write(IDsave,'(i3,8x,a7)') M_num,'M_num'
!      write(IDsave,'(f6.2,5x,a7)') Width,'width'
!write(IDsave,'(f6.3,2x,a9)') freq,'freq'
!      close(IDsave)
!
!call getID(IDlog)
! Filename = TRIM(BEdata)//'Information.log'
! Open(IDlog,file=Filename,status='unknown')
!write(IDlog,'(A25)') 'This is a NEW run'
!!	write(IDlog,*) ''
!Nf=0

!else ! to restart,
!print*,''
!print*,' now we RESTART from a previous run...'
!
!         if(.not.exists) then
!   print*,'Error: File not exist, will quit... '
!   stop 'Modify DATA_FILE or DIR or RESTART and try again'
!   endif
!102print*,''
!print*,' Choose which record # to start with:'
!print*,' Negative value: start from the END of last run'
!print*,' Positive value: you specify the record #'
!print*,' 0 to QUIT !'
!
!      WRITE(*,'(A,$)')'   your command is:> '
!read(*,*) restart
!if(restart==0) stop ' QUIT !'
!
!if(restart<0) then
!count1=0
!
!do
!write(file2,'(i5.5)') count1
!filename=trim(BEdata)//'full.'//trim(file2)
!!	print*,filename
!INQUIRE(FILE =filename, EXIST = exists)
!if(.not.exists) exit
!count1=count1+1
!enddo
!if(count1<=1) then
!print*,' The record # does NOT exist'
!print*,' change RESTART and start a NEW run'
!stop ' QUIT !'
!endif
!count1=count1-1
!
!else
!
!write(file2,'(i5.5)') restart
!filename=trim(BEdata)//'full.'//trim(file2)
!INQUIRE(FILE =filename, EXIST = exists)
!if(.not.exists) then
!print*,'The record # does NOT exist, try another record #'
!goto 102
!endif
!count1=restart
!
!endif
!Nf=count1*full_freq+1
!
!write(file2,'(i5.5)') count1
!filename =trim(BEdata)//'full.'//trim(file2)
!call getID(IDsave)
!      open(IDsave,file=filename,status='old',form='binary')
!      read(IDsave,err=101) time
!      read(IDsave,err=101) B1_con,B2_con,B3_con,E1_con,E2_con,
!     +     B1_cov,B2_cov,B3_cov,E1_cov,E2_cov
!close(IDsave)
!
!print*,''
!if(Nf>=Nt-1)	stop ' Error: change NT and try again !'
!
!print*,' now we restart from #',count1,', time=',real(time)
!
!call getID(IDlog)
!write(file2,'(i5.5)') count1
! Filename = TRIM(BEdata)//'Information.log.'//trim(file2)
! Open(IDlog,file=Filename,status='unknown')
!write(IDlog,'(A36,2x,i4,a8,f6.2,a7)') 'This is a RESTART run, from
!     + Record #',count1,', time=',time,'second'
!1write(IDlog,*) ''
!
!endif
!
!call get_time(tfrom,0,IDlog)
!
!Print*,'Output task info to ', Filename
!print*,''
!print*,'    Record #   Iteration #     Time'
!
! 
!!   *******************************************************************************************************
!!   Start Time iteration of Fields
!!   ******************************************************************************************************
!
!  D_u1(:,:)=u1Arr(IR(:),:)-u1Arr(IL(:),:)!   Set up differences on Grid
!  D_u3(:,:)=u3Arr(:,JR(:))-u3Arr(:,JL(:))
!!
!!   Time Driver on Equatorial Slice (constant u1 surface)
!!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!   Set E = div(T) for an Alfven wave T = is a guassian
!
!!    zz0    = D_out(Fix(N1/2.0+9))
!!    QQ1 = where (Lvarr GE Posn_S)
!!    zz0 = D_out(min(QQ1))
!    do ii=1,N1
!      if((Lvarr(ii)<=Posn_S).and.(Lvarr(ii+1)>Posn_S)) zz0=D_out(ii+1)
!    enddo
!
!    Do ii = 2,N1-1
! zz       = (D_out(ii)-zz0)/Hwidth
!  Tdr(ii)  = Exp(-zz**2.)
!    enddo
!    Tdr(1)=0.0 ; Tdr(N1) = 0.0
!
! Do ii = 2,N1-1
!        D_Tdr(ii) = (Tdr(ii+1) - Tdr(ii-1))/(D_out(ii+1)-D_out(ii-1))       ! Derivative of Gausssiona in U1 direction
! Enddo
!    D_Tdr(1)=0.0 ; D_Tdr(N1) = 0.0

! Do tt = Nf,Nt-1                             ! Start of Time Iteration Loop
!
! Time     = 2.d0*Dt*tt                              ! Time Counter
!
!!   Boundary Conditions
!!   ~~~~~~~~~~~~~~~~~~~~
!
!!   Boundary Conditions
!!   ***************************************************************
!!   Outer Bounday
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  E1_cov(0:N1-1,N3-1) = E1_Dr(0:N1-1)
!  E2_cov(0:N1-1,N3-1) = E2_Dr(0:N1-1)


!   Single Frequency driver (with expondential evelope)
! StartRamp = (Exp(1.d0-T_ramp/(time))/EXP(1.d0))
! EndRamp   = (Exp(1.d0-T_ramp/(Dstop-time))/EXP(1.d0))

!   Ramp on driver to mimimize transients
! Tdep   = 1.0e-2*sin(2.d0*dpi*freq/1000.d0*time)*StartRamp*EndRamp

!   Zero driver evelope
!If Time GT Dstop then Tdep = 0.0

!   Pulse driver
! Norm      = 0.2
! Tdep   = 10.0*Time/Norm*(1.0-time/Norm)*Exp(-5.0*time/Norm)


! carr    = cos(2.*dpi*(time-2.5/Width)*freq)
! Env     = exp(-((time-2.5/Width)*Width)**2.)
! Tdep    = Carr*Env
!
!**E1_outer=-B2_cov(:,N3)*Va_a(:,N3)*evens +D_Tdr*Tdep*evens*Amp        ! Allow outward travelling waves ?
!**E2_outer= B1_cov(:,N3)*Va_a(:,N3)*odds+m_num*Tdr*Tdep*odds*Amp   ! Allow outward travelling waves ?
!
!!   Conducting Northern Ionosphere
!!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  E1_cov(:,1)    = E1_N(1:N1)
!  E2_cov(:,1)    = E2_N(1:N1)                ! Conducting Northern Ionosphere using BC on Ionosphere
!
!!   Inner L shell Boundary
!!   ~~~~~~~~~~~~~~~~~~~~~~
!  E2_cov(1,:)      = 0.d0                    ! Perfectly Reflecting
!  B1_cov(1,:)      = 0.d0
!
!!   Outer L shell Boundary
!!   ~~~~~~~~~~~~~~~~~~~~~~
!   B3_cov(N1,:)   =  0.0 
!
!!   *************************************************************************************************************
!   Start time iterations
!
!!   E field iteration
!!
!!   Solve E1_con
!      E1_con(:,:) = 2.d0*Dt*V2(:,:)/Jacb(:,:)*   
!     +  (zim*M_num*B3_cov(:,:) - (B2_cov(:,JR(:)) - 
!     +  B2_cov(:,JL(:)))/D_u3(:,:)) + E1_con(:,:)
!
!!   Solve E2_con
!      E2_con(:,:) = 2.d0*Dt*V2(:,:)/Jacb(:,:)*((B1_cov(:,JR(:)) 
!     +  - B1_cov(:,JL(:)))/D_u3(:,:)  
!     +  - (B3_cov(IR(:),:) - B3_cov(IL(:),:))
!     +  /D_u1(:,:)) + E2_con(:,:)
!
!!   End of E field contravariant iteration
!
!!   Evolving Contravariant (Tangent) vector fields to Covariant vectors fields
!!   **************************************************************************
!
!!   Evolve E1_cov
!      E1_cov = E1_con/g11_con
!
!!   Evolve E2_cov
!      E2_cov = E2_con*g22_cov
!
!!   Conducting Northern Ionosphere
!!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  E1_cov(:,1)         = E1_N(:)
!  E2_cov(:,1)         = E2_N(:)               ! Conducting Northern Ionosphere using BC on Ionosphere
!
!!   Outer Boundary (with Driver)
!!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  E1_cov(:,N3)   = E1_Outer(:)
!  E2_cov(:,N3)   = E2_Outer(:) 
!!   Inner L shell Boundary
!!   ~~~~~~~~~~~~~~~~~~~~~~
!  E2_cov(1,:)= 0.d0                           ! Clean up after shifts due to odd number of points
!  E2_con(1,:)= 0.d0 ;   E2_con(:,1) = 0.d0
!  B1_cov(1,:)= 0.d0
!  E1_con(1,:)= 0.d0 ;   E1_con(:,1) = 0.d0
!
!!   B field iteration
!!   ********************************************************** *******!
!
!!   Solve B1_con
!      B1_con(:,:)= 2.d0*Dt/Jacb(:,:)*(E2_cov(:,JR(:))   
!     +  - E2_cov(:,JL(:)))/D_u3(:,:) + B1_con(:,:)
!
!!   Solve B2_con
!      B2_con(:,:)= -2.d0*Dt/Jacb(:,:)*(E1_cov(:,JR(:)) 
!     +  - E1_cov(:,JL(:)))/D_u3(:,:) + B2_con(:,:)
!
!!   Solve B3_con
!      B3_con(:,:)=2.d0*Dt/Jacb(:,:)*( zim*m_num*
!     +   E1_cov(:,:) - (E2_cov(IR(:),:) -
!     +   E2_cov(IL(:),:))/D_u1(:,:)) + B3_con(:,:)
!
!!   Outer Corner Point Ionopshere and Outer L shell
!      B3_con(N1,1) = 4.d0/3.d0*B3_con(N1-2,1) -1.d0/3.d0*B3_con(N1-4,1)
!      B3_con(N1,N3)= 4.d0/3.d0*B3_con(N1-2,N3)-1.d0/3.d0*B3_con(N1-4,N3)
!!   End of B field contravariant interation
!
!!   Evolving Contravariant Interior Grid B fields to Convariant vectors fields
!!   **************************************************************************
!
!!   Evolve B1_cov
!      Ave_B3(:,:)= (B3_con(IR(:),JR(:)) + B3_con(IR(:),
!     +   JL(:)) + B3_con(IL(:),JR(:)) + B3_con(IL(:),JL(:)))/ 4.d0
!      Ave_B3(:,1)= 0.0;    Ave_B3(:,N3) = 0.0
!      B1_cov      = B1_con*g11_cov + Ave_B3*g13_cov
!      B1_cov(:,1) = 0.d0  ;   B1_cov(:,N3)    = 0.d0         ! Clean up after shifts due to odd number of points
!
!!   Inner L shell Boundary
!!   ~~~~~~~~~~~~~~~~~~~~~~
!      B1_cov(1,:)        = 0.d0
!
!!   Evolve B2_cov
!      B2_cov = B2_con*g22_cov
!      B2_cov(:,1)     = 0.d0  ;     B2_cov(:,N3)   = 0.d0   ! Clean up after shifts due to odd number of points
!
!!   Evolve B3_cov
!      Ave_B1(:,:) = (B1_con(IR(:),JR(:)) + B1_con(IR(:),
!     +   JL(:)) + B1_con(IL(:),JR(:)) +
!     +   B1_con(IL(:),JL(:)))/ 4.d0
!      Ave_B1(:,1)  = 0.0      ;    Ave_B1(:,N3)   = 0.0      ! Clean up after shifts due to odd number of points
!      Ave_B1(:,2)  = 0.0      ;    Ave_B1(:,N3-1) = 0.0
!      B3_cov  = Ave_B1*g13_cov + B3_con*g33_cov
!
!!   Along Northern Ionospheric Boundary
!      Ave_B1(:,:) = (B1_con(IL(:),JR(:)) + B1_con(IR(:)
!+   ,JR(:)) )/2.d0
!      B3_cov(:,1) = Ave_B1(:,1)*g13_cov(:,1) + B3_con(:,1)*g33_cov(:,1)
!!   Along Southern Ionosphere Boundary
!      Ave_B1(:,:) = (B1_con(IL(:),JL(:)) + B1_con(IR(:),
!     +   JL(:)) )/2.d0
!      B3_cov(:,N3) = Ave_B1(:,N3)*g13_cov(:,N3) + B3_con(:,N3)*
!     +   g33_cov(:,N3)
!
!!   Along Inner L shell Boundary
!      Ave_B1(:,:) = (B1_con(IR(:),JR(:)) + B1_con(IR(:),
!     +   JL(:)) )/2.d0
!      B3_cov(1,:) = Ave_B1(1,:)*g13_cov(1,:) + B3_con(1,:)*g33_cov(1,:)
!
!!   Along Outer L shell Boundary
!      B3_cov(N1,:) = 0.0

!   Ionospheric Boundary Condition
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!   Using B3_con? in Ionosphere and perfectly conducting ground (B3_con = 0) solve Laplaces equation in Neutral Atmopshere
!   Calculate B1_cov and B2_cov just below the Ionospheric current sheet.
!!   Calculate jump in B1 and B2 and calculate Current in Sheet
!!   Invert to find E1 and E2 values in sheet and use as BC in next iteration.
!!   Calculate Coefficients in expansion of Bz at the ionosphere (with derivative radial component scale at r = Ri)
!
!!   At ionosphere boundary:
!!   Contravariant at ionos boundary for BB3 is radial
!!   Covariant at ionos boundary for BB1(colat) and BB2(azimuth) is in the sheet
!
!!   FOR NORTHERN HEMISPHERE
!      B3_N1        = (B3_con(:,1))*h_ra_n                         ! 1,3,5 etc. of full array
!!  B3_N1     = (Shift(B3_N1,1) + Shift(B3_N1,-1))/2.0 + B3_N1     ! Linear Interpolation
!  Do ii = 1,N1_2
! B3_N(ii) = B3_N1(2*ii)                    ! Use only point directly calculated
!  enddo
!! Fit Br data to pV/pr basis set
!
!!  DataN        = transpose(B3_N)
!  DataN=B3_N
!  BetaN        = MatMul(DataN,transpose(DConjG(Ynm_B3_dr2)))
!
!!  Coeffs_B3_N   = LU_Complex(AlphaB3,BetaN)      ! Coeffs from d(psi)/dr
!!	  call DLinct(K,AlphaB3,K,1,Alphatemp,k)
!!	  Coeffs_B3_N   =MatMul(Alphatemp,BetaN)
!!      call Dlslcg(K,AlphaB3,K,BetaN,1,Coeffs_B3_N)
!      call Dlsacg(K,AlphaB3,K,BetaN,1,Coeffs_B3_N)
!  DThi_dr_N     = MatMul(transpose(Ynm_B3_dr),Coeffs_B3_N)    ! This should be B3_N at twice the spatial resolution
!!  DThi_dr_N     = reform(DThi_dr_N)
!      Thi_N        = MatMul(transpose(Ynm_B3),Coeffs_B3_N)          ! Potential Function
!
!! Calc. pV/p_Theta
!  Do jj = 2,N1-1
! DThi_dTh_N(jj) = (Thi_N(jj+1)-Thi_N(jj-1))/(-2.d0*del_Th)
!  enddo
!! Do the end points
!  DThi_dTh_N(1)  = (-3.0*Thi_N(1)+ 4.d0*Thi_N(2)- Thi_N(3))
!+      /(-2.d0*del_th)
!  DThi_dTh_N(N1) = ( 3.0*Thi_N(N1) - 4.d0*Thi_N(N1-1) + 
!     +     Thi_N(N1-2))/(-2.d0*del_th)
!  DThi_dTh_N     = 1.d0/Ri*DThi_dTh_N                              ! B_Theta Atmos
!  DThi_dPh_N     = DCmplx(0.0,1.d0)*M_num*Thi_N/(Ri*sin(Colat_N))  ! B_phi   Atmos, CoLat_N runs from INNER L shell to OUTER L shell
!
!! Interpolate B1 and B2 to ALL points just above the ionosphere from model solution to calculate currents in Ionosphere
!  B1_N  = B1_cov(:,2)/h_th_n
!  B2_N  = B2_cov(:,2)/h_ph_n
!!  B1_N  = (Shift(B1_N,1) + Shift(B1_N,-1))/2.0 + B1_N          ! Linear Interpolation
!!  B2_N  = (Shift(B2_N,1) + Shift(B2_N,-1))/2.0 + B2_N          ! Linear Interpolation
!  tempN(:)=(B1_N(IL(:))+B1_N(IR(:)))/2.
!  B1_N=tempN+B1_N
!  tempN(:)=(B2_N(IL(:))+B2_N(IR(:)))/2.
!  B2_N=tempN+B2_N 
!	print*,Real(B1_N(1:3)),Real(B2_N(1:3))
! 
!
!!   Try second order Taylor expansion for end points
!  h          = CoLat_N(N1-1)-CoLat_N(N1-3)                            ! step size along ionosphere is uniform on B1_N grid
!  FirD  = ( B1_N(N1-5) - 4.0*B1_N(N1-3) + 3.0*B1_N(N1-1))/(2.0*h)     ! First Derivative Backwards difference O(h**2)
!  SecD       = (-B1_N(N1-7) + 4.0*B1_N(N1-5) - 5.0*B1_N(N1-3) +
!     +       2.0*B1_N(N1-1))/(h**2)                                         ! Second Derivative Backwards difference O(h**2)
!  h1         = CoLat_N(N1-1)-CoLat_N(N1)                              ! Step size on B1_n_interp grid
!  B1_N(N1)   = B1_N(N1-1) + h1*FirD + ((h1**2)/2.0)*SecD
!
!  h          = CoLat_N(4)-CoLat_N(2)                                  ! step size along ionosphere is uniform
!  FirD         = (-B2_N(6) + 4.0*B2_N(4) - 3.0*B2_N(2))/(2.0*h)       ! First Derivative Forwards difference O(h**2)
!  SecD= (-B2_N(8) + 4.0*B2_N(6) - 5.0*B2_N(4)+2.0*B2_N(2))/(h**2)     ! Second Derivative Forwards difference O(h**2)
!  h1         = CoLat_N(1)-CoLat_N(2)
!  B2_N(1)     = B2_N(2) + h1*FirD + ((h1**2)/2.0)*SecD 
!
!!   Currents in the Ionosphere
!  J_Th_N    = -(B2_N - DThi_dPh_N)/u0
!  J_Ph_N    =  (B1_N - DThi_dTh_N)/u0
!
!!   Calculate electric fields in Ionosphere from discontinuity in B's for next time step...
!  E_Th_N    = (INVSS_N(:,1,1)*J_Th_N + INVSS_N(:,2,1)*J_Ph_N)
!  E_Ph_N    = (INVSS_N(:,1,2)*J_Th_N + INVSS_N(:,2,2)*J_Ph_N)
!
!**  E1_N   = E_Th_N*h_th_n*evens
!**  E2_N   = E_Ph_N*h_ph_n*odds
!
!! perfect reflect
!! E1_N(0:N1-1)  = 0.0
!! E2_N(0:N1-1)  = 0.0
!!
!! Now do the Ground fields
!  Thi_gnd_N=MatMul(transpose(Ynm_g),Coeffs_B3_N)                   ! Potential Function
! Do jj = 2,N1-1
!DThi_dTh_N_g(jj)=(Thi_gnd_N(jj+1)-Thi_gnd_N(jj-1))/(-2.d0*del_Th)
!  enddo
!! Do the end points
!  DThi_dTh_N_g(1)  = (-3.0*Thi_gnd_N(1)  + 4.d0*Thi_gnd_N(2)
!     +     - Thi_gnd_N(3))/(-2.d0*del_th)
!  DThi_dTh_N_g(N1) = ( 3.0*Thi_gnd_N(N1) - 4.d0*Thi_gnd_N(N1-1)
!     +     + Thi_gnd_N(N1-2))/(-2.d0*del_th)
! DThi_dTh_N_g = 1.d0/Ri*DThi_dTh_N_g                               ! B_Theta Ground
! DThi_dPh_N_g = DCmplx(0.0,1.d0)*M_num*Thi_gnd_N/(Re*sin(Colat_N)) ! B_phi   Ground
!
!!   Output Routines
!
!if(mod(tt,full_freq) == 0) then	 ! to save the whole data set for future restart
!count1=tt/full_freq
!call getID(IDsave)
!write(file2,'(i5.5)') count1
!filen=trim(BEdata)//'full.'//trim(file2)
!open(IDsave,file=filen,form='binary')
!      write(IDsave) time,B1_con,B2_con,B3_con,E1_con,E2_con,
!     +     B1_cov,B2_cov,B3_cov,E1_cov,E2_cov
!      close(IDsave)
!endif
!
!if(mod(tt,PLOT_freq) == 0) then
!count=tt/Plot_freq
!
!Ave_B2(:,:)=(B2_con(IL(:),JR(:))+B2_con(IL(:),
!     + JL(:))+B2_con(IR(:),JR(:))+B2_con(IR(:),JL(:)))/4.
!Ave_B3(:,:)=(B3_cov(IL(:),JR(:))+B3_cov(IL(:),
!     + JL(:))+B3_cov(IR(:),JR(:))+B3_cov(IR(:),JL(:)))/4.
!Ave_E1(:,:)=(E1_con(IL(:),JR(:))+E1_con(IL(:),
!     + JL(:))+E1_con(IR(:),JR(:))+E1_con(IR(:),JL(:)))/4.
!Ave_E2(:,:)=(E2_con(IL(:),JR(:))+E2_con(IL(:),
     !+ JL(:))+E2_con(IR(:),JR(:))+E2_con(IR(:),JL(:)))/4.
!
!      DO ii = 1,N1_2       ! uses all cells in U1 direction
!      Do jj = 1,N3_2       ! uses all cells in U3 direction
!      b_nu(ii,jj) = h_nu(2*ii-1,2*jj)*B1_con(2*ii-1,2*jj)/Re_m     ! dB in nT        ! B in T
!      b_ph(ii,jj) = h_ph(2*ii-1,2*jj)*Ave_B2(2*ii-1,2*jj)/Re_m     ! dB in nT
!      b_mu(ii,jj) = Ave_B3(2*ii-1,2*jj)/h_mu(2*ii-1,2*jj)/Re_m     ! dB in nT
!      e_nu(ii,jj) = h_nu(2*ii-1,2*jj)*Ave_E1(2*ii-1,2*jj)          ! dE in V/m
!      e_ph(ii,jj) = h_ph(2*ii-1,2*jj)*Ave_E2(2*ii-1,2*jj) ! dE in V/m
!       enddo
!       enddo 
!
!if(mod(tt,200)==0) print*,count,tt,real(time)
!
!!save,Filename=Plot_dir+'BE'+string(count,'(i4.4)')+'.sav',b_mu,b_nu,b_ph,e_nu,e_ph
!!save,Filename=Plot_dir+'BE'+string(count,'(i4.4)')+'.sav',B1_con,Ave_B2,Ave_B3,Ave_E1,Ave_E2 !b_mu,b_nu,b_ph,e_nu,e_ph
!!endif
!
!
!!  If tt Mod PLOT_freq eq 0 then $
!!  Begin
!
!!Print*,'Time iteration Number = ',tt,' Time = ',Time
!! Calc field aligned coord field arrays
!!   Re_m=6378000.d0                ! Earth radii in m
!   Do ii = 1,N1_2    ! uses all cells in U1 direction
!    Do jj = 1,N3_2   !    uses all cells in U3 direction
!       b_nu_a(ii,jj) = 1.0e9*h_nu(2*ii-1,2*jj)*B1_con(2*ii-1,2*jj)/Re_m        ! dB in nT
!       b_ph_a(ii,jj) = 1.0e9*h_ph(2*ii,2*jj)*B2_con(2*ii,2*jj)/Re_m
!       b_mu_a(ii,jj) = 1.0e9*B3_cov(2*ii,2*jj-1)/h_mu(2*ii,2*jj-1)/Re_m
!       e_nu_a(ii,jj) = 1.0e3*h_nu(2*ii,2*jj-1)*E1_con(2*ii,2*jj-1)             !dE in mV/m
!       e_ph_a(ii,jj) = 1.0e3*h_ph(2*ii-1,2*jj-1)*E2_con(2*ii-1,2*jj-1)
!       b_nu_a(ii,jj) = h_nu(2*ii-1,2*jj)*B1_con(2*ii-1,2*jj)        ! dB in nT
!       b_ph_a(ii,jj) = h_ph(2*ii,2*jj)*B2_con(2*ii,2*jj)
!       b_mu_a(ii,jj) = B3_cov(2*ii,2*jj-1)/h_mu(2*ii,2*jj-1)
!       e_nu_a(ii,jj) = h_nu(2*ii,2*jj-1)*E1_con(2*ii,2*jj-1)             !dE in mV/m
!       e_ph_a(ii,jj) = h_ph(2*ii-1,2*jj-1)*E2_con(2*ii-1,2*jj-1)
!       enddo
!enddo
!!   For ii = 0,N1_2-1 do $   ! uses all cells in U1 direction
!!    Begin
!!     b_mu_a(ii,0) = 1.0e9*B3_cov(2*ii+1,0)/h_mu(2*ii+1,0)/Re_m
!!     e_nu_a(ii,0) = 1.0e3*h_nu(2*ii+1,0)*E1_con(2*ii+1,0)
!!     e_ph_a(ii,0) = 1.0e3*h_ph(2*ii,0)*E2_con(2*ii,0)
!
!!     b_mu_a(ii,N3_2) = 1.0e9*B3_cov(2*ii+1,N3-1)/h_mu(2*ii+1,N3-1)/Re_m
!!     e_nu_a(ii,N3_2) = 1.0e3*h_nu(2*ii+1,N3-1)*E1_con(2*ii+1,N3-1)
!!     e_ph_a(ii,N3_2) = 1.0e3*h_ph(2*ii,N3-1)*E2_con(2*ii,N3-1)
!!    end
!
!
!!
!! Contravariant at ionos boundary for BB3 is radial
!! Covariant at ionos boundary for BB1(colat) and BB2(azimuth) is in the sheet
!!
!
!   BB1_n(:,count) = B1_n
!   BB2_n(:,count) = B2_n
!   BB3_n(:,count) = Dthi_dr_n
!   EE1_n(:,count) = E_Th_N
!   EE2_n(:,count) = E_Ph_N
!   JJ1_n(:,count) = J_Th_N
!   JJ2_n(:,count) = J_Ph_N
!
!   BB1_n_gnd(:,count)  = DThi_dTh_N_g
!   BB2_n_gnd(:,count)  = DThi_dPh_N_g

!    Time_Str = StrTrim(string(format='(f7.2)',Count1),1)

!call getID(IDsave)
!write(file2,'(i5.5)') count
!filen=trim(BEdata)//trim(file2)//'.dat'
!!print*,	'BE fields data:'; print*,'  ',filen
!open(IDsave,file=filen,form='binary')
!write(IDsave) time,b_mu,b_nu,b_ph,e_nu,e_ph
!      close(IDsave)
!      close(IDsave)

!	 count1=count1+1
!	 count=count+1

!FilenameM = Plot_dir+'Sheet_IAR_'+Time_Str+'.sav'
! Save,Filename=FilenameM,T_count,B1_con,B2_con,B3_con,E1_con,E2_con,B1_cov,B2_cov,B3_cov,E1_cov,E2_cov
! Save,Filename=FilenameM,T_count,e_nu_a,e_ph_a,b_nu_a,b_ph_a,b_mu_a,$
!                  B1_n,B2_n,Dthi_dr_n,$
!                  E_Th_N,E_Ph_N,J_Th_N,J_Ph_N,DThi_dTh_N_g,DThi_dPh_N_g,$
!                  DThi_dTh_N,DThi_dPh_N,Coeffs_B3_N,time,tt

! Print*,'Time = ',Count1,' seconds'
! count1 =Fix(count1+1)
! end                !    end of modulo plot loop
!endif
!
! Enddo
!
!!   End of time iteration loop ****************************************************************************************
!
!
!! 7777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777
!
! !FilenameI = Plot_dir+Ionofile
! !Save,Filename=FilenameI,T_count,Dt,Plot_freq,N1,Colat_N,LVArr,flr_a,Sd_N,Sp_N,Sh_N,$
! !                       BB1_n,BB2_n,BB3_n,BB3_s,EE1_n,EE2_n,JJ1_n,JJ2_n,Va_a
!
! !FilenameG = Plot_dir+gnd_file
! !Save,Filename=FilenameG,T_count,Dt,Plot_freq,N1,Colat_N,LVArr,flr_a,$
! !                       BB1_n_gnd,BB2_n_gnd
!
! !FilenameP = Plot_dir+Data_file
! !Save,Filename=FilenameP,T_count,Dt,Plot_freq,N1,N3,Colat_N,LVArr,Va_a,m_num,del_Th_0_1,$
! !                 Nt,Dt,Freq,Re,Dstop,Amp,T_ramp,LMin_1,LMax_1,$
! !                 g11_con,g11_cov,g22_con,g22_cov,g13_con,g13_cov,g33_con,g33_cov,$
! !                 h_mu,h_nu,h_ph,Sd_N,Sp_N,Sh_N,$
! !                 Ynm,Ynm_S,zeros,F_line,$
! !                h_ra_n,h_th_n,h_Ph_n,$
! !                  Delta_T_FA,Delta_T_CF,courant_cond,Intstarttime,Starttime,Finishtime,$
! !                xArr,yArr,Tharr,Earth,Xe,Ye,B3x_ord,B3y_ord,B2x_ord,B2y_ord,B1x_ord,B1y_ord,$
! !                E2x_ord,E2y_ord,E1x_ord,E1y_ord,F_len,u3arr,u1arr,rho_a,Tdr
!
!Print*,'  Finished computation ....'
!
!write(IDlog,'(A90)')'Two Dimensional Time Dependent Meridional Sli
!     +ce - Ideal MHD Model with Ionosphere'
!write(IDlog,*) ''
!write(IDlog,'(A48,2X,F10.5)') 'Courant Conditon = ',Courant_cond
!write(IDlog,*) ''
!write(IDlog,'(A60)') 'Model Run Information'
!write(IDlog,'(A48,2X,i5)') 'Azimuthal Wave Number m = ',M_num
!write(IDlog,'(A48,2X,i5)') 'Number of full timesteps Nt = ',Nt
!write(IDlog,'(A48,2X,F10.5)') '1/2 Time step interval Dt = ',Dt
!write(IDlog,'(A48,2X,F10.5)') 'Freq in mHz of the Driver = ',Freq
!write(IDlog,'(A48,2X,F10.2)') 'Time Driver stops Dstop = ',Dstop
!write(IDlog,'(A48,2X,E10.5)') 'Amplitude of Driver Amp = ',Amp
!write(IDlog,'(A48,2X,F10.5)') 'Height of Ionosphere (in Re) = ',RI
!write(IDlog,*) ''
!write(IDlog,'(A55)') 'Grid Information'
!write(IDlog,'(A48,2X,i5)') 'Number of field lines, u1 coords N1
!     + = ',N1
!write(IDlog,'(A48,2X,i5)') 'Number of points along field lines N3
!     + = ',N3
!write(IDlog,*) ''
!write(IDlog,'(A48,2X,F10.5)') 'Max Colatitude (deg) = ',
!     +  Maxval(Colat_N)*180./dpi
!write(IDlog,'(A48,2X,F10.5)') 'at outer L shell (in Re) = ',LMax_1
!write(IDlog,*) ''
!write(IDlog,'(A48,2X,F10.5)') 'Min Colatitude (deg) = ',
!     +  Minval(Colat_N)*180./dpi
!write(IDlog,'(A48,2X,F10.5)') 'at inner L shell (in Re) = ',LMin_1
!write(IDlog,*) ''
!write(IDlog,'(A48,2X,F10.5)') 'Cap Angle (deg) = ',
!+  (Maxval(Colat_N)-Minval(Colat_N))*180./dpi
!write(IDlog,'(A48,2X,F10.5)') 'Angle step size = ',del_Th_0_1
!write(IDlog,'(A48,2X,F10.3)') 'Spacial resolution at Ionosphere (i
!	+n M) = ',del_Th_0_1*dpi/180.*RI*Re_m
!write(IDlog,*) ''
!cwrite(IDlog,'(A50)') 'Real Ionopsheres from MSIS IRI and IGRF'
!cwrite(IDlog,'(A48,2X,A48)') 'Northern Ionosphere File = ',
!c     +  Ionosphere_N
!cwrite(IDlog,*) ''
!write(IDlog,'(A48,2X,i5)') 'Number of time steps between plots Plo
!     +T_freq = ',Plot_freq
!write(IDlog,'(A48,2X,F10.5)') 'Frequency of Plots = ',
!     + Plot_freq*2.*Dt
!write(IDlog,*) ''
!
! call get_time(tto,1,IDlog)
!
! Print*,'Total time taken = ',Real(tto - tfrom),' seconds'
! Print*,'Average Time per time iteration = ',
!     +        Real((tto - tfrom)/(Nt-Nf))
! Print*,'Courant Conditon = ',Real(Courant_cond)

!write(IDlog,'(A48,2X,F16.9)') 'Total  Time taken = ',Finishtime - Starttime
!write(IDlog,'(A48,2X,F16.9)') 'Total  Time (in Hrs) = ',(Finishtime - Starttime)/3600.
!write(IDlog,'(A48,2X,F16.9)') 'Average Time per time iteration = ',(Finishtime - Intstarttime)/Nt
!write(IDlog,'(A30,2x,a40)')'Data were saved in:',Dir
!write(IDlog,'(A30,2x,a10)')'with the name:',Data_file
!write(IDlog,*) ''
!close(IDlog)

!	Print*,'Finished Info file'

!	Print*,''
!	Print*,'Courant Conditon = ',Courant_cond
!	Print*,'Max Colatitude (deg) = ',Maxval(Colat_N)*180./dpi
!	Print*,'at outer L shell (in Re) = ', LMax_1
!	Print*,'Min Colatitude (deg) = ',Minval(Colat_N)*180./dpi
!	Print*,'at inner L shell (in Re) = ', LMin_1
!	Print*,'Cap Angle (deg) = ',(Maxval(Colat_N)-Minval(Colat_N))*
!     +    180./dpi
!	Print*,'Angle step size (deg) = ',del_Th_0_1
!	Print*,'Spatial resolution long Ionosphere (in M) = ',
!     +     del_Th_0_1*dpi/180.*RI*Re_m
!	Print*,''
!
!!jump1:continue=1

!stop
!101 print*,' Read old data wrong! '

END                  !     End of program
!
!!
!! ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
!!
!
!
!subroutine get_time(tt,Iflag,ID)
!integer::time_array(9),Iflag,IDsave,ID
!character*10::Ttime,tempfile,file2
!	Real(4)::second
!double precision::tt
!tempfile='Temp.fil'
!
! call Date_and_Time(time=Ttime,values=time_array)
! write(file2,12) time_array(1),time_array(2),time_array(3)
!12     format(i4.4,'.',i2.2,'.',i2.2)
!       print*,''
!if(Iflag==0) then
!write(*,*) 'Task starts from: ',trim(file2)//' '//Ttime(1:2)
!     +    //':'//Ttime(3:4)//':'//Ttime(5:10)
!if(ID > 0) then
!write(ID,*)''
!write(ID,*) 'Task starts from: ',trim(file2)//' '//Ttime(1:2)
!     +    //':'//Ttime(3:4)//':'//Ttime(5:10)
!write(ID,*)''
!endif
!else
!write(*,*) 'Task ends at: ',trim(file2)//' '//Ttime(1:2)
!     +    //':'//Ttime(3:4)//':'//Ttime(5:10)
!if(ID > 0) then
!write(ID,*)''
!write(ID,*) 'Task ends at: ',trim(file2)//' '//Ttime(1:2)
!     +    //':'//Ttime(3:4)//':'//Ttime(5:10)
!write(ID,*)''
!endif
!endif
!print*,''
!
!call getID(IDsave)
!open(IDsave,file=tempfile)          !,status='unknown')
!write(IDsave,53) Ttime(5:10)
!53 format(a6)
!close(IDsave)
!open(IDsave,file=tempfile)
!read(IDsave,54) second
!54format(f6.3)
!close(IDsave)
!tt=3600.*time_array(5)+60.*time_array(6)+second
!
!return
!end subroutine get_time
!
!!
!! ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
!!
!
!
!subroutine getID(ID)
!integer::ID
!    	LOGICAL::busy
!character*8:: read,write,readwrite
! 
!ID=20
!do
!       inquire (unit=ID, read=read, write=write, readwrite=readwrite)
! inquire (unit=ID,opened=busy)
!       if (read == "UNKNOWN" .and. write == "UNKNOWN" 
!     +     .and. readwrite == "UNKNOWN" .and. .not.busy) return
!        ID = ID + 1
!enddo
!
!  return
!end subroutine getID

!  include 'dgeevx.f'
