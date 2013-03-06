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

    integer ::  ii, jj, kk
    integer ::  Num_u1, Num_u1_2, Num_u3, Num_u3_2, K
    integer(8) :: nt, tt
    integer :: plot_freq
    real(DBL) :: data_sp

    real(DBL) :: LMin_u1, LMax_u1, LMin_u3, LMax_u3, m_num
    real(DBL) :: dt, del_Th_0_u1, del_Th

! Va vars
    integer :: i1, i2
    real(DBL) :: xd, va_sp, lat, bmg, PLaw, nH, nO, nO_Ht
    real(DBL), allocatable, dimension(:,:) :: V2, Va_a

! grid arrays
    real(DBL), allocatable, dimension(:) :: LVArr, CoLat_N, CoLat_S, h_mu_outer
    real(DBL), allocatable, dimension(:,:) :: u1Arr, u3Arr, RArr, ThArr, &
                                              D_u1, D_u3
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
    real(DBL), allocatable, dimension(:,:) :: Ynm, Ynm_s, Ynm_B3, Ynm_g, &
                                              Ynm_B3_dr2, Ynm_B3_dr
    real(DBL), allocatable, dimension(:) :: zeros
    complex(DBL), allocatable, dimension(:,:) :: AlphaB3

!   loop variables
    real(DBL), allocatable, dimension(:) :: Tdr, drv_B3
    real(DBL) :: time, h_width, zz, bfa
    integer :: count1

!   shift arrays
    integer, allocatable, dimension(:) :: S1L, S1R, S3L, S3R

!   loop arrays
    integer ::LDB,  INFO
    real(DBL) :: h, h1
    complex(DBL) :: FirD, SecD
    integer, allocatable, dimension(:) :: IPIV
    complex(DBL), allocatable, dimension(:,:) :: B1_con, B2_con, B3_con, &
                                                 B1_cov, B2_cov, B3_cov, &
                                                 E1_con, E2_con, &
                                                 E1_cov, E2_cov, &
                                                 Av_B1, Av_B3
    complex(DBL), allocatable, dimension(:) :: B1_N, B2_N, B3_N, E1_N, E2_N, &
                                               B1_S, B2_S, B3_S, E1_S, E2_S, &
                                               DThi_dTh_N, DThi_dPh_N, Thi_N, &
                                               DThi_dTh_S, DThi_dPh_S, Thi_S, &
                                               DThi_dTh_N_g, DThi_dPh_N_g, &
                                               DThi_dTh_S_g, DThi_dPh_S_g, &
                                               J_Th_N, J_Ph_N, E_Th_N, E_Ph_N, &
                                               J_Th_S, J_Ph_S, E_Th_S, E_Ph_S, &
                                               Thi_gnd_N, DThi_dr_N, &
                                               Thi_gnd_S, DThi_dr_S, &
                                               BetaN, Coeffs_B3_N, &
                                               BetaS, Coeffs_B3_S


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

!   sample period
    data_sp = 0.5
    plot_freq = anint (data_sp / (2*dt))

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

!  print*,Minval(Delta_T_FA),Minval(Delta_T_CF)

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



  allocate (D_out(Num_u1))
  allocate (Out_DL(Num_u1-1))
  allocate (h_mu_outer(Num_u3))

  D_out   = 0.0

  do ii = 1, Num_u3-1
    Delta_x = (xArr(Num_u1,ii+1) - xArr(Num_u1,ii))     ! Change in X coordinate along outer Boundary
    Delta_y = (yArr(Num_u1,ii+1) - yArr(Num_u1,ii))     ! Change in Y coordinate along outer Boundary
    Out_DL(ii) = sqrt(Delta_x**2 + Delta_y**2)
    D_out(ii+1) = D_out(ii) + Out_DL(ii)
  enddo
  Bnd_L = Sum(out_DL)                       ! Calculates length of outer boundary in Re

!  print*,Bnd_L

! CHECKED
  do ii = 1,Num_u3_2+1
    h_mu_outer(2*ii-1) = h_mu(Num_u1,2*ii-1)
  enddo

!  do ii = 1,Num_u3
!    write (11,'(i4,1X,f8.3)') ii,h_mu_outer(ii)
!  enddo

! Calc basis functions, V and pV/pr
  allocate (Ynm(K,Num_u1))
  allocate (Ynm_s(K,Num_u1))
  allocate (zeros(K))


!  CHECKED
!  print*,m_num, k, Num_u1, LMin_u1, LMax_u1
  call do_basis(m_num, k, Num_u1, LMin_u1, LMax_u1, del_Th, Ynm, Ynm_s, zeros)
!  print*,del_Th

!  do ii = 1, K
!    do jj = 1, Num_u1
!      write (11,'(2(i4,1X),2(E14.6,1X))') ii, jj, Ynm(ii,jj), Ynm_s(ii,jj)
!    enddo
!  enddo

!  do ii = 1, K
!    write (12,'(i4,1X,E14.6)') ii, zeros(ii)
!  enddo

! Ynm and Ynm_S are (k,N1) arrays where the Legendres are listed from HIGH to LOW Lat
! CoLat_N goes from LOW to HIGH Latitude i.e. HIGH to LOW CoLatitude
!
!  Set up scale factors for the Ionospheric Boundary Conditions for the Spherical Harmonic Expansion
!   Need to handle the alternating grid pattern on the ionopsheric boundary....

  allocate(Ynm_B3(K,Num_u1), Ynm_g(K,Num_u1), Ynm_B3_dr(K,Num_u1), Ynm_B3_dr2(K,Num_u1_2))
  allocate(AlphaB3(K,K))

  do kk = 1,K                    ! k = Number of Basis Functions in Expansion (from Sav file)
! B3 and B2 are on the same grid in u1
! South and North are the same
    do ii = 1,Num_u1
      Ynm_B3(kk,ii) = Ynm(kk,Num_u1+1-ii)            ! Reverse order of Basis Vectors due to U1 direction (0 = low boundary N1-1 = upper boundary)
      Ynm_B3_dr(kk,ii)=Ynm_S(kk,Num_u1+1-ii)
      Ynm_g(kk,ii)=Ynm(kk,Num_u1+1-ii)*(2.d0*zeros(kk)+1.d0)/(zeros(kk)+1.d0)
    enddo
    do ii = 1,Num_u1_2
      Ynm_B3_dr2(kk,ii) = Ynm_S(kk,Num_u1+1-(2*ii))  ! Reverse order and use 2nd value of Basis Vectors due to U1 direction (0 = low boundary N1-1 = upper boundary)
    enddo
  enddo

!  AlphaB3_2=MatMul(Ynm_B3_dr2,Transpose(Ynm_B3_dr2))
  AlphaB3  =MatMul(Ynm_B3_dr, transpose(Ynm_B3_dr))

!  do ii = 1, K
!    do jj = 1, K
!      write(11,'(2(i4,1X),E14.6)') ii, jj, AlphaB3 (ii, jj)
!    enddo
!  enddo

  count1 = 0
  time=0.d0

!  call get_time(tfrom,0,IDlog)

!   *******************************************************************************************************
!   Start Time iteration of Fields
!   ******************************************************************************************************

!  D_u1(:,:)=u1Arr(IR(:),:)-u1Arr(IL(:),:)!   Set up differences on Grid
!  D_u3(:,:)=u3Arr(:,JR(:))-u3Arr(:,JL(:))

  allocate (Tdr(Num_u3))

  Tdr = 0.0
!   Time Driver on Equatorial Slice (constant u1 surface)
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  h_width = 5.0

!  print*,D_out(Num_u3_2)

  do ii = 1, Num_u3, 2  ! CHECKED
    zz = (D_out(ii)-D_out(Num_u3_2))/h_width
    Tdr(ii) = exp(-zz**2)
  enddo

!  do ii = 1, Num_u3
!    write(11, '(i4,1X,2(E18.6E3,1X))') ii, Tdr(ii), D_out(ii)
!  enddo


! generate shift index arrays
  allocate (S1L(Num_u1), S1R(Num_u1), S3L (Num_u3), S3R(Num_u3)) 
  do ii = 1, Num_u1
    S1L(ii) = ii - 1
    S1R(ii) = ii + 1
  enddo
  S1R(Num_u1) = 1
  S1L(1) = Num_u1

  do ii = 1, Num_u3
    S3L(ii) = ii - 1
    S3R(ii) = ii + 1
  enddo
  S3R(Num_u3) = 1
  S3L(1) = Num_u3

  allocate (D_u1(Num_u1,Num_u3), D_u3(Num_u1,Num_u3))
  D_u1 = u1Arr(S1L(:),:) - u1Arr(S1R(:),:)
  D_u3 = u3Arr(:,S3L(:)) - u3Arr(:,S3R(:))

! ---------------------------------------------------------------------
!                       Start Iteration loop
! --------------------------------------------------------------------

  allocate (drv_B3(Num_u3))
  allocate (B1_con(Num_u1,Num_u3), B2_con(Num_u1,Num_u3), B3_con(Num_u1,Num_u3))
  allocate (B1_cov(Num_u1,Num_u3), B2_cov(Num_u1,Num_u3), B3_cov(Num_u1,Num_u3))
  allocate (E1_con(Num_u1,Num_u3), E2_con(Num_u1,Num_u3))
  allocate (E1_cov(Num_u1,Num_u3), E2_cov(Num_u1,Num_u3))

  allocate (Av_B1(Num_u1,Num_u3), Av_B3(Num_u1,Num_u3))

  allocate (B1_N(Num_u1), B2_N(Num_u1), B3_N(Num_u1_2))
  allocate (B1_S(Num_u1), B2_S(Num_u1), B3_S(Num_u1_2))
  allocate (E1_N(Num_u1), E2_N(Num_u1))
  allocate (E1_S(Num_u1), E2_S(Num_u1))

  allocate (Thi_N(Num_u1))
  allocate (Thi_S(Num_u1))
  allocate (DThi_dPh_N(Num_u1), DThi_dth_N(Num_u1))
  allocate (DThi_dPh_S(Num_u1), DThi_dth_S(Num_u1))
  allocate (j_Th_N(Num_u1), j_Ph_N(Num_u1))
  allocate (j_Th_S(Num_u1), j_Ph_S(Num_u1))
  allocate (E_Th_N(Num_u1), E_Ph_N(Num_u1))
  allocate (E_Th_S(Num_u1), E_Ph_S(Num_u1))
  allocate (Thi_gnd_N(Num_u1), DThi_dTh_N_g(Num_u1))
  allocate (Thi_gnd_S(Num_u1), DThi_dTh_S_g(Num_u1))

  allocate (BetaN(K), Coeffs_B3_N(K))
  allocate (BetaS(K), Coeffs_B3_S(K))
  allocate (IPIV(K))

  LDB = K
  IPIV = 0

  B1_N = 0.d0
  B2_N = 0.d0
  B3_N = 0.d0
  E1_N = 0.d0
  E2_N = 0.d0

  B1_S = 0.d0
  B2_S = 0.d0
  B3_S = 0.d0
  E1_S = 0.d0
  E2_S = 0.d0

  B1_con = 0.d0
  B2_con = 0.d0
  B3_con = 0.d0
  E1_con = 0.d0
  E2_con = 0.d0

  B1_cov = 0.d0
  B2_cov = 0.d0
  B3_cov = 0.d0
  E1_cov = 0.d0
  E2_cov = 0.d0

  do tt = 1,Nt-1                             ! Start of Time Iteration Loop
    time = 2.d0*Dt*tt                              ! Time Counter
    print*,tt
    bfa = 2.8e-9*sin(2*pi*time/30.0)
    drv_B3 = Tdr*bfa*Re_m*h_mu_outer

    B3_cov(Num_u1,1:Num_u3) = drv_B3
    B3_cov(Num_u1,1) = 0.d0

!   Conducting Northern Ionosphere
    E1_cov(1:Num_u1,1) = E1_N(1:Num_u1)
    E2_cov(1:Num_u1,1) = E2_N(1:Num_u1)

!   Conducting Southern Ionosphere
    E1_cov(1:Num_u1,Num_u3) = E1_S(1:Num_u1)
    E2_cov(1:Num_u1,Num_u3) = E2_S(1:Num_u1)

!   Inner L Shell Boundary
    E2_cov(1,1:Num_u3) = 0.d0
    B1_cov(1,1:Num_u3) = 0.d0

!   E field iteration
    E1_con = 2.d0 * dt * V2 / Jacb * (im * m_num * B3_cov - &
        (B2_cov(:,S3L(:)) - B2_cov(:,S3R(:))) / D_u3) + E1_con
    E2_con = 2.d0 * dt * V2 / Jacb * ((B1_cov(:,S3L(:)) - B1_cov(:,S3R(:))) / &
        D_u3 - B3_cov(S1L(:),:) - B3_cov(S1R(:),:)) + E2_con

!   Evolve Contravariant (Tangent) vector fields to Covariant vectors fields
    E1_cov = E1_con / g11_con
    E2_cov = E2_con * g22_cov

!   Northern Ionosphere
    E1_cov(1:Num_u1,1)= E1_N(1:Num_u1)
    E2_cov(1:Num_u1,1)= E2_N(1:Num_u1)      ! Conducting Northern Ionosphere using BC on Ionosphere

!   Southern Ionosphere
    E1_cov(1:Num_u1,Num_u3) = E1_S(1:Num_u1)
    E2_cov(1:Num_u1,Num_u3) = E2_S(1:Num_u1)      ! Conducting Southern Ionosphere using BC on Ionosphere

!   Inner L shell Boundary
    E2_cov(1,1:Num_u3)= 0.d0  ! Perfectly Reflecting
    E2_con(1,:) = 0.d0          ! inner L shell field line
    E2_con(:,1) = 0.d0          ! E2_con is populated as E2_CON(*,evens)
    E2_con(:,Num_u3) = 0.d0   ! other end

    B1_cov(1,1:Num_u3)= 0.d0
    E1_con(1,:) = 0.d0
    E1_con(:,1) = 0.d0
    E1_con(:,Num_u3) = 0.d0

!   B field iteration
    B1_con = 2.d0 * dt / Jacb * ((E2_cov(:,S3L(:)) - E2_cov(:,S3R(:))) / &
         D_u3) + B1_con
    B1_con = -2.d0 * dt / Jacb * ((E1_cov(:,S3L(:)) - E1_cov(:,S3R(:))) / &
         D_u3) + B2_con
    B3_con = 2.d0 * dt / Jacb * (im * m_num * E1_cov - (E2_cov(S1L(:),:) - &
        E2_cov(S1R(:),:)) / D_u1) + B3_con

!   Outer Corner Point Ionopshere and Outer L shell
    B3_con(Num_u1,1) = 4.d0 / 3.d0 * B3_con(Num_u1-2,1) - &
        1.d0 / 3.d0 * B3_con(Num_u1-4,1)
    B3_con(Num_u1,Num_u3) = 4.d0/3.d0*B3_con(Num_u1-2,Num_u3) - &
        1.d0 / 3.d0 * B3_con(Num_u1-4,Num_u3)

!   Evolve Contravariant Interior Grid B fields to Convariant vectors fields
    Av_B3 = (B3_con(S1R(:),S3R(:)) + B3_con(S1R(:),S3L(:)) + &
        B3_con(S1L(:),S3R(:)) + B3_con(S1L(:),S3L(:))) / 4.d0
    Av_B3(:,1) = 0.0d0
    Av_B3(:,Num_u3) = 0.0d0
    B1_cov = B1_con * g11_cov + Av_B3 * g13_cov
    B1_cov(:,1) = 0.d0
    B1_cov(:,Num_u3) = 0.d0                    ! Clean up after shifts due to odd number of points

!   Inner L shell Boundary
    B1_cov(1,:) = 0.d0

!   Evolve B2_cov
    B2_cov= B2_con * g22_cov
    B2_cov(:,1) = 0.0d0
    B2_cov(:,Num_u3) = 0.0d0   ! Clean up after shifts due to odd number of points

!   Evolve B3_cov
    Av_B1 = (B1_con(S1R(:),S3R(:)) + B1_con(S1R(:),S3L(:)) + &
        B1_con(S1L(:),S3R(:)) + B1_con(S1L(:),S3L(:))) / 4.d0
    Av_B1(:,1) = 0.0
    Av_B1(:,Num_u3) = 0.0d0                     ! Clean up after shifts due to odd number of points
    Av_B1(:,2) = 0.0d0
    Av_B1(:,Num_u3-1)=0.0d0
    B3_cov = Av_B1 * g13_cov + B3_con * g33_cov

!   Along Northern Ionospheric Boundary
    Av_B1 = (B1_con(S1R(:),S3L(:)) + B1_con(S1L(:),S3L(:))) / 2.0d0
    B3_cov(:,1) = Av_B1(:,1) * g13_cov(:,1) + B3_con(:,1) * g33_cov(:,1)

!   Along Southern Ionosphere Boundary
    Av_B1 = (B1_con(S1R(:),S3R(:)) + B1_con(S1L(:),S3R(:))) / 2.0d0
    B3_cov(:,Num_u3) = Av_B1(:,Num_u3) * g13_cov(:,Num_u3) + B3_con(:,Num_u3) * &
        g33_cov(:,Num_u3)

!   Along Inner L shell Boundary
    Av_B1 = (B1_con(S1L(:),S3R(:)) + B1_con(S1L(:),S3L(:))) / 2.0d0
    B3_cov(1,:) = Av_B1(1,:) * g13_cov(1,:) + B3_con(1,:) * g33_cov(1,:)

!   Along Outer L shell Boundary
    B3_cov(Num_u1,:) = drv_B3

! Ionosphere Boundary Condition
!       Using B3_con? in Ionosphere and perfectly conducting ground (B3_con = 0)
!       solve Laplaces equation in Neutral Atmopshere
!       Calculate B1_cov and B2_cov just below the Ionospheric current sheet.
!       Calculate jump in B1 and B2 and calculate Current in Sheet
!       Invert to find E1 and E2 values in sheet and use as BC in next
!       iteration.
!       Calculate Coefficients in expansion of Bz at the ionosphere (with
!       derivative radial component scale at r = Ri)
! At ionosphere boundary:
! Contravariant at ionos boundary for BB3 is radial
! Covariant at ionos boundary for BB1(colat) and BB2(azimuth) is in the sheet
!
!       FOR NORTHERN HEMISPHERE
!==============================================================================
!==============================================================================
!==============================================================================
!    B3_N1 = (B3_con(:,1)) * h_ra_n ! 1,3,5 etc. of full array
!   Interpolation
    do ii = 1,Num_u1_2
      B3_N(ii)= B3_con(2*ii,1) * h_ra_n(2*ii) ! Use only those points directly calculated
    enddo
! B3_N is a complex array
! Fit Br data to pV/pr basis set
    BetaN = matmul(B3_N,transpose(Ynm_B3_dr2))
    Coeffs_B3_N(:) = BetaN(:)
!    call zgesv(K, K, AlphaB3, K, IPIV, Coeffs_B3_N, K, INFO)
    DThi_dr_N = matmul(transpose(Ynm_B3_dr),Coeffs_B3_N) ! This should be B3_N at 2x the spatial resolution
    Thi_N = matmul(Ynm_B3,Coeffs_B3_N)          ! Potential Function

! Calc. pV/p_Theta
    do jj=2,Num_u1-1 
      DThi_dTh_N(jj) = (Thi_N(jj+1) - Thi_N(jj-1)) / (-2.d0 * del_Th)
    enddo

! Do the end points
    DThi_dTh_N(1) = (-3.0 * Thi_N(1) + 4.d0 * Thi_N(2) - Thi_N(3)) / &
        (-2.d0 * del_th)
    DThi_dTh_N(Num_u1) = (3.0 * Thi_N(Num_u1) - 4.d0 * Thi_N(Num_u1-1) + &
        Thi_N(Num_u1-2)) / (-2.d0 * del_th)
    DThi_dTh_N = 1.d0 / RI_s * DThi_dTh_N        ! B_Theta Atmos
    ! B_phi Atmos, CoLat_N runs from INNER L shell to OUTER L shell
    DThi_dPh_N = CMPLX(0.0,1.d0) * m_num * Thi_N / (RI_s * sin(Colat_N))

! Interpolate B1 and B2 to ALL points just above the ionosphere from model
! solution to calculate currents in Ionosphere
    B1_N = B1_cov(:,2) / h_th_n
    B2_N = B2_cov(:,2) / h_ph_n
    B1_N = (B1_N(S1R(:)) + B1_N(S1L(:))) / 2.0 + B1_N ! Linear Interpolation
    B2_N = (B2_N(S1R(:)) + B2_N(S1L(:))) / 2.0 + B2_N ! Linear Interpolation

! Try second order Taylor expansion for end points
    h = CoLat_N(Num_u1 - 1) - CoLat_N(Num_u1 - 3) ! step size along ionosphere is uniform on B1_N grid
    FirD = ( B1_N(Num_u1 - 5) - 4.0 * B1_N(Num_u1 - 3) + &
        3.0 * B1_N(Num_u1 - 1)) / (2.0 * h) ! First Derivative Backwards difference O(h^2)
    SecD = (-B1_N(Num_u1 - 7) + 4.0 * B1_N(Num_u1 - 5) - &
        5.0 * B1_N(Num_u1 - 3) + 2.0 *  B1_N(Num_u1 - 1)) / (h ** 2)       ! Second Derivative Backwards difference O(h^2)
    h1 = CoLat_N(Num_u1-2)-CoLat_N(Num_u1-1) ! Step size on B1_n_interp grid
    B1_N(Num_u1) = B1_N(Num_u1 - 1) + h1 * FirD + ((h1 ** 2) / 2.0) * SecD

    h = CoLat_N(4)-CoLat_N(2) !; step size along ionosphere is uniform
    FirD = (-B2_N(6) + 4.0 * B2_N(4) - 3.0 * B2_N(2)) / (2.0 * h)       ! First Derivative Forwards difference O(h^2)
    SecD = (-B2_N(8) + 4.0 * B2_N(6) - 5.0 * B2_N(4) + 2.0 * B2_N(2)) / (h ** 2)   ! Second Derivative Forwards difference O(h^2)
    h1 = CoLat_N(1) - CoLat_N(2)
    B2_N(1)= B2_N(2) + h1 * FirD + ((h1 ** 2) / 2.0) * SecD

! Currents in the Ionosphere
    J_Th_N = -(B2_N - DThi_dPh_N) / mu0_s
    J_Ph_N =  (B1_N - DThi_dTh_N) / mu0_s

! Calculate electric fields in Ionosphere from discontinuity in B's for
! next time step...
    E_Th_N = (INVSS_N(:,1,1) * J_Th_N + INVSS_N(:,2,1) * J_Ph_N)
    E_Ph_N = (INVSS_N(:,1,2) * J_Th_N + INVSS_N(:,2,2) * J_Ph_N)

    E1_N = E_Th_N * h_th_n * odds
    E2_N = E_Ph_N * h_ph_n * evens

!   Now do the Ground fields
    Thi_gnd_N = matmul(Ynm_g, Coeffs_B3_N)                            ! Potential Function
    do jj = 2,Num_u1-1 
      DThi_dTh_N_g(jj) = (Thi_gnd_N(jj + 1) - Thi_gnd_N(jj - 1)) / (-2.d0 * del_Th)
    enddo

! Do the end points
    DThi_dTh_N_g(1) = (-3.0 * Thi_gnd_N(1) + 4.d0 * Thi_gnd_N(2) - Thi_gnd_N(3))/(-2.d0 * del_th)
    DThi_dTh_N_g(Num_u1) = (3.0 * Thi_gnd_N(Num_u1) - 4.d0 * Thi_gnd_N(Num_u1 - 1) + Thi_gnd_N(Num_u1 - 2)) / (-2.d0 * del_th)
    DThi_dTh_N_g = 1.d0/RI_s*DThi_dTh_N_g ! B_Theta Ground
    DThi_dPh_N_g = Cmplx(0.0,1.d0) * m_num * Thi_gnd_N / (Re_s * sin(Colat_N)) !B_phi   Ground    

!------------------------------------------------------------------------------
!   SOUTHERN HEMISPHERE
!-----------------------------------------------------------------------------
!    B3_S1 = (B3_con(*,Num_u3))*h_ra_s ! 1,3,5 etc. of full array
!    B3_S1 = (Shift(B3_S1,1) + Shift(B3_S1,-1))/2.0 + B3_S1
!   Linear Interpolation
    do ii = 1,Num_u1_2
      B3_S(ii) = B3_con(2*ii,Num_u3) * h_ra_s(2*ii)
    enddo

!   Fit Br data to pV/pr basis set
    BetaS = matmul(B3_S,transpose(Ynm_B3_dr2))
    Coeffs_B3_S(:) = BetaS(:)
!    call zgesv(K, K, AlphaB3, K, IPIV, Coeffs_B3_S, K, INFO)
    DThi_dr_S = matmul(Ynm_B3_dr,Coeffs_B3_S)             ! This should be B3_S at twice the spatial resolution
    Thi_S = matmul(Ynm_B3,Coeffs_B3_S)                    ! Potential Function

!   Calc. pV/p_Theta
    do jj=2,Num_u1-1 
      DThi_dTh_S(jj) = (Thi_S(jj+1) - Thi_S(jj-1)) / (2.d0 * del_Th)
    enddo

! Do the end points
    DThi_dTh_S(1) = (-3.d0 * Thi_S(1) + 4.d0 * Thi_S(2) - Thi_S(3)) / &
        (2.d0 * del_th)
    DThi_dTh_S(Num_u1-1) = (3.d0 * Thi_S(Num_u1-1) - 4.d0 * Thi_S(Num_u1-2) + &
        Thi_S(Num_u1-3)) / (2.d0 * del_Th)
    DThi_dTh_S = 1.d0 / RI_s * DThi_dTh_S     ! B_Theta Atmos
    DThi_dPh_S = Cmplx(0.0,1.d0) * m_num * Thi_S / (RI_s * sin(Colat_S))  ! B_phi Atmos

!   Extrapolate  B1 and B2 to ALL points just above ionosphere from model
!   solution to calculate currents in Ionosphere
    B1_S= B1_cov(:,Num_u3-1) / h_th_s
    B2_S= B2_cov(:,Num_u3-1) / h_ph_s
    B1_S = (B1_S(S1R(:)) + B1_S(S1L(:))) / 2.0 + B1_S ! Linear Interpolation
    B2_S = (B2_S(S1R(:)) + B2_S(S1R(:))) / 2.0 + B2_S ! Linear Interpolation

!   Try second order Taylor expansion for end points
    h = CoLat_S(Num_u1-1) - CoLat_S(Num_u1-3) ! step size along ionosphere is uniform on B1_N grid
!   First Derivative Backwards difference O(h^2)
    FirD =( B1_S(Num_u1-5) - 4.0 * B1_S(Num_u1-3) + 3.0 * B1_S(Num_u1-1)) / &
        (2.0 * h)
!   Second Derivative Backwards difference O(h^2)
    SecD =(-B1_S(Num_u1-7) + 4.0 * B1_S(Num_u1-5) - 5.0 * B1_S(Num_u1-3) + &
        2.0 * B1_S(Num_u1-2)) / (h ** 2)
    h1 = CoLat_S(Num_u1-1) - CoLat_S(Num_u1) ! Step size on B1_n_interp grid
    B1_S(Num_u1-1) = B1_S(Num_u1-1) + h1 * FirD + ((h1 ** 2) / 2.0) * SecD

    h = CoLat_S(4) - CoLat_S(2) ! step size along ionosphere is uniform
!   First Derivative Forwards difference O(h^2)
    FirD = (-B2_S(6) + 4.0 * B2_S(4) - 3.0 * B2_S(2)) / (2.0 * h)
!   Second Derivative Forwards difference O(h^2)
    SecD = (-B2_S(8) + 4.0 * B2_S(6) - 5.0 * B2_S(4) + 2.0 * B2_S(2)) / (h ** 2)
    h1 = CoLat_S(1) - CoLat_S(2)
    B2_S(1) = B2_S(2) + h1 * FirD + ((h1 ** 2) / 2.0) * SecD

!       Currents in the Ionosphere
    J_Th_S = -(B2_S - DThi_dPh_S) / mu0_s
    J_Ph_S =  (B1_S - DThi_dTh_S) / mu0_s

!       Calculate electric fields in Ionosphere from discontinuity in B's for
!       next time step...
    E_Th_S = (INVSS_S(:,1,1) * J_Th_S + INVSS_S(:,2,1) * J_Ph_S)
    E_Ph_S = (INVSS_S(:,1,2) * J_Th_S + INVSS_S(:,2,2) * J_Ph_S)
    E1_S = E_Th_S * h_th_s * odds
    E2_S = E_Ph_S * h_ph_s * evens

! Now do the Ground fields
    Thi_gnd_S = matmul(Ynm_g, Coeffs_B3_S)  ! Potential Function
    do jj=2,Num_u1-1 
      DThi_dTh_S_g(jj) = (Thi_gnd_S(jj+1) - Thi_gnd_S(jj-1)) / (-2.d0 * del_Th)
    enddo

! Do the end points
    DThi_dTh_S_g(1) = (-3.0 * Thi_gnd_S(1) + 4.d0 * Thi_gnd_S(2) - Thi_gnd_S(3)) / &
        (-2.d0 * del_th)
    DThi_dTh_S_g(Num_u1) = (3.0 * Thi_gnd_S(Num_u1) - 4.d0 * Thi_gnd_S(Num_u1-1) + &
        Thi_gnd_S(Num_u1-2)) / (-2.d0 * del_th)
    DThi_dTh_S_g = 1.d0 / RI_s * DThi_dTh_S_g ! B_Theta Ground
    DThi_dPh_S_g = Cmplx(0.0, 1.d0) * m_num * Thi_gnd_S / (Re_s * sin(Colat_S)) ! B_phi   Ground
  enddo  ! ============== End of Time Loop ====================

