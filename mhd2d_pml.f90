!   Simulation of ULF wave propagation through the magnetosphere
!   with realistic ionospheric boundaries
!   Version 2.0     5 Aug 2006
!
!   C.L. Waters, M.D. Sciffer
!   Centre for Space Physics
!   University of Newcastle
!   New South Wales, Australia
!
! Mods:
! Changed the basis function routine to scaleup by 4 [Mar 2013]
! Added PML [Apr 2013]
!
! ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
!
  program mhd2d_pml

    use mhd2d_constants
    use mhd2d_va_cond
    use mhd2d_grid

    implicit none

    integer ::  ii, jj, kk
    integer ::  Num_u1, Num_u1_2, Num_u3, Num_u3_2, K
    integer(8) :: Nt, tt, tcnt
    integer :: plot_freq
    real(DBL) :: data_sp

    real(DBL) :: LMin_u1, LMax_u1, LMin_u3, LMax_u3, m_num
    real(DBL) :: dt, del_Th_0_u1, del_Th, drv_amp, mxDrv
    real(DBL) :: strt_ramp, t_ramp, drv_norm

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

! PML vars
    integer, allocatable, dimension(:) :: idx_tmp
    integer :: iL_src, PML_Loc
    real(DBL) :: PVaPML, PML_Pow, PMLL, L_src, scale_PML, Nu0, NuPML
    real(DBL), allocatable, dimension(:,:) :: NE1, NE2, NB1, NB2, NB3, &
        NuE1, NuE12, NuE2, NuE22, NuB1, NuB12, NuB2, NuB22, NuB3, NuB32

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
    integer :: use_va_file

! Courant condition arrays
    real(DBL), allocatable, dimension(:,:) :: Delta_T_CF, Delta_T_FA, &
                                              Delta_L_CF, Delta_L_FA, Eperp
    real(DBL), allocatable, dimension(:) :: D_out, FLR_Per
    real(DBL) :: Av_Va, Delta_x, Delta_y, Length, Courant_cond, &
                 Freq

! Basis set arrays
    real(DBL), allocatable, dimension(:,:) :: Ynm, Ynm_s, Ynm_B3, Ynm_g, &
                                              Ynm_B3_dr2, Ynm_B3_dr, &
                                              AlphaB3, inv_arr, b3coef_mult
    real(DBL), allocatable, dimension(:) :: zeros, WORK

!   loop variables
    real(DBL), allocatable, dimension(:,:) :: emult, bmult
    real(DBL), allocatable, dimension(:) :: Tdr, drv_B3
    real(DBL) :: time, hwidth, zz, bfa
    integer :: count1

!   shift arrays
    integer, allocatable, dimension(:) :: S1L, S1R, S3L, S3R, pS1L, pS1R

!   loop arrays
    integer ::LDB,  INFO
    real(DBL) :: h, h1
    complex(DBL) :: FirD, SecD
    integer, allocatable, dimension(:) :: IPIV
    complex(DBL), allocatable, dimension(:,:) :: B1_con, B2_con, B3_con, &
                                                 B1_cov, B2_cov, B3_cov, &
                                                 E1_con, E2_con, &
                                                 E1_cov, E2_cov, &
                                                 Av_B1, Av_B3, &
                                         b_nu_a, b_ph_a, b_mu_a, &
                                         e_nu_a, e_ph_a
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
                                               Coeffs_B3_N, Coeffs_B3_S
    character(len=100) :: out_dir, mk_out_dir, out_file
    real :: t1, t2, t3

! - - - - - - - - - start of code - - - - - - - - - - -

! Azimuthal wave number, m
    m_num = 2.0

! Number of iterations
    Nt = 400000
    use_va_file = 0        ! 1=constant ionos conductance
! dt time step
    dt = 0.00005

! Output data dir
    out_dir='/data/mhd2d/output'

    if (use_va_file.eq.1) then
! Parameters for run - to come from file
      print*,'Reading Grid parameters and Va data file ...'
!      iriva_file='/home/mdw/mhd2d/iri_va_sig.mhd'
      iriva_file='/home/mdw/mhd2d/iri_va_sig_pml.mhd'
      open(unit=10,file=iriva_file,status='old',action='read')
      read(10,*) yr,doy,ut_hr,geolon,Num_u1,Num_u3
      read(10,*) f107,f107a,ap
      read(10,*) LMin_u1, LMax_u1, LMin_u3,LMax_u3
    else
      Num_u1 = 210
      Num_u3 = 151
      LMin_u1 = 1.2d0
      LMax_u1 = 15.0d0
      LMin_u3 = RI_s
      LMax_u3 = 95.0d0
    endif

!   sample period
    data_sp = 1.0
    plot_freq = anint(data_sp/(2.0*dt))

! driver parameters for sharp pulse - CLW
    drv_amp = 1.0d-8       ! amplitude of sine
    hwidth = 5.0           ! inc to spread excitation further along field line
    t_ramp = 30.0
    drv_norm = 35.0

! PML parameters and driver (b_mu) field line
    L_src = 9.8            ! Re of driver
    PVaPML = 2.0           ! Power law on Va in PML
    PML_Pow = 2.0          ! Power law on loss terms in PML
    PMLL = 10.0            ! start L of PML

! For additional density near ionosphere
    PLaw = 4.0
    nO = 1.0d7             ! O2 density (kg/Re^3)
    nO_Ht = 250.0d3/rE_m   ! scale height of oxygen
  
! - - - - - end parameters - - - - -

! Check grid size and numbers
    Num_u1_2 = int(Num_u1/2.0)
    Num_u1 = int(2*Num_u1_2)

    Num_u3_2 = int(Num_u3/2.0)
    Num_u3 = int(2*Num_u3_2)+1

! spherical harmonic soln order for atmosphere
    K = int(Num_u1/4)

! setup solution grid arrays
    allocate (LVArr(Num_u1), CoLat_N(Num_u1), CoLat_S(Num_u1) )
    allocate (u1Arr(Num_u1,Num_u3), u3Arr(Num_u1,Num_u3) )
    allocate (RArr(Num_u1,Num_u3), ThArr(Num_u1,Num_u3) )
    allocate (xArr(Num_u1,Num_u3), yArr(Num_u1,Num_u3) )
! Contravariant geometric factors
    allocate (g11_con(Num_u1,Num_u3), g13_con(Num_u1,Num_u3) )
    allocate (g22_con(Num_u1,Num_u3), g33_con(Num_u1,Num_u3) )
! Covariant geomateric factors
    allocate (g11_cov(Num_u1,Num_u3), g13_cov(Num_u1,Num_u3) )
    allocate (g22_cov(Num_u1,Num_u3), g33_cov(Num_u1,Num_u3) )
! Jacobian to allow conversion CON<=>COV
    allocate (Jacb(Num_u1,Num_u3) )
! Spherical coord geomatric factors for ionos boundary
    allocate (h_nu(Num_u1,Num_u3) )
    allocate (h_ph(Num_u1,Num_u3) )
    allocate (h_mu(Num_u1,Num_u3) )

! Call grid generation routine
    print*,'Generating grid...'
! Inputs:
!   Num_u1, Num_u1_2
!   Num_u3, Num_u3_2
!   LMin_u1, LMax_u1
!   LMin_u3, LMax_u3
! Outputs:
!   del_Th, del_Th_0_u1               (real)
!   u1Arr, u3Arr                      (real[Num_u1,Num_u3])
!   RArr, ThArr, LVArr, xArr, yArr
!   g11_con, g13_con, g22_con, g33_con
!   g11_cov, g13_cov, g22_cov, g33_cov, Jacb
!   CoLat_N, CoLat_S                  (real[Num_u1])
!   h_nu, h_ph, h_mu
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
!    print*,'Stopped after grid generation.'
!    stop
! # # #

! Find index of LVArr of driver field line
    iL_src = 0
    allocate(idx_tmp(Num_u1))
    idx_tmp = 0
    jj = 1
    do ii=1,Num_u1
      if (LVArr(ii).ge.L_src) then
        idx_tmp(jj) = ii
        jj=jj+1
      endif
    enddo
    iL_src = nint(idx_tmp(1)/2.0)*2           ! idx of PML
    print*,'iL_src=',iL_src

! Find index of start of PML
    idx_tmp = 0
    jj = 1
    do ii=1,Num_u1
      if (LVArr(ii).ge.PMLL) then
        idx_tmp(jj) = ii
        jj=jj+1
      endif
    enddo
    PML_Loc = nint(idx_tmp(1)/2.0)*2           ! idx of PML
    deallocate(idx_tmp)
    print*,'PML_Loc=',PML_Loc

! Plasma mass density
    allocate (rho_a(Num_U1,Num_u3) )
! Va and conductance arrays
    allocate(Va_a(Num_u1,Num_u3))    
    allocate (Sd_N(PML_Loc))
    allocate (Sp_N(PML_Loc))
    allocate (Sh_N(PML_Loc))
    allocate (Sd_S(PML_Loc))
    allocate (Sp_S(PML_Loc))
    allocate (Sh_S(PML_Loc))

! If we are using IRI generated Va and conductance data file
    if (use_va_file.eq.1) then
! Read in Va data
      print*,'Reading Va and ionosphere data...'
      do ii=1,Num_u1
        do jj=1,Num_u3
          read(10,*) i1, i2, lat, va_sp
          Va_a(i1+1,i2+1) = va_sp
        enddo
      enddo
! Read in conductance data
      do ii=1,Num_u1
        read(10,*) lat, s0n,s1n,s2n,s0s,s1s,s2s
        if (ii.le.PML_Loc) then
          Sd_N(ii)=s0n
          Sp_N(ii)=s1n
          Sh_N(ii)=s2n
          Sd_S(ii)=s0s
          Sp_S(ii)=s1s
          Sh_S(ii)=s2s
        endif
      enddo
      close(unit=10)

      print*,'Calculating Va...'
      do ii=1,Num_u1                   ! chooses which field line
        xd = LVArr(ii)
        lat = 0.0d0                    ! set to equator
        call Get_B0(Lat,xd,bmg)        ! Given xd and Latitude, calc B_0 at equator
        va_sp=Va_a(ii,Num_u3_2)        ! use data from input file
        nH = bmg**2/(mu0*va_sp**2)   ! Calc Hydrogen plasma number density at the equator

        do jj=1,Num_u3                 ! Loop over points along field line
          Lat=pi/2.0-ThArr(ii,jj)
          call Get_B0(Lat,xd,bmg)
!          rho_a(ii,jj)=bmg**2/(mu0*Va_a(ii,jj)**2)
          rho_a(ii,jj) = nH*(xd/RArr(ii,jj))**PLaw + nO/rE_m**3* &
            exp(-(Rarr(ii,jj)-Re_s)/nO_Ht)
          Va_a(ii,jj) = bmg/sqrt(mu0*rho_a(ii,jj))/Re_m
!          Va_a(ii,jj) = Va_a(ii,jj)/rE_m      ! modify Va
        enddo                                   ! Num_u3 loop

! call eigen_solver to check FLR frequencies
!   call Find_flr(xd,N3,rho_a(ii,:),har)
!   flr_a(ii,0:5)=har(0:5)             ! Fundamental FLR
      enddo                            ! Num_u1 loop
    else
! If we are generating CLW Va profile and constant conductance case
      print*,'Calculating Va...'
      do ii=1,Num_u1
        xd = LVArr(ii)
        call clw_va(xd,va_sp)          ! Given xd, calc Va at equator
        Lat = 0.0d0
        call Get_B0(Lat,xd,bmg)        ! Given xd and Latitude, calc B_0 at equator
        nH   = bmg**2/(mu0*va_sp**2)   ! Calc Hydrogen plasma number density at the equator
        nO = 1.0d7
        nO_Ht = 250.0d3/Re_m
        do jj=1,Num_u3                 ! Loop over points along field line
          Lat=pi/2.0-ThArr(ii,jj)
          call Get_B0(Lat,xd,bmg)
!          rho_a(ii,jj)=bmg**2/(mu0*Va_a(ii,jj)**2)
          rho_a(ii,jj) = nH*(xd/RArr(ii,jj))**PLaw + nO/rE_m**3* &
            exp(-(Rarr(ii,jj)-Re_s)/nO_Ht)
          Va_a(ii,jj) = bmg/sqrt(mu0*rho_a(ii,jj))/Re_m      ! modify Va
        enddo                                                ! Num_u3 loop
      enddo                                                  ! Num_u1 loop
    endif                                                    ! if use_va_file

! Get memory for scale factors and conductance matrices
    allocate (h_ph_n(PML_Loc))
    allocate (h_th_n(PML_Loc))
    allocate (h_ra_n(PML_Loc))
    allocate (h_ph_s(PML_Loc))
    allocate (h_th_s(PML_Loc))
    allocate (h_ra_s(PML_Loc))

    allocate (INVSS_N(PML_Loc,2,2))
    allocate (INVSS_S(PML_Loc,2,2))
!
    print*,'Calculating ionosphere conductance...'
! Inputs
!   PML_Loc (idx of outer soln space field line)
!   CoLat_N, CoLat_S               (real[Num_u1])
!   use_va_file (0 or 1)
! Output
!   Sd_N, Sp_N, Sh_N               (real[PML_Loc])
!   Sd_S, Sp_S, Sh_S
!   INVSS_N, INVSS_S               (real[PML_Loc,2,2])]
!   h_ph_n, h_th_n, h_ra,n         (real[PML_Loc])
!   h_ph_s, h_th_s, h_ra,s
    call Calc_cond_ionos(PML_Loc, CoLat_N, CoLat_S, &
                 Sd_N,Sp_N,Sh_N,INVSS_N, &
                 Sd_S,Sp_S,Sh_S,INVSS_S, &
                 h_ph_n,h_th_n,h_ra_n, &
                 h_ph_s,h_th_s,h_ra_s,use_va_file)

! =========================================================
! TEST - output to test calc_cond_ionos
!do ii=1,Num_u1
!  write(10,*) sd_n(ii),' ',sp_n(ii),' ',sh_n(ii),' ', &
!              sd_s(ii),' ',sp_s(ii),' ',sh_s(ii)
!  write(11,*) h_ph_n(ii),' ',h_th_n(ii),' ',h_ra_n(ii),' ',&
!              h_ph_s(ii),' ',h_th_s(ii),' ',h_ra_s(ii)
!  do jj=1,Num_u3
!    write(12,*) Va_a(ii,jj)
!  enddo
!enddo
!print*,'Stopped after calc_conductance'
!stop
!=========================================================
!
! Discretisation Scheme
! The number of full grid cells is (Num_u1_2, Num_u3_2)
! In this programing section ii and jj indentifies the field values in the (ii,jj) cell
! and within the cell the location of parameters (such as Alfven speed etc ) associated
! with a particular position in space are given by ....
!
!       (2*ii,   2*jj-1)      location of E1
!       (2*ii-1, 2*jj-1)      location of E2
!       (2*ii-1, 2*jj)        location of B1
!       (2*ii,   2*jj)        location of B2
!       (2*ii,   2*jj-1)      location of B3   ( located at the same position in space as E1)
!
! Set up array for output plots of fields
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
!write(13,*) ii,' ',jj,' ',b1x_ord(ii,jj),' ',b1y_ord(ii,jj),' ', &
!                          b2x_ord(ii,jj),' ',b2y_ord(ii,jj),' ', &
!                          b3x_ord(ii,jj),' ',b3y_ord(ii,jj),' ', &
!                          e1x_ord(ii,jj),' ',e1y_ord(ii,jj),' ', &
!                          e2x_ord(ii,jj),' ',e2y_ord(ii,jj)
    enddo
  enddo 
!   Addition Half cell in u1 direction
!jj = num_u3_2+1
  do ii = 1,Num_u1_2
    E1x_ord(ii,Num_u3_2+1) = Xarr(2*ii  ,Num_u3)
    E2x_ord(ii,Num_u3_2+1) = Xarr(2*ii-1,Num_u3)
    B3x_ord(ii,Num_u3_2+1) = Xarr(2*ii  ,Num_u3)
    E1y_ord(ii,Num_u3_2+1) = Yarr(2*ii  ,Num_u3)
    E2y_ord(ii,Num_u3_2+1) = Yarr(2*ii-1,Num_u3)
    B3y_ord(ii,Num_u3_2+1) = Yarr(2*ii  ,Num_u3)
!write(14,*) ii,' ',jj,' ',b3x_ord(ii,jj),' ',b3y_ord(ii,jj),' ', &
!                          e1x_ord(ii,jj),' ',e1y_ord(ii,jj),' ', &
!                          e2x_ord(ii,jj),' ',e2y_ord(ii,jj)
  enddo 

!print*,'stopped at plot grid'
!stop

! Calc Va^2 used in time iteration which is time invariant over the full grid
  allocate (V2(Num_u1, Num_u3))
  allocate (Eperp(Num_u1,Num_u3))
  do ii=1,Num_u1
    do jj=1,Num_u3
      Eperp(ii,jj)  = eps0_s*(1.0 + c_s**2/Va_a(ii,jj)**2)   ! in SI units
      V2(ii,jj) = 1.0/(mu0_s*Eperp(ii,jj))  ! Va^2 in SI
    enddo
  enddo

! Calculate Info about run
  allocate (Delta_L_FA(Num_u1,Num_u3))
  allocate (Delta_T_FA(Num_u1,Num_u3))
  allocate(FLR_Per(Num_u1))

  do ii = 1, Num_u1
    do jj = 2, Num_u3-1
      Delta_x = (xArr(ii,jj+1) - xArr(ii,jj-1))         !     Change in X coordinate along field line
      Delta_y = (yArr(ii,jj+1) - yArr(ii,jj-1))         !     Change in Y coordinate along field line
      Delta_L_FA(ii,jj)= sqrt(Delta_x**2 + Delta_y**2)  !     Line segment Length (in Re)
      Av_Va =max(Va_a(ii,jj+1),Va_a(ii,jj-1))
      Delta_T_FA(ii,jj)= Delta_L_FA(ii,jj)/Av_Va        !    Transit time across Delta_L segment
      FLR_Per(ii) = FLR_Per(ii) + Delta_T_FA(ii,jj)     !    Transit time along field line
    enddo
  enddo

  Delta_L_FA(:,Num_u3) = Maxval(Delta_L_FA)
  Delta_T_FA(:,Num_u3) = Maxval(Delta_T_FA)
  Delta_L_FA(:,1) = Maxval(Delta_L_FA)
  Delta_T_FA(:,1) = Maxval(Delta_T_FA)

  allocate (Delta_L_CF(Num_u1,Num_u3))
  allocate (Delta_T_CF(Num_u1,Num_u3))

  do ii = 2, Num_u1-1
    do jj = 1, Num_u3
!   Across field line
      Delta_x = (xArr(ii+1,jj) - xArr(ii-1,jj))  !  Change in X coordinate along field line
      Delta_y = (yArr(ii+1,jj) - yArr(ii-1,jj))  !  Change in Y coordinate along field line
      Length = sqrt(Delta_x**2 + Delta_y**2)     !  Line segment Length (in Re)
      Av_Va =Max(Va_a(ii+1,jj),Va_a(ii-1,jj))    !  Average Alfven Speed across Delta_L segment (in m/s)
      Delta_L_CF(ii,jj)= Length                  !  Transit time across Delta_L segment
      Delta_T_CF(ii,jj)= Delta_L_CF(ii,jj)/Av_Va !  Transit time across Delta_L segment
    enddo
  enddo

  Delta_L_CF(Num_u1,:) = Maxval(Delta_L_CF)
  Delta_T_CF(Num_u1,:) = Maxval(Delta_T_CF)
  Delta_L_CF(1,:) = Maxval(Delta_L_CF)
  Delta_T_CF(1,:) = Maxval(Delta_T_CF)

  Courant_cond = min(Minval(Delta_T_FA),Minval(Delta_T_CF))/2
  print*,'Courant cond = ',courant_cond

  if (dt.GT.Courant_cond) then
    print*,"WARNING: Courant condition problem. Try reducing 'dt'"
    stop
  endif

  deallocate(Delta_T_CF)
  deallocate(Delta_L_CF)
  deallocate(FLR_Per)
  deallocate(Delta_T_FA)
  deallocate(Delta_L_FA)

! reduce Va in PML
  do ii=PML_Loc+1,Num_u1
    do jj=1,Num_u3
      Va_a(ii,jj) = Va_a(PML_Loc,jj)*(u1Arr(ii,jj)/u1Arr(PML_Loc,jj))**PVaPML
    enddo
  enddo
! Recalc V2 after PML modification
  do ii=1,Num_u1
    do jj=1,Num_u3
      Eperp(ii,jj)  = eps0_s*(1.0 + c_s**2/Va_a(ii,jj)**2)   ! in SI units
      V2(ii,jj) = 1.0/(mu0_s*Eperp(ii,jj))  ! Va^2 in SI
    enddo
  enddo

! Width of PML in u1 units
  scale_PML = (u1Arr(Num_u1,1) - u1Arr(PML_Loc,1))

  allocate(NE1(Num_u1,Num_u3))
  allocate(NE2(Num_u1,Num_u3))
  allocate(NB1(Num_u1,Num_u3))
  allocate(NB2(Num_u1,Num_u3))
  allocate(NB3(Num_u1,Num_u3))

  allocate(NuE1(Num_u1,Num_u3))
  allocate(NuE12(Num_u1,Num_u3))
  allocate(NuE2(Num_u1,Num_u3))
  allocate(NuE22(Num_u1,Num_u3))
  allocate(NuB1(Num_u1,Num_u3))
  allocate(NuB12(Num_u1,Num_u3))
  allocate(NuB2(Num_u1,Num_u3))
  allocate(NuB22(Num_u1,Num_u3))
  allocate(NuB3(Num_u1,Num_u3))
  allocate(NuB32(Num_u1,Num_u3))

! Fill PML arrays
  do jj=1,Num_u3
    Nu0 = sqrt(1.0d0/(mu0*Eperp(PML_Loc,jj)))
    NuPML = Nu0/(1.0d2)**PML_Pow
    NE1(PML_Loc+1:Num_u1,jj) = NuPML*((u1Arr(PML_Loc+1:Num_u1,jj)-&
                      u1Arr(PML_Loc,jj))/scale_PML)**PML_Pow
    NE2(PML_Loc+1:Num_u1,jj) = NuPML*((u1Arr(PML_Loc+1:Num_u1,jj)-&
                      u1Arr(PML_Loc,jj))/scale_PML)**PML_Pow
    NB1(PML_Loc+1:Num_u1,jj) = mu0_s/Eperp(PML_Loc,jj)*NuPML* &
         ((u1Arr(PML_Loc+1:Num_u1,jj)-u1Arr(PML_Loc,jj))/scale_PML)**PML_Pow
    NB2(PML_Loc+1:Num_u1,jj) = mu0_s/Eperp(PML_Loc,jj)*NuPML* &
         ((u1Arr(PML_Loc+1:Num_u1,jj)-u1Arr(PML_Loc,jj))/scale_PML)**PML_Pow
    NB3(PML_Loc+1:Num_u1,jj) = mu0_s/Eperp(PML_Loc,jj)*NuPML* &
         ((u1Arr(PML_Loc+1:Num_u1,jj)-u1Arr(PML_Loc,jj))/scale_PML)**PML_Pow
  enddo
  NuE1  = exp(-NE1*dt*1.d0)
  NuE12 = exp(-NE1*dt*2.d0)
  NuE2  = exp(-NE2*dt*1.d0)
  NuE22 = exp(-NE2*dt*2.d0)

  NuB1  = exp(-NB1*dt*1.d0)
  NuB12 = exp(-NB1*dt*2.d0)
  NuB2  = exp(-NB2*dt*1.d0)
  NuB22 = exp(-NB2*dt*2.d0)
  NuB3  = exp(-NB3*dt*1.d0)
  NuB32 = exp(-NB3*dt*2.d0)

! Calc D_out for spatial function variable on driver field line
  allocate (D_out(Num_u3))
  allocate (h_mu_outer(Num_u3))

  D_out   = 0.0
  do ii = 1, Num_u3-1
    Delta_x = (xArr(Num_u1,ii+1) - xArr(Num_u1,ii))     ! Change in X coordinate along outer Boundary
    Delta_y = (yArr(Num_u1,ii+1) - yArr(Num_u1,ii))     ! Change in Y coordinate along outer Boundary
    D_out(ii+1) = D_out(ii) + sqrt(Delta_x**2+Delta_y**2)
  enddo

! scale factors for compressional dB on outer L shell
  do ii = 1,Num_u3_2+1
    h_mu_outer(2*ii-1) = h_mu(iL_src,2*ii-1)
  enddo

  print*,'Generate atmosphere soln basis functions...'
! Calc basis functions, V and pV/pr
  allocate (Ynm(K,PML_Loc))
  allocate (Ynm_s(K,PML_Loc))
  allocate (zeros(K))
  call do_basis(m_num, k, PML_Loc, LMin_u1, LMax_u1, del_Th, Ynm, Ynm_s, zeros)
!
!============================================================================
! TEST atmosphere basis function routine
!print*,'Stopped at basis function routine'
!stop

!
! Ynm and Ynm_S are real(k,PML_Loc) arrays where the Legendre basis functions
! are stored from high to low Latitude
! CoLat_N goes from low to high latitude i.e. high to low co_latitude
!
! Calc scale factors for ionospheric boundary from the Spherical Harmonic expansion
!   Need to handle the alternating grid pattern on the ionopsheric boundary.
!   This is done using Ynm_B3_dr2.
  allocate(Ynm_B3(K,PML_Loc), Ynm_g(K,PML_Loc), Ynm_B3_dr(K,PML_Loc), Ynm_B3_dr2(K,PML_Loc/2))
  allocate(AlphaB3(K,K))

  do kk = 1,K                    ! k = Number of Basis Functions in Expansion (from Sav file)
! B3 and B2 are on the same grid in u1. South and north hemis are the same
    do ii = 1,PML_Loc
      Ynm_B3(kk,ii) = Ynm(kk,PML_Loc+1-ii)            ! Reverse order of Basis Vectors due to U1 direction (0 = low boundary N1-1 = upper boundary)
      Ynm_B3_dr(kk,ii)=Ynm_S(kk,PML_Loc+1-ii)
      Ynm_g(kk,ii)=Ynm(kk,PML_Loc+1-ii)*(2.d0*zeros(kk)+1.d0)/(zeros(kk)+1.d0)
    enddo
    do ii = 1,PML_Loc/2
      Ynm_B3_dr2(kk,ii) = Ynm_S(kk,PML_Loc+1-(2*ii))  ! Reverse order and use 2nd value of Basis Vectors due to U1 direction (0 = low boundary N1-1 = upper boundary)
    enddo
  enddo

!  call get_time(tfrom,0,IDlog)

! Time Driver function on outer field line
! Set to gaussian spatial function
  allocate (Tdr(Num_u3))
  Tdr = 0.0
  do ii = 1, Num_u3, 2
    zz = (D_out(ii)-D_out(Num_u3_2))/hwidth
!    Tdr(ii) = 0.5*zz*exp(-zz**2) + exp(-zz**2)   ! odd + even function driver
    Tdr(ii) = exp(-zz**2)   ! even function (gauss) driver
  enddo
  mxDrv = maxval(abs(Tdr))
  Tdr = Tdr/mxDrv

! generate shift index arrays
  allocate (S1L(Num_u1), S1R(Num_u1), S3L (Num_u3), S3R(Num_u3)) 
  do ii = 1, Num_u1
    S1L(ii) = ii + 1
    S1R(ii) = ii - 1
  enddo
  S1R(1) = Num_u1
  S1L(Num_u1) = 1

  do ii = 1, Num_u3
    S3L(ii) = ii + 1
    S3R(ii) = ii - 1
  enddo
  S3R(1) = Num_u3
  S3L(Num_u3) = 1

  allocate (pS1L(PML_Loc), pS1R(PML_Loc)) 
  do ii = 1, PML_Loc
    pS1L(ii) = ii + 1
    pS1R(ii) = ii - 1
  enddo
  pS1R(1) = PML_Loc
  pS1L(PML_Loc) = 1

  print*,'Allocate E and B field soln arrays...'
  allocate (D_u1(Num_u1,Num_u3), D_u3(Num_u1,Num_u3))
  D_u1 = u1Arr(S1L(:),:) - u1Arr(S1R(:),:)
  D_u3 = u3Arr(:,S3L(:)) - u3Arr(:,S3R(:))

  allocate (drv_B3(Num_u3))
  allocate (B1_con(Num_u1,Num_u3), B2_con(Num_u1,Num_u3), B3_con(Num_u1,Num_u3))
  allocate (B1_cov(Num_u1,Num_u3), B2_cov(Num_u1,Num_u3), B3_cov(Num_u1,Num_u3))
  allocate (E1_con(Num_u1,Num_u3), E2_con(Num_u1,Num_u3))
  allocate (E1_cov(Num_u1,Num_u3), E2_cov(Num_u1,Num_u3))

  allocate (Av_B1(Num_u1,Num_u3), Av_B3(Num_u1,Num_u3))

  allocate (B1_N(PML_Loc), B2_N(PML_Loc), B3_N(PML_Loc/2))
  allocate (B1_S(PML_Loc), B2_S(PML_Loc), B3_S(PML_Loc/2))
  allocate (E1_N(PML_Loc), E2_N(PML_Loc))
  allocate (E1_S(PML_Loc), E2_S(PML_Loc))

  allocate (Thi_N(PML_Loc))
  allocate (Thi_S(PML_Loc))
  allocate (DThi_dTh_N(PML_Loc), DThi_dPh_N(PML_Loc))
  allocate (DThi_dTh_S(PML_Loc), DThi_dPh_S(PML_Loc))
  allocate (DThi_dr_N(PML_Loc), DThi_dr_S(PML_Loc))
  allocate (DThi_dTh_N_g(PML_Loc), DThi_dPh_N_g(PML_Loc))
  allocate (DThi_dTh_S_g(PML_Loc), DThi_dPh_S_g(PML_Loc))
  allocate (J_Th_N(PML_Loc), J_Ph_N(PML_Loc))
  allocate (J_Th_S(PML_Loc), J_Ph_S(PML_Loc))
  allocate (E_Th_N(PML_Loc), E_Ph_N(PML_Loc))
  allocate (E_Th_S(PML_Loc), E_Ph_S(PML_Loc))
  allocate (Thi_gnd_N(PML_Loc))
  allocate (Thi_gnd_S(PML_Loc))

  allocate (Coeffs_B3_N(K))
  allocate (Coeffs_B3_S(K))
  allocate (IPIV(K))

  allocate (b_nu_a(Num_u1_2,Num_u3_2))
  allocate (b_ph_a(Num_u1_2,Num_u3_2))
  allocate (b_mu_a(Num_u1_2,Num_u3_2+1))
  allocate (e_nu_a(Num_u1_2,Num_u3_2+1))
  allocate (e_ph_a(Num_u1_2,Num_u3_2+1))

  LDB = K
  IPIV = 0
! Ensure solution arrays are zero
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
  allocate(inv_arr(K,K))
  allocate(b3coef_mult(K,PML_Loc/2))
  allocate(WORK(K*120))

! Calc matrix to multiply dB_r fields at ionosphere for
! atmosphere solution
! Do matrix inverse outside time loop
  AlphaB3 = MatMul(Ynm_B3_dr, transpose(Ynm_B3_dr))
  inv_arr = AlphaB3
  call dgetrf(K, K, inv_arr, K, IPIV, INFO)
  if (info.eq.0) then
    call dgetri(K,inv_arr,K,IPIV,WORK,K*120,INFO)
  else
    print*,'## ERR: Inverse matrix prob: AlphaB3'
    stop
  endif
  deallocate(WORK)

  b3coef_mult = matmul(inv_arr,Ynm_B3_dr2)

!===================================================
! TEST
!do ii=1,k
!  do jj=1,num_u1
!     write(1,*) ynm_b3_dr(ii,jj)
!  enddo
!enddo
!
!do ii=1,k
!  do jj=1,num_u1
!     write(2,*) ynm(ii,jj)
!  enddo
!enddo
!
!do ii=1,k
!  do jj=1,k
!     write(3,*) alphab3(ii,jj),' ',inv_arr(ii,jj)
!  enddo
!enddo
!
!do ii=1,k
!  do jj=1,num_u1_2
!     write(4,*) b3coef_mult(ii,jj)
!  enddo
!enddo
!print*,'stopped'
!stop
! ================================================

  deallocate (inv_arr)

  allocate(emult(num_u1,num_u3))
  allocate(bmult(num_u1,num_u3))
  emult = 2.0*dt*V2/Jacb         ! E field multiplier
  bmult = 2.0*dt/Jacb            ! B field multiplier

! ----------------------------------------------------
!              Start Time Iteration Loop
! ----------------------------------------------------
  print*,'Start time iteration loop...'
  count1 = 0
  do tcnt = 1,Nt-1                          ! Start of time iteration loop
    call cpu_time(t1)

    tt = tcnt-1
    time = 2.d0*Dt*tt                     ! Time counter

    strt_ramp = exp(1.d0-t_ramp/(time+dt))/exp(1.d0)
    bfa = drv_amp*257.d0*time/drv_norm*(1.d0-time/drv_norm)* &
          exp(-5.d0*time/drv_norm)*strt_ramp

! Apply driver on outer field line on B3
    drv_B3 = Tdr*bfa*Re_m*h_mu_outer
    B3_cov(iL_src,1:Num_u3) = -drv_B3 + B3_cov(iL_src,1:Num_u3)

!    B3_cov(Num_u1,1) = 0.d0
!    B3_cov(Num_u1,Num_u3) = 0.d0

!   Conducting Northern Ionosphere
    E1_cov(1:PML_Loc,1) = E1_N(1:PML_Loc)
    E2_cov(1:PML_Loc,1) = E2_N(1:PML_Loc)

    E1_cov(PML_Loc+1:Num_u1,1) = 0.0d0
    E2_cov(PML_Loc+1:Num_u1,1) = 0.0d0

!   Conducting Southern Ionosphere
    E1_cov(1:PML_Loc,Num_u3) = E1_S(1:PML_Loc)
    E2_cov(1:PML_Loc,Num_u3) = E2_S(1:PML_Loc)

    E1_cov(PML_Loc+1:Num_u1,Num_u3) = 0.0d0
    E2_cov(PML_Loc+1:Num_u1,Num_u3) = 0.0d0

!   Inner L Shell Boundary - perfect reflect condition
! B1 is essentially 'radial'
    E2_cov(1,1:Num_u3) = 0.d0
    B1_cov(1,1:Num_u3) = 0.d0

!   E contravariant field iteration
!    E1_con = 2.d0*dt*V2/Jacb*(im*m_num*B3_cov - &
!        (B2_cov(:,S3R(:))-B2_cov(:,S3L(:)))/D_u3) + E1_con
!    E2_con = 2.d0*dt*V2/Jacb*((B1_cov(:,S3R(:)) - B1_cov(:,S3L(:)))/ &
!        D_u3 - (B3_cov(S1R(:),:)-B3_cov(S1L(:),:))/D_u1) + E2_con
    E1_con = NuE1*emult*(im*m_num*B3_cov - &
        (B2_cov(:,S3L(:))-B2_cov(:,S3R(:)))/D_u3) + NuE12*E1_con
    E2_con = NuE2*emult*((B1_cov(:,S3L(:)) - B1_cov(:,S3R(:)))/ &
         D_u3 - (B3_cov(S1L(:),:)-B3_cov(S1R(:),:))/D_u1) + NuE22*E2_con

!   Evolve Contravariant (Tangent) vector fields to Covariant vectors fields
    E1_cov = E1_con/g11_con
    E2_cov = E2_con*g22_cov

!   Northern Ionosphere
    E1_cov(1:PML_Loc,1)= E1_N(1:PML_Loc)
    E2_cov(1:PML_Loc,1)= E2_N(1:PML_Loc)      ! Conducting Northern Ionosphere using BC on Ionosphere

!   Southern Ionosphere
    E1_cov(1:PML_Loc,Num_u3) = E1_S(1:PML_Loc)
    E2_cov(1:PML_Loc,Num_u3) = E2_S(1:PML_Loc)      ! Conducting Southern Ionosphere using BC on Ionosphere

!   Inner L shell Boundary
    E2_cov(1,:) = 0.d0  ! Perfectly Reflecting
    E2_con(1,:) = 0.d0        ! inner L shell field line

! clean up after shift, this means the ionos field cannot come from here
    E2_con(:,1) = 0.d0        ! E2_con is populated as E2_CON(*,evens)
    E2_con(:,Num_u3) = 0.d0   ! other end

    B1_cov(1,:)= 0.d0
    E1_con(1,:) = 0.d0
    E1_con(:,1) = 0.d0
    E1_con(:,Num_u3) = 0.d0

!   B contravariant field iteration
!    B1_con = 2.d0*dt/Jacb*((E2_cov(:,S3R(:))-E2_cov(:,S3L(:)))/ &
!         D_u3)+B1_con
!    B2_con = -2.d0*dt/Jacb*((E1_cov(:,S3R(:))-E1_cov(:,S3L(:)))/ &
!         D_u3)+B2_con
!    B3_con = 2.d0*dt/Jacb*(im*m_num*E1_cov-(E2_cov(S1R(:),:)- &
!         E2_cov(S1L(:),:))/D_u1) + B3_con
    B1_con = NuB1*bmult*((E2_cov(:,S3L(:))-E2_cov(:,S3R(:)))/ &
         D_u3) + NuB12*B1_con
    B2_con = -NuB2*bmult*((E1_cov(:,S3L(:))-E1_cov(:,S3R(:)))/ &
         D_u3) + NuB22*B2_con
    B3_con = NuB3*bmult*(im*m_num*E1_cov-(E2_cov(S1L(:),:)- &
         E2_cov(S1R(:),:))/D_u1) + NuB32*B3_con

!  Outer Corner Point Ionopshere and Outer L shell
!  Set 1st deriv=0, assuming basis expansion has 1st derv=0 functions
    B3_con(PML_Loc,1) = 4.d0/3.d0*B3_con(PML_Loc-2,1) - &
        1.d0/3.d0*B3_con(PML_Loc-4,1)
    B3_con(PML_Loc,Num_u3) = 4.d0/3.d0*B3_con(PML_Loc-2,Num_u3) - &
        1.d0/3.d0*B3_con(PML_Loc-4,Num_u3)

!   Evolve Contravariant Interior Grid B fields to Convariant vectors fields
    Av_B3 = (B3_con(S1R(:),S3R(:)) + B3_con(S1R(:),S3L(:)) + &
        B3_con(S1L(:),S3R(:)) + B3_con(S1L(:),S3L(:))) / 4.d0
    Av_B3(:,1) = 0.0d0
    Av_B3(:,Num_u3) = 0.0d0
    B1_cov = B1_con*g11_cov + Av_B3*g13_cov
    B1_cov(:,1) = 0.d0
    B1_cov(:,Num_u3) = 0.d0                    ! Clean up after shifts due to odd number of points

!   Inner L shell Boundary
    B1_cov(1,:) = 0.d0

!   Evolve B2_cov
    B2_cov = B2_con*g22_cov
    B2_cov(:,1) = 0.0d0
    B2_cov(:,Num_u3) = 0.0d0   ! Clean up after shifts due to odd number of points

!   Evolve B3_cov
    Av_B1 = (B1_con(S1R(:),S3R(:)) + B1_con(S1R(:),S3L(:)) + &
             B1_con(S1L(:),S3R(:)) + B1_con(S1L(:),S3L(:))) / 4.d0
    Av_B1(:,1) = 0.0
    Av_B1(:,Num_u3) = 0.0d0                     ! Clean up after shifts due to odd number of points
    Av_B1(:,2) = 0.0d0
    Av_B1(:,Num_u3-1)=0.0d0
    B3_cov = Av_B1*g13_cov + B3_con*g33_cov

!   Along Northern Ionospheric Boundary
    Av_B1 = (B1_con(S1R(:),S3L(:)) + B1_con(S1L(:),S3L(:))) / 2.0d0
    B3_cov(:,1) = Av_B1(:,1)*g13_cov(:,1) + B3_con(:,1)*g33_cov(:,1)

!   Along Southern Ionosphere Boundary
    Av_B1 = (B1_con(S1R(:),S3R(:)) + B1_con(S1L(:),S3R(:))) / 2.0d0
    B3_cov(:,Num_u3) = Av_B1(:,Num_u3)*g13_cov(:,Num_u3) + B3_con(:,Num_u3)* &
                       g33_cov(:,Num_u3)

! TEST
!    If (mod(tcnt,plot_freq) == 0) then
!      print*,B3_cov(206,num_u3-3:num_u3)
!    endif

!   Along Inner L shell Boundary
    Av_B1 = (B1_con(S1L(:),S3R(:)) + B1_con(S1L(:),S3L(:))) / 2.0d0
    B3_cov(1,:) = Av_B1(1,:)*g13_cov(1,:) + B3_con(1,:)*g33_cov(1,:)

!   Along Outer L shell Boundary
    B3_cov(Num_u1,:) = 0.0d0    ! should be zero after PML

    call cpu_time (t2)

    print*, "prelims: ", (t2-t1), "s"

! Ionosphere Boundary Condition
!  Using B3_con in the ionosphere and perfectly conducting ground (B3_con = 0)
!  solve Laplaces equation in Neutral Atmopshere
!  Calculate B1_cov and B2_cov just below the Ionospheric current sheet.
!  Calculate jump in B1 and B2 and calculate Current in Sheet
!  Invert to find E1 and E2 values in sheet and use as BC in next iteration.
!  Calculate Coefficients in expansion of Bz at the ionosphere (with
!    derivative radial component scale at r = Ri)
!  At ionosphere boundary:
!   Contravariant at ionos boundary for BB3 is radial
!   Covariant at ionos boundary for BB1(colat) and BB2(azimuth) is in the sheet
!
! FOR NORTHERN HEMISPHERE
    do ii = 1,PML_Loc/2
      B3_N(ii)= B3_con(2*ii,1) * h_ra_n(2*ii) ! Use only those points directly calculated
    enddo
! B3_N is a complex array
! Now fit Br data to pV/pr basis set
! inv_arr comes from alphaB3, before the loop
    Coeffs_B3_N = matmul(b3coef_mult,B3_N)
    DThi_dr_N = matmul(Coeffs_B3_N,Ynm_B3_dr) ! This should be B3_N at 2x the spatial resolution
    Thi_N = matmul(Coeffs_B3_N,Ynm_B3)          ! Potential Function
! Calc pV/p_Theta
    do jj=2,PML_Loc-1 
      DThi_dTh_N(jj) = (Thi_N(jj+1)-Thi_N(jj-1))/(-2.d0*del_Th)
    enddo

! Do the end points
!    DThi_dTh_N(1) = (-3.0*Thi_N(1)+4.d0*Thi_N(2)-Thi_N(3))/(-2.d0*del_th)
!    DThi_dTh_N(PML_Loc) = (3.0*Thi_N(PML_Loc)-4.d0*Thi_N(PML_Loc-1) + &
!        Thi_N(PML_Loc-2))/(-2.d0*del_th)
! Use 1st deriv=0 conditions - assuming 1st deriv=0 for basis functions
    DThi_dTh_N(1)=4.0/3.0*Thi_N(2) - 1.0/3.0*Thi_N(3)
    DThi_dTh_N(PML_Loc)=4.0/3.0*Thi_N(PML_Loc-1) - 1.0/3.0*Thi_N(PML_Loc-2)

    DThi_dTh_N = 1.d0/RI_s*DThi_dTh_N        ! B_Theta Atmos
    DThi_dPh_N = CMPLX(0.0,1.d0)*m_num*Thi_N/(RI_s*sin(Colat_N))

! Interpolate B1 and B2 to ALL points just above the ionosphere from model
! solution to calculate currents in Ionosphere
    B1_N = B1_cov(1:PML_Loc,2)/h_th_n
    B2_N = B2_cov(1:PML_Loc,2)/h_ph_n
    B1_N = (B1_N(pS1R(1:PML_Loc))+B1_N(pS1L(1:PML_Loc)))/2.0 + B1_N ! Linear Interpolation
    B2_N = (B2_N(pS1R(1:PML_Loc))+B2_N(pS1L(1:PML_Loc)))/2.0 + B2_N ! Linear Interpolation

! Try second order Taylor expansion for end points
    h = CoLat_N(PML_Loc-1) - CoLat_N(PML_Loc-3) ! step size along ionosphere is uniform on B1_N grid
    FirD = (B1_N(PML_Loc-5) - 4.0*B1_N(PML_Loc-3) + &
        3.0*B1_N(PML_Loc-1))/(2.0*h)            ! First Derivative Backwards difference O(h^2)
    SecD = (-B1_N(PML_Loc-7) + 4.0*B1_N(PML_Loc-5) - &
        5.0*B1_N(PML_Loc-3) + 2.0*B1_N(PML_Loc-1))/(h**2)       ! Second Derivative Backwards difference O(h^2)
!    h1 = CoLat_N(PML_Loc-1)-CoLat_N(PML_Loc)    ! Step size on B1_n_interp grid
! B1_N is defined at all points
!    B1_N(PML_Loc) = B1_N(PML_Loc-1) + h1*FirD + ((h1**2) / 2.0)*SecD
    B1_N(PML_LOc) = B1_N(PML_Loc-1) + h*FirD + ((h**2)/2.0)*SecD

    h = CoLat_N(4)-CoLat_N(2) !; step size along ionosphere is uniform
    FirD = (-B2_N(6) + 4.0*B2_N(4) - 3.0*B2_N(2)) / (2.0*h)       ! First Derivative Forwards difference O(h^2)
    SecD = (-B2_N(8) + 4.0*B2_N(6) - 5.0*B2_N(4) + 2.0*B2_N(2)) / (h**2)   ! Second Derivative Forwards difference O(h^2)
!    h1 = CoLat_N(1) - CoLat_N(2)
! B2_N is defined at all points
!    B2_N(1)= B2_N(2) + h1*FirD + ((h1**2) / 2.0)*SecD
    B2_N(1)= B2_N(2) + h*FirD + ((h**2)/2.0)*SecD

! Currents in the Ionosphere
    J_Th_N = -(B2_N - DThi_dPh_N) / mu0_s
    J_Ph_N =  (B1_N - DThi_dTh_N) / mu0_s

! Calculate electric fields in Ionosphere from discontinuity in B's for
! next time step...
    E_Th_N = (INVSS_N(1:PML_Loc,1,1) * J_Th_N + INVSS_N(1:PML_Loc,2,1) * J_Ph_N)
    E_Ph_N = (INVSS_N(1:PML_Loc,1,2) * J_Th_N + INVSS_N(1:PML_Loc,2,2) * J_Ph_N)

    E1_N = E_Th_N*h_th_n*evens(1:PML_Loc)
    E2_N = E_Ph_N*h_ph_n*odds(1:PML_Loc)

! Perfect reflection condition
!    E1_N(1:num_u1) = 0.d0
!    E2_N(1:num_u1) = 0.d0

!   Now do the Ground fields
!    Thi_gnd_N = matmul(transpose(Ynm_g), Coeffs_B3_N)                 ! Potential Function
    Thi_gnd_N = matmul(Coeffs_B3_N,Ynm_g)                 ! Potential Function
    do jj = 2,PML_Loc-1 
      DThi_dTh_N_g(jj) = (Thi_gnd_N(jj+1) - Thi_gnd_N(jj-1)) / (-2.d0*del_Th)
    enddo

! Do the end points
!    DThi_dTh_N_g(1) = (-3.0*Thi_gnd_N(1) + 4.d0*Thi_gnd_N(2) - Thi_gnd_N(3))/(-2.d0*del_th)
!    DThi_dTh_N_g(PML_Loc) = (3.0*Thi_gnd_N(PML_Loc) - 4.d0*Thi_gnd_N(PML_Loc-1) + Thi_gnd_N(PML_Loc-2)) / (-2.d0*del_th)
! 1st deriv=0 condition
    DThi_dTh_N_g(1) = 4.0/3.0*Thi_gnd_N(2) - 1.0/3.0*Thi_gnd_N(3)
    DThi_dTh_N_g(PML_Loc) = 4.0/3.0*Thi_gnd_N(PML_Loc-1) - 1.0/3.0*Thi_gnd_N(PML_Loc-2)
    DThi_dTh_N_g = 1.d0/RI_s*DThi_dTh_N_g                                ! B_Theta Ground
    DThi_dPh_N_g = Cmplx(0.0,1.d0)*m_num*Thi_gnd_N / (Re_s*sin(Colat_N)) !B_phi Ground

    call cpu_time(t3)

    print*, "Northern Hemisphere: ", (t3-t1), "s"

!  SOUTHERN HEMISPHERE
!   Linear Interpolation
    do ii = 1,PML_Loc/2
      B3_S(ii) = B3_con(2*ii,Num_u3)*h_ra_s(2*ii)
    enddo

!   Fit Br data to pV/pr basis set
    Coeffs_B3_S = matmul(b3coef_mult,B3_S)
    DThi_dr_S = matmul(Coeffs_B3_S,Ynm_B3_dr)        ! This should be B3_S at twice the spatial resolution
    Thi_S = matmul(Coeffs_B3_S,Ynm_B3)                    ! Potential Function

!   Calc. pV/p_Theta
    do jj=2,PML_Loc-1 
      DThi_dTh_S(jj) = (Thi_S(jj+1)-Thi_S(jj-1))/(2.d0*del_Th)
    enddo

! Do the end points
!    DThi_dTh_S(1) = (-3.d0*Thi_S(1)+4.d0*Thi_S(2)-Thi_S(3))/(2.d0*del_th)
!    DThi_dTh_S(PML_Loc) = (3.d0*Thi_S(PML_Loc) - 4.d0*Thi_S(PML_Loc-1) + &
!        Thi_S(PML_Loc-2))/(2.d0*del_Th)
    DThi_dTh_S(1) = 4.0/3.0*Thi_S(2) - 1.0/3.0*Thi_S(3)
    DThi_dTh_S(PML_Loc) = 4.0/3.0*Thi_S(PML_Loc-1) - 1.0/3.0*Thi_S(PML_Loc-2)
    DThi_dTh_S = 1.d0/RI_s*DThi_dTh_S                             ! B_Theta Atmos
    DThi_dPh_S = Cmplx(0.0,1.d0)*m_num*Thi_S/(RI_s*sin(Colat_S))  ! B_phi Atmos

!   Extrapolate B1 and B2 to ALL points just above ionosphere from model
!   solution to calculate currents in Ionosphere
    B1_S= B1_cov(1:PML_Loc,Num_u3-1)/h_th_s
    B2_S= B2_cov(1:PML_Loc,Num_u3-1)/h_ph_s
    B1_S = (B1_S(pS1R(1:PML_Loc)) + B1_S(pS1L(1:PML_Loc)))/2.0 + B1_S ! Linear Interpolation
    B2_S = (B2_S(pS1R(1:PML_Loc)) + B2_S(pS1L(1:PML_Loc)))/2.0 + B2_S ! Linear Interpolation

!   Try second order Taylor expansion for end points
    h = CoLat_S(PML_Loc-1)-CoLat_S(PML_Loc-3)  ! step size along ionosphere is uniform on B1_N grid
!   First Derivative Backwards difference O(h^2)
    FirD =(B1_S(PML_Loc-5) - 4.0*B1_S(PML_Loc-3) + 3.0*B1_S(PML_Loc-1))/(2.0*h)
!   Second Derivative Backwards difference O(h^2)
    SecD =(-B1_S(PML_Loc-7) + 4.0*B1_S(PML_Loc-5) - 5.0*B1_S(PML_Loc-3)+ &
        2.0*B1_S(PML_Loc-1))/(h**2)
!    h1 = CoLat_S(PML_Loc-1)-CoLat_S(PML_Loc) ! Step size on B1_n_interp grid
!    B1_S(PML_Loc) = B1_S(PML_Loc-1)+h1*FirD+((h1**2)/2.0)*SecD
    B1_S(PML_Loc) = B1_S(PML_Loc-1)+h*FirD+((h**2)/2.0)*SecD

    h = CoLat_S(4)-CoLat_S(2) ! step size along ionosphere is uniform
!   First Derivative Forwards difference O(h^2)
    FirD = (-B2_S(6) + 4.0*B2_S(4) - 3.0*B2_S(2))/(2.0*h)
!   Second Derivative Forwards difference O(h^2)
    SecD = (-B2_S(8) + 4.0*B2_S(6) - 5.0*B2_S(4) + 2.0*B2_S(2))/(h**2)
!    h1 = CoLat_S(1)-CoLat_S(2)
!    B2_S(1) = B2_S(2) + h1*FirD + ((h1**2)/2.0)*SecD
    B1_S(1) = B1_S(2)+h*FirD+((h**2)/2.0)*SecD

! Currents in the Ionosphere
    J_Th_S = -(B2_S - DThi_dPh_S) / mu0_s
    J_Ph_S =  (B1_S - DThi_dTh_S) / mu0_s

! Calculate electric fields in Ionosphere from discontinuity in B's for
! next time step.
    E_Th_S = (INVSS_S(1:PML_Loc,1,1)*J_Th_S + INVSS_S(1:PML_Loc,2,1)*J_Ph_S)
    E_Ph_S = (INVSS_S(1:PML_Loc,1,2)*J_Th_S + INVSS_S(1:PML_Loc,2,2)*J_Ph_S)
    E1_S = E_Th_S*h_th_s*evens(1:PML_Loc)
    E2_S = E_Ph_S*h_ph_s*odds(1:PML_Loc)

!If (mod(tcnt,plot_freq) == 0) then
!  print*,Thi_S(205:209)
!  print*,' '
!  print*,DThi_dr_S(205:209)
!  print*,E1_S(205:209)
!   print*,J_Th_S(205:209)
!print*,''
!  print*,E2_S(205:209)
!   print*,J_Ph_S(205:209)
!print*,' '
!endif

! Perfect reflection condition for testing
!    E1_S(1:num_u1) = 0.d0
!    E2_S(1:num_u1) = 0.d0

! Now do the Ground fields
!    Thi_gnd_S = matmul(transpose(Ynm_g), Coeffs_B3_S)  ! Potential Function
    Thi_gnd_S = matmul(Coeffs_B3_S,Ynm_g)  ! Potential Function
    do jj=2,PML_Loc-1 
      DThi_dTh_S_g(jj) = (Thi_gnd_S(jj+1)-Thi_gnd_S(jj-1))/(-2.d0*del_Th)
    enddo

! Do the end points
!    DThi_dTh_S_g(1) = (-3.0*Thi_gnd_S(1)+4.d0*Thi_gnd_S(2)-Thi_gnd_S(3))/ &
!        (-2.d0*del_th)
!    DThi_dTh_S_g(PML_Loc) = (3.0*Thi_gnd_S(PML_Loc)-4.d0*Thi_gnd_S(PML_Loc-1) + &
!        Thi_gnd_S(PML_Loc-2))/(-2.d0*del_th)
    DThi_dTh_S_g(1) = 4.0/3.0*Thi_gnd_S(2) - 1.0/3.0*Thi_gnd_S(3)
    DThi_dTh_S_g(PML_Loc) = 4.0/3.0*Thi_gnd_S(PML_Loc-1) - 1.0/3.0*Thi_gnd_S(PML_Loc-2)
    DThi_dTh_S_g = 1.d0 / RI_s * DThi_dTh_S_g                           ! B_Theta Ground
    DThi_dPh_S_g = Cmplx(0.0, 1.d0)*m_num*Thi_gnd_S/(Re_s*sin(Colat_S)) ! B_phi Ground

    call cpu_time(t2)

    print*, "Southern Hemisphere: ", (t2-t3), "s"

! =================================================
! Output data arrays
    If (mod(tcnt,plot_freq) == 0) then
      do ii=1,Num_u1_2
        do jj=1,Num_u3_2
          b_nu_a(ii,jj) = 1.0e9*h_nu(2*ii-1,2*jj  )*B1_con(2*ii-1,2*jj  )/Re_m
          b_ph_a(ii,jj) = 1.0e9*h_ph(2*ii,  2*jj  )*B2_con(2*ii  ,2*jj  )/Re_m
          b_mu_a(ii,jj) = 1.0e9*B3_cov(2*ii,2*jj-1)/h_mu(2*ii    ,2*jj-1)/Re_m
          e_nu_a(ii,jj) = 1.0e3*h_nu(2*ii  ,2*jj-1)*E1_con(2*ii  ,2*jj-1)
          e_ph_a(ii,jj) = 1.0e3*h_ph(2*ii-1,2*jj-1)*E2_con(2*ii-1,2*jj-1)
        enddo
      enddo
! do field line end points
      do ii=1,Num_u1_2
        b_mu_a(ii,1) = 1.0e9*B3_cov(2*ii,1)/h_mu(2*ii,1)/Re_m
        e_nu_a(ii,1) = 1.0e3*h_nu(2*ii,1)*E1_con(2*ii,1)
        e_ph_a(ii,1) = 1.0e3*h_ph(2*ii-1,1)*E2_con(2*ii-1,1)
        b_mu_a(ii,num_u3_2+1) = 1.0e9*B3_cov(2*ii,num_u3)/h_mu(2*ii,num_u3)/Re_m
        e_nu_a(ii,num_u3_2+1) = 1.0e3*h_nu(2*ii  ,num_u3)*E1_con(2*ii  ,num_u3)
        e_ph_a(ii,num_u3_2+1) = 1.0e3*h_ph(2*ii-1,num_u3)*E2_con(2*ii-1,num_u3)
      enddo

! Output data section
      if (count1.eq.0) then
!        mk_out_dir = 'mkdir ' // out_dir
!        call system(mk_out_dir)
        write(out_file,'(a,a)') trim(out_dir),'/gridxy.dat'
        open(unit=1,file=out_file,form='unformatted')
        write(1) xarr,yarr, b1x_ord,b1y_ord, &
                            b2x_ord,b2y_ord, &
                            b3x_ord,b3y_ord, &
                            e1x_ord,e1y_ord, &
                            e2x_ord,e2y_ord
        close(unit=1)
      endif           ! If 1st time through here

      count1 = count1 + 1
      write(out_file,'(a,a,i4.4,a4)') trim(out_dir),'/mag_fields_',count1,'.dat'
      open(unit=1,file=out_file,form='unformatted')
      write(1) b_nu_a, b_ph_a, b_mu_a, e_nu_a, e_ph_a
      close(unit=1)

      print*,'count and bfa(nT) = ',count1,' ',bfa*1.d9
      print*,'Max of abs(bmu) = ',maxval(abs(b_mu_a))
!stop
    endif       ! if outout arrays

  enddo  ! ============== End of Time Loop ====================

  print*,'Finished'

end 

