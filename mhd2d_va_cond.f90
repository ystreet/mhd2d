module mhd2d_va_cond
! Two-dimensional MHD code using FDTD algorithm
! C. L. Waters and M. D. Sciffer
! Centre for Space Physics
! University of Newcastle
! New South Wales 
! Australia
!
! March 2012
!
! Contains the following subroutines
! clw_va : default Va profile routine
! get_B0 : dipole field subroutine
! calc_cond_ionos : 

  use mhd2d_constants

  implicit none

  contains

  subroutine clw_va(xd,va_sp)
! Calculates the Alfven speed in the equatorial plane as a function of R
! Input:
!    xd, radial distance value in Re
! Output:
!    va_sp: output Va in m/s
!
! C. L. Waters - see Waters et al., JGR, 2000
!
    real(DBL), intent(in) :: xd
    real(DBL), intent(out) :: va_sp
    real(DBL) :: tanhm, tterm

!    tanhm=4.0*(xd-5.6)                    !   Plasmapause at L = 5.6
!    tterm=(157./xd**2.6)*(tanh(tanhm)+1.0)/1.6

    tanhm=2.0*(xd-4.5)               !   Plasmapause at L = 4.5 (Used for 1st 2D Model Paper)
    tterm=(55.0/xd**2.6)*(tanh(tanhm)+1.0)/2.0

!    tanhm=3.0*(xd-3.5)                    !   Plasmapause at L = 3.5
!    tterm=(27.0/xd**2.6)*(tanh(tanhm)+1.0)/1.6

!    tanhm=2.0*(xd-14.5)                   !   Plasmapause at L = 14.5 (ie. No Plasmapause)
!    tterm=(55.0/xd**2.6)*(tanh(tanhm)+1.0)/2.0

    va_sp=(2.09/xd**1.1+tterm)*990.0-50.0
    va_sp=va_sp*1000.0

    return
  end subroutine clw_va
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  subroutine Get_B0(Lat,LVal,b_0)
! Dipole magnetic field versus latitude
! Inputs:
!   Lat   ! Latitude in radians (not Co_Lat)
!   L     ! L Value of field line
! Outputs:
!   B_0   ! Magnetic field in Tesla
!
    use mhd2d_constants
    
    real(DBL), intent(in) :: Lat, LVal
    real(DBL), intent(out) :: B_0
    real(DBL) :: K_0   ! Earth magnetic moment

    K_0 = 8.0d15
    b_0=sqrt(4.d0-3.d0*cos(Lat)**2.)/(cos(Lat)**6.)*K_0/(LVal*rE_m)**3

    return
  end subroutine Get_B0
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  Subroutine calc_cond_ionos(Num_u1,CoLat_N,CoLat_S, &
                     Sd_N,Sp_N,Sh_N,INVSS_N, &
                     Sd_S,Sp_S,Sh_S,INVSS_S, &
                     h_ph_n,h_th_n,h_ra_n, &
                     h_ph_s,h_th_s,h_ra_s)
!
!   Calc height integrated ionosphere conductivities
! Inputs:
!   Num_u1 : num points used for latitude grid
!   CoLat_N, CoLat_S : co-latitudes (rad) array(Num_u1)
! Outputs:
!   Sd_N,Sp_N,sh_N : conductance values (direct, Pedersen, Hall) for North
!   Sd_S,Sp_S,Sh_S : conductance values for south hemisphere
!   INVSS_N,INVSS  : inverse of conductance matrix
!   h_ph, h_th, h_ra : scale factors to get to spherical coords
!
    use mhd2d_constants

    implicit none
    save

    integer, intent(in) :: Num_u1
    real(DBL), dimension(Num_u1), intent(in) :: CoLat_N, CoLat_S
    real(DBL), dimension(Num_u1), intent(out) :: h_ph_n, h_th_n,h_ra_n
    real(DBL), dimension(Num_u1), intent(out) :: h_ph_s, h_th_s,h_ra_s
    real(DBL), dimension(Num_u1,2,2), intent(out) :: invss_n, invss_s
    real(DBL), dimension(Num_u1), intent(inout) :: Sd_N, Sp_N, Sh_N
    real(DBL), dimension(Num_u1), intent(inout) :: Sd_S, Sp_S, Sh_S

    integer :: ii, M,N, LDA, INFO, LWORK, mm, nn
    real(DBL), dimension(2,2) :: SS, Inverse_SS, mat
    double precision, dimension(2) :: WORK
    integer, dimension(2) :: IPIV
!    complex(8) :: Sig0, Sig1, Sig2
    real(DBL) :: Theta_0, Theta, I, Alpha, Sz_N, Sz_S

!    Sig0 = dcmplx(0.0,0.0)
!    Sig1 = Sig0
!    Sig2 = Sig0
!
! Variables used for matrix inverse routine
    M = 2
    N = 2
    LDA = 2
    LWORK = 2

    do ii = 1, Num_u1 
      Theta  = Colat_N(ii)                    ! Co-latitude
      I = pi/2.d0-dacos(-2.d0*cos(Theta)/(sqrt(1.d0+3.d0*cos(theta)**2)))    ! Dip Angle (as per sciffer and waters 2002)
! TEST conductances
!      Sd_N(ii) = Double(Sig0)*Re_m**2                    ! Height Integrated Direct Conductivity in Northern Hemipshere
!      Sp_N(ii) =  Double(Sig1)*Re_m**2                    ! Height Integrated Pederson Conductivity in Northern Hemipshere
!      Sh_N(ii) =  Double(Sig2)*Re_m**2                    ! Height Integrated Hall Conductivity in Noprthern Hemipshere

      Sd_N(ii) = Sd_N(ii)*Re_m**2       ! Height Integrated Direct Conductivity, Nth Hemisphere
      Sp_N(ii) = Sp_N(ii)*Re_m**2       ! Height Integrated Pedersen Conductivity in Northern Hemipshere
      Sh_N(ii) = Sh_N(ii)*Re_m**2       ! Height Integrated Hall Conductivity in Noprthern Hemipshere

      Theta_0 = Colat_N(ii)
      h_ph_n(ii) =  RI_s*sin(Theta)
      h_th_n(ii) = -RI_s/(2.d0*sin(Theta_0)*cos(Theta_0))
      h_ra_n(ii) = -(2.d0*RI_s*(cos(Theta_0))**2)/(1.d0+3.d0*(cos(theta_0))**2)

      Alpha = dacos(-2.d0*cos(Theta)/(sqrt(1.d0+3.d0*cos(theta)**2)))
      Sz_N = Sd_N(ii)*cos(Alpha)*cos(Alpha) + Sp_N(ii)*sin(Alpha)*sin(Alpha)

      SS(1,1) =  Sd_N(ii)*Sp_N(ii)/Sz_N                       ! Lysak 2004
      SS(2,1) = -Sd_N(ii)*Sh_N(ii)*cos(Alpha)/Sz_N
      SS(1,2) =  Sd_N(ii)*Sh_N(ii)*cos(Alpha)/Sz_N
      SS(2,2) =  Sp_N(ii)+(Sh_N(ii)**2*sin(Alpha)**2)/Sz_N

! Inverse of Conductivity Tensor used to compute E's in Ionosphere
!      call DLincg(2,SS,2,Inverse_SS,2)
      mat = SS
      call dgetrf(M,N,mat,LDA,IPIV,INFO)
      call dgetri(N,mat,LDA,IPIV,WORK,LWORK,INFO)
      do mm=1,2
        do nn=1,2
          InvSS_N(ii,mm,nn) = mat(mm,nn)
        enddo
      enddo
    enddo            ! end of Num_u1 loop
!
! South Hemisphere
    do ii = 1, Num_u1 
      Theta  = Colat_S(ii)                    ! Co-latitude
      I = pi/2.d0-dacos(-2.d0*cos(Theta)/(sqrt(1.d0+3.d0*cos(theta)**2)))    ! Dip Angle (as per sciffer and waters 2002)
! TEST conductances
!      Sd_S(ii) = Double(Sig0)*Re_m**2                    ! Height Integrated Direct Conductivity in Northern Hemipshere
!      Sp_S(ii) =  Double(Sig1)*Re_m**2                    ! Height Integrated Pederson Conductivity in Northern Hemipshere
!      Sh_S(ii) =  Double(Sig2)*Re_m**2                    ! Height Integrated Hall Conductivity in Noprthern Hemipshere

      Sd_S(ii) = Sd_S(ii)*Re_m**2       ! Height Integrated Direct Conductivity, Nth Hemisphere
      Sp_S(ii) = Sp_S(ii)*Re_m**2       ! Height Integrated Pedersen Conductivity in Northern Hemipshere
      Sh_S(ii) = Sh_S(ii)*Re_m**2       ! Height Integrated Hall Conductivity in Noprthern Hemipshere

      Theta_0 = Colat_N(ii)
      h_ph_s(ii) =  RI_s*sin(Theta)
      h_th_s(ii) = -RI_s/(2.d0*sin(Theta_0)*cos(Theta_0))
      h_ra_s(ii) = -(2.d0*RI_s*(cos(Theta_0))**2)/(1.d0+3.d0*(cos(theta_0))**2)

      Alpha = dacos(-2.d0*cos(Theta)/(sqrt(1.d0+3.d0*cos(theta)**2)))
      Sz_S = Sd_S(ii)*cos(Alpha)*cos(Alpha) + Sp_S(ii)*sin(Alpha)*sin(Alpha)

      SS(1,1) =  Sd_S(ii)*Sp_S(ii)/Sz_S                       ! Lysak 2004
      SS(2,1) = -Sd_S(ii)*Sh_S(ii)*cos(Alpha)/Sz_S
      SS(1,2) =  Sd_S(ii)*Sh_S(ii)*cos(Alpha)/Sz_S
      SS(2,2) =  Sp_S(ii)+(Sh_S(ii)**2*sin(Alpha)**2)/Sz_S

! Inverse of Conductivity Tensor used to compute E's in Ionosphere
!      call DLincg(2,SS,2,Inverse_SS,2)
!      InvSS_S(ii,:,:) = Inverse_SS(:,:)
      mat = SS
      call dgetrf(M,N,mat,LDA,IPIV,INFO)
      call dgetri(N,mat,LDA,IPIV,WORK,LWORK,INFO)
      do mm=1,2
        do nn=1,2
          InvSS_S(ii,mm,nn) = mat(mm,nn)
        enddo
      enddo

    enddo            ! end of Num_u1 loop

  end Subroutine calc_cond_ionos
!
end module mhd2d_va_cond
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
