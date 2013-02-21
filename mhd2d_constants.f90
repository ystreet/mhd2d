module mhd2d_constants
!
! Two-dimensional MHD solution using FDTD algorithm
! C. L. Waters & M. D. Sciffer
! Centre for Space Physics
! University of Newcastle
! Australia
!
! Constants used in the 2D MHD code
! Last reviewed: feb 2013
!
  implicit none
  
  integer, parameter :: SGL = selected_real_kind ( p=6, r=37)
  integer, parameter :: DBL = selected_real_kind ( p=13, r=200)

  real(kind=DBL), parameter :: rE_m = 6371.0e3    ! metres
  real(kind=DBL), parameter :: pi = 3.1415926535897932384626
  real(kind=DBL), parameter :: DegToRad = pi / 180.0d0
  real(kind=DBL), parameter :: RadToDeg = 180.0d0/pi
  real(kind=DBL), parameter :: mu0 = 4.0d0*pi*1.0d-7
  real(kind=DBL), parameter :: eps0 = 8.8541878d-12
  real(kind=DBL), parameter :: c = sqrt(1.0d0/(mu0*eps0))
  real(kind=DBL), parameter :: im = dcmplx(0.0d0,1.0d0)

  real(kind=DBL), parameter :: mu0_s = mu0/rE_m
  real(kind=DBL), parameter :: eps0_s = eps0*rE_m**3
  real(kind=DBL), parameter :: c_s = c/rE_m
  real(kind=DBL), parameter :: rE_s = 1.0d0
  real(kind=DBL), parameter :: rI_s = 1.02d0   ! ionos height in Re

end module mhd2d_constants
