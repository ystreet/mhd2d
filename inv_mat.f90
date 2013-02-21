! test inverse of a matrix (square and real) using La_pack
! C.L. Waters
! March 2012
!
 program inv_mat

   use mhd2d_constants

   implicit none

!   double precision, dimension(2,2) ::  mat, invmat
   real(DBL), dimension(2,2) :: mat, invmat
   double precision, dimension(2) :: WORK
   integer, dimension(2) :: IPIV
   integer :: M, N, LDA, INFO, LWORK

   mat(1,1) = 1.0
   mat(1,2) = 2.0
   mat(2,1) = 3.0
   mat(2,2) = 4.0

   M = 2
   N = 2
   LDA = 2
   invmat = mat
   print*,'MAT '
   print*,invmat(1,1:2)
   print*,invmat(2,1:2)

   print*,'calling DGETRF...'
   call dgetrf(M,N,invmat,LDA,IPIV,INFO)
   print*,'INFO=',INFO
   print*,'MAT '
   print*,invmat(1,1:2)
   print*,invmat(2,1:2)

   LWORK = 2
   print*,'calling DGETRI...'
   call dgetri (N,invmat,LDA,IPIV,WORK,LWORK,INFO)
   print*,'inv_MAT '
   print*,invmat(1,1:2)
   print*,invmat(2,1:2)

end program inv_mat

