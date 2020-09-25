!---------------------------------------------------------------------
! fmin
!	f = 1/2 F.F needed for a globally convergent
!	NR nonlinear solver 
! - Requires access to BLAS libraries
!---------------------------------------------------------------------
! n	: int, number of variables
! fvec  : 1D real*8, current vector of F 
real(kind=8) function fmin(n,Fvec)
  implicit none
  integer, intent(in) :: n
  real(kind=8), intent(in) :: Fvec(1:n)
  real(kind=8) :: ddot
  fmin = 0.5d0*ddot(n,Fvec(1:n),1,Fvec(1:n),1)
end function fmin
!---------------------------------------------------------------------
