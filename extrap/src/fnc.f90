!---------------------------------------------------------------------
! fnc
!	- user input subroutine that calculates the F vector
!	  of a function we wish to solve with a given set of 
!	  variables (x) and parameters (prm) 
!
! Here:
!	F = E_X - E_Y + a(1/Y^b - 1/X^b)
!       x = [a,b]
!     prm = [X,Y,Z,E_X,E_Y,E_Z]
!---------------------------------------------------------------------
! n		: int, number of variables
! m		: int, number of parameters
! x		: 1D real*8, variables
! prm		: 1D real*8, parameters
! Fvec		: 1D real*8, F vector
subroutine fnc(n,m,x,prm,Fvec)
  implicit none
  integer, intent(in) :: n,m
  real(kind=8), intent(in) :: x(1:n),prm(1:m)
  real(kind=8), intent(inout) :: Fvec(1:n)
  Fvec(1) = prm(4) - prm(5) + x(1)*(1.d0/prm(2)**x(2) - 1.d0/prm(1)**x(2))
  Fvec(2) = prm(5) - prm(6) + x(1)*(1.d0/prm(3)**x(2) - 1.d0/prm(2)**x(2))
end subroutine fnc
!---------------------------------------------------------------------
