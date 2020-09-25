!---------------------------------------------------------------------
! jac
!	- user input subroutine that calculates the jacobian 
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
! Fjac		: 2D real*8, Jacobian, dFi/dxj
subroutine jac(n,m,x,prm,Fjac)
  implicit none
  integer, intent(in) :: n,m
  real(kind=8), intent(in) :: x(1:n),prm(1:m)
  real(kind=8), intent(inout) :: Fjac(1:n,1:n)
  Fjac(1,1) = (1.d0/prm(2)**x(2) - 1.d0/prm(1)**x(2)) 
  Fjac(1,2) = x(1)*(log(prm(1))/prm(1)**x(2) - log(prm(2))/prm(2)**x(2)) 
  Fjac(2,1) = (1.d0/prm(3)**x(2) - 1.d0/prm(2)**x(2)) 
  Fjac(2,2) = x(1)*(log(prm(2))/prm(2)**x(2) - log(prm(3))/prm(3)**x(2)) 
end subroutine jac
!---------------------------------------------------------------------
