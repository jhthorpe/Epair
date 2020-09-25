!---------------------------------------------------------------------
! nlslv
!	- module containing subrotuines for calculating solutions
!	  to nonlinear equations. Requires BLAS and LAPACK libraries
!---------------------------------------------------------------------
module nlslv

contains
   
  include 'fmin.f90'
  include 'jac.f90'
  include 'fnc.f90'
  include 'lnsrch.f90' 
  include 'gcnr.f90'

end module nlslv
!---------------------------------------------------------------------
