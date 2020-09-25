!---------------------------------------------------------------------
! gcnr
!	- globally convergent NR nonlinear solver
!	  uses the user defined subroutines in fnc.f90 and jac.f90
!	  to compute the F vector to be made zero (fnc) and it's 
!	  jacobian (jac).
!	  
!---------------------------------------------------------------------
! n		: int, number of variables
! m		: int, number of parameters
! x		: 1D real*8, variables 
! prm		: 1D real*8, parameters
! errF		: 1D real*8, function error estimate
! errx		: 1D real*8, variable error estimate
! tolF		: real*8, F convergence tolerance
! tolx		: real*8, varaible convergence tolerance
! maxit		: int, maximum number of iterations
! stat		: int, on exit, 0 if converged
!         	                1 if local converge
!				2 if no convergence
!                               3 if error 
subroutine gcnr(n,m,x,prm,errF,errx,tolF,tolx,maxit,stat)
  implicit none
  integer, intent(in) :: n,m,maxit
  integer, intent(inout) :: stat
  real(kind=8), intent(in) :: prm(1:m),tolF,tolx
  real(kind=8), intent(inout) :: x(1:n)
  real(kind=8), dimension(1:n), intent(inout) :: errF,errx
  integer :: i,itr,ipiv(1:n),itmp
  real(kind=8) :: f,fold,rtmp,rtmp2,den,stpmax
  real(kind=8) :: Fjac(1:n,1:n)
  real(kind=8), dimension(1:n) :: Fvec,g,xold,p
  real(kind=8) :: STPMX,TOLMIN
  parameter(STPMX=100.,TOLMIN=1.d-6)
  real(kind=8) :: ddot
  stat = 2

  !test if initial guess is a root
  f = fmin(n,x(1:n))
  call fnc(n,m,x(1:n),prm(1:m),Fvec(1:n))
  rtmp = 0.d0
  do i=1,n
    if (abs(Fvec(i)) .gt. rtmp) rtmp = abs(Fvec(i))
  end do 
  if (rtmp .lt. 1.d-2*tolF) then
    stat = 0
    errx = 0.d0
    errf = 1.d-2*tolF
    return
  end if

  !calculate max stepsize 
  rtmp = 0.d0 
  rtmp = ddot(n,Fvec(1:n),1,Fvec(1:n),1)
  stpmax = STPMX*max(sqrt(rtmp),float(n))

  !iterations
  do itr=1,maxit 

!    write(*,*) "-------------------------------------"
!    write(*,*) "GC iteration #",itr

    !generate new F vector and Jacobian
    call fnc(n,m,x(1:n),prm(1:m),Fvec(1:n))
    call jac(n,m,x(1:n),prm(1:m),Fjac(1:n,1:n))

    !calculate del f for line search
    call dgemv('T',n,n,1.d0,Fjac(1:n,1:n),n,Fvec(1:n),1,0.d0,g(1:n),1)

    !store x and f
    xold = x
    fold = f

    !solve linear equations
    p(1:n) = -1.d0*Fvec(1:n) !RHS of lineq
    call dgesv(n,1,Fjac(1:n,1:n),n,ipiv(1:n),p(1:n),n,itmp)
    if (itmp .ne. 0) then
      write(*,*) "ERROR ERROR ERROR"
      write(*,*) "gcnr : dgesv exited with status",itmp
      stat = 3
      errX = 100000.
      errF = 100000.
      return
    end if

    !call line search
    ! line search will output if it thinks the convergence 
    ! needs to be checked in the stat variable
    call lnsrch(n,m,xold(1:n),prm(1:m),fold,&
                g(1:n),p(1:n),x(1:n),stpmax,tolx,stat)
    !TESTING TESTING
!    call daxpy(n,1.d0,p(1:n),1,x(1:n),1)
!    write(*,*) "x old =", xold(1:n)
!    write(*,*) "x new =", x(1:n)
    if (stat .ne. 1 .and. stat .ne. 0) then
      write(*,*) "There was an error from lnsrch"
      stat = 3
      errF = 100000.
      errX = 100000.
      return 
    end if
     
    !check convergence
    rtmp = 0.d0
    do i=1,n
      if (abs(Fvec(i)) .gt. rtmp) rtmp = abs(Fvec(i))
    end do
    if (rtmp .lt. tolf) then !converged
      errX(1:n) = abs(x(1:n) - xold(1:n))
      stat = 0
      return
    end if

    !check for spurious (local) convegence, indicated from 
    ! lnsrch
    if (stat .eq. 1) then
      rtmp = 0.d0
      den = max(f,0.5d0*n)
      do i=1,n
        rtmp2 = abs(g(i))*max(abs(x(i)),1.d0)/den 
        if (rtmp2 .gt. rtmp) rtmp = rtmp2
      end do  

      if (rtmp .lt. TOLMIN) then
        stat = 0
      else
        stat = 1
      end if
      return
    end if

    !test convergence on δx
    rtmp = 0.d0
    do i=1,n
      rtmp2 = (abs(x(i)-xold(i)))/max(abs(x(i)),1.d0) 
      if (rtmp2 .gt. rtmp) rtmp = rtmp2
    end do 

    errX(1:n) = abs(x(1:n) - xold(1:n))
    if (rtmp .lt. tolX) return

  end do !iteration loop
  stat = 2

end subroutine gcnr
!---------------------------------------------------------------------
