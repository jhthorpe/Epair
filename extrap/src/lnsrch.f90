!---------------------------------------------------------------------
! lnsrch
!	- subroutine that takes in a set of variables (x) and 
!	  parameters (prm) of some function vector (F) to be solved,
!	  which has some initial vector (p). lnsrch requires the 
!	  descent direction (g), which points down the function 
!	  f = 1/2 F.F   
!
!	- the calling function should check that the convegence is
! 	  not spurious
!---------------------------------------------------------------------
! n		: int, number of variables
! m		: int, number of parameters
! xold		: 1D real*8, input vector of variables
! prm		: 1D real*8, input vector of parameters
! fold		: real*8, input f function
! g		: 1D real*8, input step along f 
! p		: 1D real*8, input direction from xold
! x		: 1D real*8, output new x
! stpmax	: real*8, input limit of stepsize
! tolx		: real*8, convergence criteria of x
! stat		: int, on exit, 0 if converged
! 			 	1 if locally conveged
!				2 if not converged
!				3 if error
subroutine lnsrch(n,m,xold,prm,fold,g,p,x,stpmax,tolx,stat)
  implicit none
  integer, intent(in) :: n,m
  integer, intent(inout) :: stat
  real(kind=8), intent(in) :: xold(1:n),prm(1:m),fold,stpmax,tolx 
  real(kind=8), intent(inout) :: g(1:n),p(1:n),x(1:n)
  integer :: i,itr,MAXITR 
  real(kind=8) :: Fvec(1:n),ALF,rhs1,rhs2,alam2,a,b,f2,disc
  real(kind=8) :: rtmp,rtmp2,slope,alam,alamin,f,test
  real(kind=8) :: ddot
  parameter(ALF=1.d-4,MAXITR=100)

  stat = 2

  !get initial stepsize
  rtmp = ddot(n,p(1:n),1,p(1:n),1)
  rtmp = sqrt(rtmp)
 
  !scale if attempted step is too large
  if (rtmp .gt. stpmax) then
    p(1:n) = (/ (p(i)*stpmax/rtmp, i=1,n) /)
  end if 

  !get slope
  slope = ddot(n,g(1:n),1,p(1:n),1)
  if (slope .ge. 0.d0) then
    stat = 3
    write(*,*) "ERROR ERROR ERROR"
    write(*,*) "lnsrch : There is roundoff error in lnsrch slope"
    return
  end if
 
  !determine initial test step
  test = 0.d0
  do i=1,n
    rtmp = abs(p(i))/max(abs(xold(i)),1.d0)
    if (rtmp .gt. test) test = rtmp
  end do  
  alamin = tolx/test 
  alam = 1.d0

  !loop 
  do itr=1,MAXITR
!    write(*,*) "itr is",itr
    x = xold
    call daxpy(n,alam,p(1:n),1,x(1:n),1) !x = λ*p + xold
    call fnc(n,m,x(1:n),prm(1:m),Fvec(1:n))
    f = fmin(n,Fvec(1:n))
    
    !delta x convergence
    if (alam .lt. alamin) then
!      if (itr .eq. 1) write(*,*) "Took NR"
      x = xold
      stat = 1
      return

    !sufficient function decrease
    else if (f .le. fold+ALF*alam*slope) then
      stat = 0
!      if (itr .eq. 1) write(*,*) "adapted"
      return

    !we need to backtrack
    else
      if (alam .eq. 1.d0) then    !first backtrack
        rtmp = -1.d0*slope/(2.d0*(f-fold-slope))

      else
        rhs1 = f - fold - alam*slope
        rhs2 = f2 - fold - alam2*slope
        a = (rhs1/alam**2.d0 - rhs2/alam2**2.d0)/(alam-alam2)
        b = (-1.d0*alam2*rhs1/alam**2.d0 + alam*rhs2/alam2**2.d0)/(alam-alam2)

        if (a .eq. 0.d0) then
          rtmp = -1.d0*slope/(2.d0*b)

        else
          disc = b*b-3.d0*a*slope

          if (disc .lt. 0.d0) then
            rtmp = 0.5d0*alam

          else if (b .le. 0.d0) then
            rtmp = (-1.d0*b + sqrt(disc))/(3.d0*a)

          else
            rtmp = -1.d0*slope/(b + sqrt(disc))
          end if !disc val if
        end if !a zero if

        if (rtmp .gt. 0.5d0*alam) rtmp = 0.5d0*alam !λ < 0.5λ1
      end if !backtrack iter if
    end if !backtrack if
    alam2 = alam
    f2 = f
    alam = max(rtmp, 0.1d0*alam) !λ > 0.1λ1
  end do 
  write(*,*) "BACKTRACK DID NOT CONVERGE IN its: ",MAXITR

end subroutine lnsrch
!---------------------------------------------------------------------
