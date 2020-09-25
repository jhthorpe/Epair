! E_X = E_inf + a/X^n
program extrap
  use nlslv
  implicit none
  integer :: zmin,zmax,npair,nocc,ij
  integer :: fid,idum,i,j
  real(kind=8) :: rdum,ztol
  integer, dimension(1:3) :: zeta
  integer, dimension(:), allocatable :: statAA,statAB
  real(kind=8), dimension(:), allocatable :: exptAA,fctrAA,ExAA
  real(kind=8), dimension(:), allocatable :: exptAB,fctrAB,ExAB
  real(kind=8), dimension(:,:), allocatable :: AA,AB,errAA,errAB
  real(kind=8) :: var(1:2),prm(1:6),errF(1:2)
  real(kind=8) :: Fjac(1:3,1:3)

  fid = 100
  ztol = 1.d-15

  !read in data
  open(file='Epair',unit=fid,status='old')
  read(fid,*) nocc,zmin,zmax 
  read(fid,*)
  npair = nocc*nocc
  write(*,*) "npair,zmin,zmax",npair,zmin,zmax
  allocate(AA(zmin:zmax,1:npair))
  allocate(AB(zmin:zmax,1:npair))
  allocate(exptAA(1:npair))
  allocate(exptAB(1:npair))
  allocate(fctrAA(1:npair))
  allocate(fctrAB(1:npair))
  allocate(ExAA(1:npair))
  allocate(ExAB(1:npair))
  allocate(errAA(1:3,1:npair))
  allocate(errAB(1:3,1:npair))
  allocate(statAA(1:npair))
  allocate(statAB(1:npair))
  errAA = 0.d0
  errAB = 0.d0
  statAB = 0
  statAA = 0

  !read AA
  do ij=1,npair 
    read(fid,*) idum,AA(zmin:zmax,ij)
  end do  
  read(fid,*)

  !read AB
  do ij=1,npair 
    read(fid,*) idum,AB(zmin:zmax,ij)
  end do  
  close(unit=fid)

  !printing
  write(*,*) "AA pairs"
  call Epair_print(npair,zmin,zmax,AA(zmin:zmax,1:npair))
  write(*,*) "AB pairs"
  call Epair_print(npair,zmin,zmax,AB(zmin:zmax,1:npair))
  write(*,*)
  write(*,*)

  !initial guess
  write(*,*) "enter initial guess for AA exponent"
  read(*,*) rdum
  exptAA = rdum
  write(*,*) "enter 2 zeta to use"
  read(*,*) zeta(1:2) 
  idum = maxval(zeta(1:2))
  do ij=1,npair 
    fctrAA(ij) = (AA(zeta(2),ij) - AA(zeta(1),ij)) / (1.0/zeta(2)**exptAA(ij) - 1.d0/zeta(1)**exptAA(ij))  
    ExAA(ij) = Einf(idum,AA(idum,ij),fctrAA(ij),exptAA(ij)) 
  end do  
  call Extrp_print(npair,ExAA(1:npair),fctrAA(1:npair),exptAA(1:npair),errAA(1:3,1:npair),statAA(1:npair))

  write(*,*) "enter initial guess for AB exponent"
  read(*,*) rdum
  exptAB = rdum
  write(*,*) "enter 2 zeta to use"
  read(*,*) zeta(1:2) 
  idum = maxval(zeta(1:2))
  do ij=1,npair 
    fctrAB(ij) = (AB(zeta(2),ij) - AB(zeta(1),ij))/ (1.d0/zeta(2)**exptAB(ij) - 1.d0/zeta(1)**exptAB(ij))  
    ExAB(ij) = Einf(idum,AB(idum,ij),fctrAB(ij),exptAB(ij)) 
  end do  
  call Extrp_print(npair,ExAB(1:npair),fctrAB(1:npair),exptAB(1:npair),errAB(1:3,1:npair),statAB(1:npair))
  
  !iterative procedure per pair
  ! each iterative procudure needs to determine a new fctr from the 
  ! guessed expt and Ex, otherwise hell breaks loose in the optimization
  !AA
  write(*,*) "Solving for AA block"
  write(*,*) "Enter zeta to use"
  read(*,*) zeta(1:3)
  idum = maxval(zeta(1:3))
  do ij=1,npair
    var = [fctrAA(ij),exptAA(ij)]
    prm = [1.d0*zeta(1),1.d0*zeta(2),1.d0*zeta(3),AA(zeta(1),ij),AA(zeta(2),ij),AA(zeta(3),ij)]
    call gcnr(2,6,var,prm(1:6),errF(1:2),errAA(1:2,ij),1.d-11,1.d-11,100,statAA(ij))
    fctrAA(ij) = var(1)
    exptAA(ij) = var(2)
    
!    call nr2(100,zeta(1:3),[AA(zeta(1),ij),AA(zeta(2),ij),AA(zeta(3),ij)],&
!             fctrAA(ij),exptAA(ij),1.d-7,1.d-7)

    ExAA(ij) = Einf(idum,AA(idum,ij),fctrAA(ij),exptAA(ij))
    errAA(3,ij) = abs(Einf(idum,AA(idum,ij),fctrAA(ij)+errAA(1,ij),exptAA(ij)+errAA(2,ij)) - ExAA(ij))
  end do 
  call Extrp_print(npair,ExAA(1:npair),fctrAA(1:npair),exptAA(1:npair),errAA(1:3,1:npair),statAA(1:npair))

  !AB
  write(*,*) "Solving for AB block"
  write(*,*) "Enter zeta to use"
  read(*,*) zeta(1:3)
  idum = maxval(zeta(1:3))
  do ij=1,npair

    var = [fctrAB(ij),exptAB(ij)]
    prm = [1.d0*zeta(1),1.d0*zeta(2),1.d0*zeta(3),AB(zeta(1),ij),AB(zeta(2),ij),AB(zeta(3),ij)]
    call gcnr(2,6,var(1:2),prm(1:6),errF(1:2),errAB(1:2,ij),1.d-11,1.d-11,100,statAB(ij))
    fctrAB(ij) = var(1)
    exptAB(ij) = var(2)

!    call nr2(100,zeta(1:3),[AB(zeta(1),ij),AB(zeta(2),ij),AB(zeta(3),ij)],&
!             fctrAB(ij),exptAB(ij),1.d-7,1.d-7)

    ExAB(ij) = Einf(idum,AB(idum,ij),fctrAB(ij),exptAB(ij))
    errAB(3,ij) = abs(Einf(idum,AB(idum,ij),fctrAB(ij)+errAB(1,ij),exptAB(ij)+errAB(2,ij)) - ExAB(ij))
  end do 
  call Extrp_print(npair,ExAB(1:npair),fctrAB(1:npair),exptAB(1:npair),errAB(1:3,1:npair),statAB(1:npair))

  deallocate(AA)
  deallocate(AB)
  deallocate(exptAA)
  deallocate(exptAB)
  deallocate(fctrAA)
  deallocate(fctrAB)
  deallocate(ExAA)
  deallocate(ExAB)
  deallocate(errAA)
  deallocate(errAB)
  
contains

real(kind=8) function Einf(zeta,Ezeta,fctr,expt)
  integer, intent(in) :: zeta
  real(kind=8), intent(in) :: Ezeta,expt,fctr
  Einf = Ezeta - fctr/(1.d0*zeta)**expt
end function Einf

subroutine Epair_print(npair,zmin,zmax,Epair)
  integer, intent(in) :: npair,zmin,zmax
  real(kind=8), dimension(zmin:zmax,npair), intent(in) :: Epair
  integer :: ij
  write(*,'(1x,A,999(19x,I1))') "pair",(ij,ij=zmin,zmax)
  write(*,*) "--------------------------------------------------"
  do ij=1,npair
    write(*,'(1x,I4,999(1x,F20.10))') ij,Epair(zmin:zmax,ij)
  end do 
end subroutine Epair_print

subroutine Extrp_print(npair,Ex,fctr,expt,errx,stat)
  integer, intent(in) :: npair
  integer, dimension(1:npair), intent(in) :: stat
  real(kind=8), dimension(1:npair), intent(in) :: Ex,fctr,expt
  real(kind=8), dimension(1:3,npair), intent(in) :: errx
  integer :: ij
  write(*,'(1x,A,15x,A,24x,A,22x,A)') "pair","Einf","factor","exponent"
  write(*,*) "--------------------------------------------------------------------------------------------------"
  do ij=1,npair
    write(*,'(1x,I4,3(1x,F15.10,1x,A,F9.6,A))',advance='no') ij,Ex(ij),"(",errx(3,ij),")",&
            fctr(ij),"(",errx(1,ij),")",expt(ij),"(",errx(2,ij),")"
    if (stat(ij) .eq. 1) then
      write(*,'(1x,A)') "LOCAL"
    else if (stat(ij) .eq. 2) then
      write(*,'(1x,A)') "ERROR"
    else
      write(*,*)
    end if
!    write(*,'(1x,I4,3(1x,F20.10))') ij,Ex(ij),fctr(ij),expt(ij)
    
  end do 
  write(*,*) "--------------------------------------------------------------------------------------------------"
  write(*,'(1x,A,14x,A,15x,A,15x,A)') "   ","Einf","factor","exponent"
  write(*,'(1x,A,1x,3(1x,F20.10))') "avg  ", nzavg(npair,Ex(1:npair)),&
                                             nzavg(npair,fctr(1:npair)),&
                                             nzavg(npair,expt(1:npair))
  write(*,'(1x,A,1x,3(1x,F20.10))') "rmsd ", nzrmsd(npair,Ex(1:npair)),&
                                             nzrmsd(npair,fctr(1:npair)),&
                                             nzrmsd(npair,expt(1:npair))
  write(*,'(1x,A,2x,F20.10)') "total", sum(Ex(1:npair))
end subroutine Extrp_print

real(kind=8) function nzavg(n,x)
  integer, intent(in) :: n
  real(kind=8), dimension(1:n), intent(in) :: x
  integer :: i,nz
  real(kind=8) :: rdum,ztol
  rdum = 0.d0
  ztol = 1.d-15
  nz = 0
  do i=1,n
    if (abs(x(i)) .gt. ztol) then
      rdum = rdum + x(i)
      nz = nz + 1
    end if 
  end do  
  nzavg = rdum/(1.d0*nz)
end function nzavg

real(kind=8) function nzrmsd(n,x)
  integer, intent(in) :: n
  real(kind=8), dimension(1:n), intent(in) :: x
  integer :: i,nz
  real(kind=8) :: rdum,ztol,avg
  rdum = 0.d0
  ztol = 1.d-15
  nz = 0
  avg = nzavg(n,x(1:n))
  do i=1,n
    if (abs(x(i)) .gt. ztol) then
      rdum = rdum + (x(i)-avg)**2.d0
      nz = nz + 1
    end if
  end do 
  nzrmsd = sqrt(rdum/(1.d0*nz))
end function nzrmsd

subroutine fvec2(zeta,Ezeta,fctr,expt,fv)
  !generates fv for nr3
  implicit none
  integer, intent(in) :: zeta(1:3)
  real(kind=8),intent(in) :: fctr,expt,Ezeta(1:3)
  real(kind=8), intent(inout) :: fv(1:2)
  integer :: j
  fv(1) = Ezeta(1) - Ezeta(2) - fctr*(1.d0/zeta(1)**(expt) - 1.d0/zeta(2)**(expt)) 
  fv(2) = Ezeta(2) - Ezeta(3) - fctr*(1.d0/zeta(2)**(expt) - 1.d0/zeta(3)**(expt))
end subroutine fvec2

subroutine fvec3(zeta,Ezeta,Ex,fctr,expt,fv)
  !generates fv for nr3
  implicit none
  integer, intent(in) :: zeta(1:3)
  real(kind=8),intent(in) :: Ex,fctr,expt,Ezeta(1:3)
  real(kind=8), intent(inout) :: fv(1:3)
  integer :: j
  fv(1) = Ex - Ezeta(1) + fctr/(1.d0*zeta(1))**(expt)
  fv(2) = Ex - Ezeta(2) + fctr/(1.d0*zeta(2))**(expt)
  fv(3) = Ex - Ezeta(3) + fctr/(1.d0*zeta(3))**(expt)
  !fv(1:3) = (/ (Ex - Ezeta(j) + fcrt/(1.d0*zeta(j))**(expt)), j=1,3/)
end subroutine fvec3

subroutine anjac2(zeta,fctr,expt,Fjac)
  !analytical partial-deriv of fnc3 (assuming Ex,fctr,and expt are seperable)
  !f1 = EX - EY + a(1/Y^n - 1/X^n)
  !f2 = EY - EZ + a(1/Z^n - 1/Y^n)
  implicit none
  integer, intent(in) :: zeta(1:3)
  real(kind=8), intent(in) :: fctr,expt
  real(kind=8), dimension(1:2,1:2), intent(inout) :: Fjac
  integer :: i
  Fjac(1,1) = (1.d0/zeta(2)**expt - 1.d0/zeta(1)**expt)
  Fjac(1,2) = fctr*(log(1.d0*zeta(1))/zeta(1)**expt - log(1.d0*zeta(2))/zeta(2)**expt) 
  Fjac(2,1) = (1.d0/zeta(3)**expt - 1.d0/zeta(2)**expt)
  Fjac(2,2) = fctr*(log(1.d0*zeta(2))/zeta(2)**expt - log(1.d0*zeta(3))/zeta(3)**expt)
end subroutine anjac2

subroutine anjac3(zeta,fctr,expt,Fjac)
  !analytical partial-deriv of fnc3 (assuming Ex,fctr,and expt are seperable)
  ! F = Ex - Ezeta + a/zeta^n
  implicit none
  integer, intent(in) :: zeta(1:3)
  real(kind=8), intent(in) :: fctr,expt
  real(kind=8), dimension(1:3,1:3), intent(inout) :: Fjac
  integer :: i
  Fjac(1,1) = 1.d0 !dF/dEx = 1
  Fjac(1,2) = (1.d0*zeta(1))**(-1.d0*expt) !dF/da
  Fjac(1,3) = -1.d0*fctr*log(1.d0*zeta(1))/(zeta(1)**expt)!dF/dn 
  Fjac(2,1) = 1.d0 !dF/dEx = 1
  Fjac(2,2) = (1.d0*zeta(2))**(-1.d0*expt)
  Fjac(2,3) = -1.d0*fctr*log(1.d0*zeta(2))/(zeta(2)**expt)!dF/dn 
  Fjac(3,1) = 1.d0 !dF/dEx = 1
  Fjac(3,1) = (1.d0*zeta(3))**(-1.d0*expt)
  Fjac(3,3) = -1.d0*fctr*log(1.d0*zeta(3))/(zeta(3)**expt)!dF/dn 
end subroutine anjac3

subroutine nr3(ntrial,zeta,Ezeta,Ex,fctr,expt,tolx,tolf)
  !uses newton-rapsom to solve NL eqs for Ex,fctr,expt
  implicit none
  integer, intent(in) :: ntrial,zeta(1:3)
  real(kind=8), intent(in) :: Ezeta(1:3),tolx,tolf
  real(kind=8), intent(inout) :: Ex,fctr,expt
  integer :: i,k,p(1:3),stat
  real(kind=8) :: errf,errx,Fjac(1:3,1:3),fv(1:3)
  real(kind=8) :: x(1:3),b(1:3)
  write(*,*) "hi hi hi"
  write(*,*) "input values:"
  write(*,*) "zeta:",zeta(1:3)
  write(*,*) "Ezeta:",Ezeta(1:3)
  write(*,*) "Ex,fctr,expt",Ex,fctr,expt
  do k=1,ntrial
    call fvec3(zeta(1:3),Ezeta(1:3),Ex,fctr,expt,fv(1:3))
    call anjac3(zeta(1:3),fctr,expt,Fjac(1:3,1:3))
    errf = 0.d0
    do i=1,3
      errf = errf+abs(fv(i))
    end do 
   if(errf .le. tolf .and. k .gt. 1) return
    x(1:3) = [Ex,fctr,expt]
    b = -1.d0*fv
!    write(*,*) "bbefore",b
!    write(*,*) "xbefore",x
    call dgesv(3,1,Fjac(1:3,1:3),3,p(1:3),b(1:3),3,stat)
!    write(*,*) "stat is",stat
    Ex = Ex + b(1) 
    fctr = fctr + b(2) 
    expt = expt + b(3) 
    errx = 0.d0
    do i=1,3
      errx = errx + abs(b(i))
    end do 
!    write(*,*) "b =", b(1:3)
!    write(*,*) 
    write(*,*) "iter :", k
    write(*,*) "x =", Ex,fctr,expt
    write(*,*) "errf, errx",errf,errx
    if (errx .le. tolx .and. k .gt. 1) return
  end do 
end subroutine nr3

subroutine nr2(ntrial,zeta,Ezeta,fctr,expt,tolx,tolf)
  !uses newton-rapsom to solve NL eqs for fctr,expt
  implicit none
  integer, intent(in) :: ntrial,zeta(1:3)
  real(kind=8), intent(in) :: Ezeta(1:3),tolx,tolf
  real(kind=8), intent(inout) :: fctr,expt
  integer :: i,k,p(1:3),stat
  real(kind=8) :: errf,errx,Fjac(1:2,1:2),fv(1:2)
  real(kind=8) :: x(1:2),b(1:2)
  write(*,*) "hi hi hi"
  write(*,*) "input values:"
  write(*,*) "zeta:",zeta(1:3)
  write(*,*) "Ezeta:",Ezeta(1:3)
  write(*,*) "fctr,expt",fctr,expt
  do k=1,ntrial
    write(*,*) "-------------------------------------"
    write(*,*) "NR2 iter",k
!    call fvec2(zeta(1:3),Ezeta(1:3),fctr,expt,fv(1:2))
!    call anjac2(zeta(1:3),fctr,expt,Fjac(1:2,1:2))
    call fnc(2,6,[fctr,expt],[1.d0*zeta(1),1.d0*zeta(2),1.d0*zeta(3),Ezeta(1),Ezeta(2),Ezeta(3)],&
             fv(1:2))
    call jac(2,6,[fctr,expt],[1.d0*zeta(1),1.d0*zeta(2),1.d0*zeta(3),Ezeta(1),Ezeta(2),Ezeta(3)],&
             Fjac(1:2,1:2))
    errf = 0.d0
    do i=1,2
      errf = errf+abs(fv(i))
    end do 
    if(errf .le. tolf .and. k .gt. 1) return
    x(1:2) = [fctr,expt]
    b = -1.d0*fv
    write(*,*) "x old =",x
    call dgesv(2,1,Fjac(1:2,1:2),2,p(1:2),b(1:2),2,stat)
!    write(*,*) "stat is",stat
    fctr = fctr + b(1) 
    expt = expt + b(2) 
    errx = 0.d0
    do i=1,2
      errx = errx + abs(b(i))
    end do 
!    write(*,*) "b =", b(1:3)
!    write(*,*) 
    write(*,*) "x new =", fctr,expt
    if (errx .le. tolx .and. k .gt. 1) return
  end do 
end subroutine nr2

end program extrap
