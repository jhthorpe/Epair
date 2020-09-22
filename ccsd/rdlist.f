      program rdlist
      implicit integer (a-z)
      integer nvrt,nocc,ab,ia,ib,ja,jb
      double precision x,DNRM2,zcore,abenergy,pairenergy,ddot,aaenergy
      double precision pairab,pairaa,pairaat1,pairaat2,pairabt1,pairabt2
      double precision aat1energy,aat2energy,abt1energy,abt2energy,tote
      common // zcore(1)
      common /machsp/ iintln,ifltln,iintfp,ialone,ibitwd
      common /syminf/ nstart,nirrep,irreps(255,2),dirprd(8,8)
      common /sympop/ irpdpd(8,22),ISYTYP(2,500),id(18)
      common /istart/ i0
      common /iopos/ maxcor,ichcsz,ioff(2),ijunk
      common /lists/ moio(10,500),moiowd(10,500),moiosz(10,500),
     &               moiods(10,500),moiofl(10,500)
C
      call crapsi(zcore,iuhf,0)
C
C     list 14 : alpha integrals stored <AB|IJ> for all A<B, I<J triangular 
C     list 15 : beta integrals stored  <ab|ij> for all a<b, i<j triangular 
C     list 16 : RHF integrals stored <AB|IJ> for all A,B and I,J unrestricted
C     list 44 : alpha amplitudes stored t_IJ_AB for all A<B, I<J triangular
C     list 45 : beta  amplitudes stored t_ij_ab for all a<b, i<j triangular 
C     list 46 : RHF   amplitudes stored t_IJ_AB for all A,B and I,J unrestricted
C     list 90 : RHF   amplitudes stored t_I_A for all A,I virt by occ
      listint=16
      listt2=46
      listt1=90 
      aaenergy=0.d0
      abenergy=0.d0
      write(6,*) "number of irrep",nirrep
      do irrep=1,nirrep
       numij=moiods(irrep,listint)
       numab=moiosz(irrep,listt2)
       nocc = nint(sqrt(1.d0*numij))
       nvrt = nint(sqrt(1.d0*numab))
       write(6,*) "num occ",nocc,"num vrt",nvrt
       iint=i0
       t2=iint+numab
       tmp=t2+numab
       tau=tmp+numab
       t1=tau+numab
      
       !read in amplitudes
       call getlst(zcore(t1),1,1,1,1,listt1) 
C       write(6,*) "t1 amplidues"
C       write(6,*) "i a t_i^a"
C       idx=t1
C       do i=1,nocc
C        do a=1,nvrt
C          ia = nvrt*(i-1) + a-1 
C          write(6,'(i2,1x,i2,1x,f8.5,1x,f8.5)') i,a+nocc,zcore(idx),zcore(t1+ia)
C          idx=idx+1
C        enddo 
C
C       enddo
C
C       write(*,*) "i,j,a,b,tijab,tia,tjb"
       ijblock = 0
       aaenergy = 0.d0
       abenergy = 0.d0
       aat1energy = 0.d0
       aat2energy = 0.d0
       abt1energy = 0.d0
       abt2energy = 0.d0
       tote=0.d0
       write(*,*) "ij,i,j,AA,AB,tot"
       do j=1,nocc
        do i=1,nocc
         ijblock=ijblock+1
         call getlst(zcore(iint),ijblock,1,1,irrep,listint)
         call transp(zcore(iint),zcore(tmp),nvrt,nvrt)
         call getlst(zcore(t2),ijblock,1,1,irrep,listt2)
         ab=0
         pairaa=0.d0
         pairab=0.d0
         pairaat1=0.d0
         pairaat2=0.d0
         pairabt1=0.d0
         pairabt2=0.d0
C         write(6,*) "i,j,a,b,t2"
C         write(*,*) "i,j,a,b,tijab,tia,tjb"
         do b=1,nvrt
          do a=1,nvrt
           ia = nvrt*(i-1)+ a-1
           ib = nvrt*(i-1)+ b-1
           jb = nvrt*(j-1)+ b-1
           ja = nvrt*(j-1)+ a-1
C           write(6,'(4(1x,i2),3(1x,f8.5))') i,j,a+nocc,b+nocc,zcore(t2+ab),zcore(t1+ia),zcore(t1+jb)

            pairabt2=pairabt2+zcore(iint+ab)*zcore(t2+ab)
            pairabt1=pairabt1+zcore(iint+ab)*zcore(t1+ia)*zcore(t1+jb)
            pairaat1=pairaat1+(zcore(iint+ab)-zcore(tmp+ab))*(zcore(t1+ia)*zcore(t1+jb))
!            pairaat1=pairaat1+(zcore(iint+ab))*(zcore(t1+ia)*zcore(t1+jb)-zcore(t1+ib)*zcore(t1+ja))
            pairaat2=pairaat2+(zcore(iint+ab)-zcore(tmp+ab))*zcore(t2+ab)
          
           ab=ab+1
          enddo !a loop
         enddo !b loop

         pairaa=pairaat1+pairaat2
         pairab=pairabt1+pairabt2
         pairenergy=pairaa+pairab
         aaenergy=aaenergy+pairaa
         abenergy=abenergy+pairab
         aat1energy=aat1energy+pairaat1
         aat2energy=aat2energy+pairaat2
         abt1energy=abt1energy+pairabt1
         abt2energy=abt2energy+pairabt2
C         write(6,'(3(i3),5(f20.10))')ijblock,i,j,pairaat1,pairaat2,pairabt1,pairabt2,pairenergy
         write(6,'(3(i3),5(f20.10))') ijblock,i,j,pairaat1+pairaat2,pairabt1+pairabt2,pairenergy

        enddo !i loop
       enddo !j loop

      enddo !irrep loop 
      write(*,*) "RESULTS"
      write(6,'(A,F20.10)') "aa  energy",aaenergy
      write(6,'(A,F20.10)') "ab  energy",abenergy
      write(6,'(A,F20.10)') "aa  t1 energy",aat1energy
      write(6,'(A,F20.10)') "ab  t1 energy",abt1energy
      write(6,'(A,F20.10)') "aa  t2 energy",aat2energy
      write(6,'(A,F20.10)') "ab  t2 energy",abt2energy
      write(6,'(A,F20.10)') "tot energy",aaenergy+abenergy
     

      stop
      end
