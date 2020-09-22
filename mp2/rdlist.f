      program rdlist
      implicit integer (a-z)
      integer nvrt,nocc,abint
      double precision x,DNRM2,zcore,abenergy,pairenergy,ddot,aaenergy
      double precision pairab,pairaa
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
C     list 14 : beta integrals stored  <ab|ij> for all a<b, i<j triangular 
C     list 16 : RHF integrals stored <AB|IJ> for all A,B and I,J unrestricted
C     list 44 : alpha amplitudes stored t_IJ_AB for all A<B, I<J triangular
C     list 45 : beta  amplitudes stored t_ij_ab for all a<b, i<j triangular 
C     list 46 : RHF   amplitudes stored t_IJ_AB for all A,B and I,J unrestricted
      listint=16
      listamp=46
      write(6,*) "number of irrep",nirrep
      abenergy=0.0d0
      do irrep=1,nirrep
       numdis=moiods(irrep,listint)
       dissiz=moiosz(irrep,listint)
       nocc = nint(sqrt(1.d0*numdis))
       nvrt = nint(sqrt(1.d0*dissiz))
       write(6,*) "numdis",numdis,"nocc",nocc
       write(6,*) "dissiz",dissiz,"nvrt",nvrt
       iint=i0
       iamp=iint+dissiz
       itemp=iamp+dissiz
c read <ab||ij>
c
c mp2 = t(ab,ij) * [2*<ab|ij> - <ba|ij> ]
c
       do ijblock=1,numdis
        call getlst(zcore(iint),ijblock,1,1,irrep,listint)
        call getlst(zcore(iamp),ijblock,1,1,irrep,listamp)
        call transp(zcore(iint),zcore(itemp),nvrt,nvrt) 
        call dscal(dissiz,2.0D0,zcore(iint),1)
        call daxpy(dissiz,-1.0d0,zcore(itemp),1,zcore(iint),1)
        pairenergy=ddot(dissiz,zcore(iint),1,zcore(iamp),1)
        write(6,'(i5,f20.10)')ijblock,pairenergy
        abenergy=abenergy+pairenergy
       enddo
      enddo
      write(6,*)' total ab energy is ',abenergy
C
C
C     JAMES' STUFF
      write(6,*) 'Testing seperate AA,AB'
      aaenergy = 0.d0
      abenergy = 0.d0
      do irrep=1,nirrep 
        numdis = moiods(irrep,listint)
        dissiz = moiosz(irrep,listint)
        nocc = nint(sqrt(1.d0*numdis))
        nvrt = nint(sqrt(1.d0*dissiz))
        write(6,*) "numdis",numdis,"nocc",nocc
        write(6,*) "dissiz",dissiz,"nvrt",nvrt
        iint = i0
        iamp = iint+dissiz
        itemp = iamp+dissiz
    
        ijblock = 0
        do j=1,nocc
          do i=1,nocc
            ijblock = ijblock + 1
            call getlst(zcore(iint),ijblock,1,1,irrep,listint)
            call getlst(zcore(iamp),ijblock,1,1,irrep,listamp)
            call transp(zcore(iint),zcore(itemp),nvrt,nvrt)

            pairenergy = 0.d0
            pairaa = 0.d0
            pairab = 0.d0 

            abint = 0
            do b=1,nvrt  
              do a=1,nvrt
C               AA <= 2*<ab||ij>                  
C               AB <=  <ab|ij>
                pairaa=pairaa+(zcore(iint+abint)-zcore(itemp+abint))*zcore(iamp+abint)
                pairab=pairab+zcore(iint+abint)*zcore(iamp+abint)
                abint = abint + 1  
              end do !loop over b 
            end do !loop over a

            pairenergy=pairaa+pairab
            aaenergy=aaenergy+pairaa
            abenergy=abenergy+pairab
            write(6,'(3(i3),3(f20.10))')ijblock,i,j,pairaa,pairab,pairenergy

          end do !loop over i
        end do !loop over j

      end do !loop over irrep 
      write(6,'(A,F20.10)') "aa  energy",aaenergy
      write(6,'(A,F20.10)') "ab  energy",abenergy
      write(6,'(A,F20.10)') "tot energy",aaenergy+abenergy
      stop
      end 

