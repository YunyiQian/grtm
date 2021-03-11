
      subroutine green_spectra_r(nf1, nf2, integf)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Note: This is the core subroutine of the main program. The purpose c
!c         of this subroutine is to calculate the Green spectrum for    c
!c         the array-case.                                              c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      use         MSIMSL       
      include     "green_com.inc"

      integer    :: nf
      complex*16 :: integf(n/2,10)
      real*8     :: kmin, kmax
      complex*16 :: sums(10)
      complex*16 :: om,cs,cp,omax,otemp
!      integer    :: min,sec
      real*8     ::  time_begin, time_end
      real*8     :: f0
      real*8 :: phi

!      open(2, file='CON', carriagecontrol='FORTRAN')
!c      call timdy(hour1,minute1,second1)
      call cpu_time(time_begin)

      f0 = fmax
      do nf=nf1,nf2
         o   = pi2*df*(nf-1) - aj*oi       ! complex*16 frequency
         om  = cmplx(pi2)                  ! referenced frequency
         do i=1,nly                        ! complex*16 velocities
! according to Aki & Richard, 1980, P182
!            cs = (1.0+cdlog(o/pi2/f0)/(pi*Qs(i))-aj/(2.0*Qs(i)))
!            cp = (1.0+cdlog(o/pi2/f0)/(pi*Qp(i))-aj/(2.0*Qp(i)))
            cs = (1.0+cdlog(o/pi2)/(pi*Qs(i))+aj/(2.0*Qs(i)))
            cp = (1.0+cdlog(o/pi2)/(pi*Qp(i))+aj/(2.0*Qp(i)))
            vs(i) = vs0(i)*cs
            vp(i) = vp0(i)*cp
         end do           

         if (nf.le.3) then
            otemp = pi2*df*3d0 - aj*oi
            call kmaxcalculate(otemp,kmax,3,nf2)
         else 
            call kmaxcalculate(o,kmax,nf,nf2)
         end if
	   kmax = 1.1*kmax
                     
         call DWIM(o, nf, 0d0, kmax, sums)
         !integf(nf,:) = sums*exp(pi2*aj*shift*dt*o)
         phi = pi2*df*(nf-1)*shift
         integf(nf,:) = sums*cmplx(cos(phi),sin(phi))

         write(2,110) ' ',nf,'/',nf2,' = ',float(nf)*100./float(nf2),   
     &              '% ;   max(k)=',kmax

      end do

!c      call timdy(hour2,minute2,second2)
      call cpu_time(time_end)
      print *, 'Time cost = ',time_end-time_begin, ' seconds'
!      sec=(hour2-hour1)*3600+(minute2-minute1)*60+second2-second1
!      sec = time_end-time_begin
!      min=sec/60
!      sec=sec-min*60

!      if (min.eq.0)   then
!         print*,' *****************************************'
!         write(*,300) ' **    Time cost = ', sec, ' seconds.   **'
!      else  
!         print*,' **************************************************'
!         write(*,310) ' **   Time cost = ',min,' minutes ',sec,
!     &              ' seconds.  **' 
!      end if

110   format('+',2(a,i5),a,f7.2,a,f7.4)      
300   format(1x,a,i3,a)
310   format(1x,a,i5,a,i3,a)

!      close(2)
      return
      end
