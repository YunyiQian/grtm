!c
        subroutine green_basic

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Note: The purpose of this subroutine is to assign values to  c
!c         some basic parameters, and calclulate values for some  c
!c         other basic parameters.                                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        include    "green_com.inc"	
	  real*8     amax, amin, vave
        real*8 maxr0

!c  SETTING BASIC PARAMETERS:
!c  -------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        m1   = 2**(m-1)         ! Half of sampling points
        mt   = 2*m1             ! Total sampling points
	  Twin = dt*mt
        print*, 'Twin',Twin
	  fmax = 1d0/(2d0*dt)	      ! Nyquist frequency
	  !df   = fmax/(m1-1)	  ! frequency interval
      df = 1/(mt*dt)
	  !fmax = 1d0/dt
	  !df = fmax/mt
 	  !oi   = pi/Twin   ! imaginary frequency
      oi = 2/Twin
	  
	  vmax = amax(nly,vp0)	  ! max velocity
	  vmin = amin(nly,vs0)	  ! min velocity
        vave = 0.5*(vmax+vmin)
        maxr0 = amax(ndist,r0)
        L = (vmax*Twin+maxr0+vmax/vave*dsqrt(maxr0*maxr0+(zs-z0)**2))+100        
        L =2*L
	  dk   = pi2/L         ! k-integral step

	  print*,'!!!!!!!!!!!!!!!!!!! Basic Parameters !!!!!!!!!!!!!!!!!!'
	  write(*,*)  '              Time window width : Twin =',Twin
	  write(*,*)  '     Sample points(time domain) : mt   =',mt
	  write(*,*)  '       number distance : R    =',ndist
	  write(*,*)  '         Increment of time(sec) : dt   =',dt
	  write(*,*)  '           Cutoff frequency(Hz) : fmax =',fmax
	  write(*,*)  '     Increment of frequency(Hz) : df   =',df
	  write(*,*)  'The period structure length(km) : L    =',L
	  write(*,*)  ' Increment of wavenumber(km^-1) : dk   =',dk
	  print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

100     format(1x,a,f12.4)
200     format(1x,a,i5)

!c	  z(nly)=z(nly-1)+5.0

        do i = 1, nly 
           mu(i) = rho(i)*vs0(i)**2
        end do

!c       Identity matrix:
        unit(1,1) = 1d0
        unit(1,2) = 0d0
        unit(2,1) = 0d0
        unit(2,2) = 1d0

        RETURN
        END


