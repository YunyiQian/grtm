
        subroutine kmaxcalculate(o,kmax,nf,nfmax)
  
        include "green_com.inc"
        real*8  pmin,dzmax,kmax
        integer nf,nfmax

        if(lo.eq.ls-1) then
          pmin=dmax1(cdabs(o/vs(ls)),cdabs(o/vs(ls-1)))
          dzmax=dmin1(zs-z(ls-1),z(ls-1)-z0)
        else if(lo.lt.ls-1) then
          pmin=dmax1(cdabs(o/vs(ls)),cdabs(o/vs(ls-1)))
          dzmax=dmin1(zs-z(ls-1),z(ls-1)-z0)
          do lay=lo,ls-2
             if(pmin.lt.cdabs(o/vs(lay+1)))  pmin=cdabs(o/vs(lay+1))
             if(dzmax.gt.z(lay+1)-z(lay))   dzmax=z(lay+1)-z(lay)
          end do
        else if(lo.eq.ls) then
          pmin=cdabs(o/vs(ls))
          dzmax=dabs(zs-z0)
        else if(lo.eq.ls+1) then
          pmin=dmax1(cdabs(o/vs(ls)),cdabs(o/vs(ls+1)))
          dzmax=dmin1(z(ls)-zs,z0-z(ls))
        else if(lo.gt.ls+1) then
          pmin=dmax1(cdabs(o/vs(ls)),cdabs(o/vs(ls+1)))
          dzmax=dmin1(z(ls)-zs,z0-z(ls))
          do lay=ls+1,lo-1
             if(pmin.lt.cdabs(o/vs(lay)))  pmin=cdabs(o/vs(lay))
             if(dzmax.gt.z(lay)-z(lay-1))  dzmax=z(lay)-z(lay-1)
          end do
        end if
        if(dzmax.lt.0.2d0) then
          if(nf.lt.10) then
             kmax=dmin1(pmin+2.0,20.*pmin)
          else
             kmax=dmin1(pmin+2.0,2.*pmin)
          end if
        else      
          kmax=dsqrt((3.0/dzmax)**2+pmin**2)    
        end if
        if(kmax-pmin.lt.0.5d0) kmax=pmin+0.5d0

        return
        end 

!c----------------------------------------------------
!c       function to compute the layer-number of an  l
!c       arbitrary depth varible "z"                 l
!c       ---------------------------------------------

        function layernumb(z0,n,z)
!c       ........
        real*8    z0,z(0:n)
        integer layernumb, n, i
        do i = n-1,0,-1
           if (z0.gt.z(i)) then
              layernumb = i+1
              return
           end if
        end do
        layernumb=1
!c       print*,'Error: info@layernumb'  

        return
        end
           
!c       .........................................
        function amax(n,a)
        real*8     a(n), amax
        integer  n, i
!c       ..............
        amax=a(1)
        do i=1,n
           if (a(i).gt.amax) amax=a(i)
        end do
        return
        end 
           

!c       .........................................
        function amin(n,a)
        real*8     a(n), amin
        integer  n, i
!c       ..............
        amin=a(1)
        do i=1,n
           if (a(i).lt.amin) amin=a(i)
        end do
        return
        end 

        
!c---------------------------------------------------------
!c       function to compute the source-time function     l
!c       -------------------------------------------      l
!c       
        complex*16 function sw(w)

        include "green_com.inc"
        complex*16  w, c
        real*8      wc
!        common      /B_shift/ shift
        
!	fc=2.0
        wc = 2.*pi*fc

!c   CASE 1. Ricker wavelet: 
!c       r(t) = (u**2 - .5)*sqrt(2.)*.5*exp(-u**2/4.),
!c       where, u  = 2*sqrt(6.)*(t-to)/tb, 
!c              to -- arrival time,
!c              tb -- width of wavelet.
!c       Its spectral is:
!c       R(w) =  (w/wc)**2 * exp(-(w/wc)**2),
!c        wc  -- central frequency.
        if ( Type_STF.eq.'Ricker' ) then
           c  = (w/wc)**2
           sw = c * cdexp(-c) / wc
        return
        end if

!c   CASE 2. Exponential function:
!c       e(t) = 0.0                      for t < 0;
!c       e(t) = (1. - exp( -t/tou))      for t >= 0.
!c       tou -- rise time.
!c       Its spectral is:
!c       E(w) = 1./(i*w*(1+i*w*tou))/sqrt(2.*pi).
        if ( Type_STF.eq.'Exp' ) then
           sw = 1./(aj*w*(aj*w*tou-1))/sqrt(2.*pi)
        return
        end if

!c   CASE 3. Smooth step-function:
!c       sms(t) = (2/pi)*atan(t/tou).
!c       Its spectral is:
!c       SMS(w) = -i*sqrt(pi/2)*exp(-tou*abs(w))/w.
        if ( Type_STF.eq.'SSTEP' ) then
           sw = -aj*sqrt(pi/2.)*exp(-tou*cdabs(w))/w
        return
        end if

!c   CASE 4. Ramp function:
!c       ra(t) = 0.0                     for t < 0;
!c       ra(t) = t/tou                   for 0 < t < tou;
!c       ra(t) = 1.0                     for t > tou.
!c       Its spectal is:
!c       RA(w) = (exp(-i*w*tou)-1.)/(tou*w**2).
        if ( Type_STF.eq.'Ramp' ) then
           sw = (cdexp(-aj*w*tou) - 1.)/(tou*w*w)      !!!!!!!!!!!!
!c           sw = -ai*(cdexp(-ai*w*tou) - 1.)/(tou*w)
        return                                                          
        end if

!c   CASE 5. Bouchon function:
!c     ra(t) = tanh(t/to);
!c       Its spectal is:
!c       RA(w) = -i*tou*pi*cosech(w*pi*tou/2).
        if ( Type_STF.eq.'Bouchon' ) then
           sw =-aj*tou*pi*2./(cdexp(w*pi*tou/2)-cdexp(-w*pi*tou/2))
        return
        end if

        if ( Type_STF.eq.'Pitaka' ) then
!c          ur_z(i)=-ai*to*pi*cdexp(-ai*w*shift*to)
!c       &            /(cdexp(w*pi*to/2)-cdexp(-w*pi*to/2))*dw*mt/(2*pi)                
          sw =pi2**2*(cdexp(-aj*w*tou) - 1.)/(tou*w*w)
     &        /(pi2**2-w**2*tou*tou)*exp(-aj*w*shift)
        return
        end if
        
        if ( Type_STF.eq.'B(shift)' ) then
!c          ur_z(i)=-ai*to*pi*cdexp(-ai*w*shift*to)
!c       &            /(cdexp(w*pi*to/2)-cdexp(-w*pi*to/2))*dw*mt/(2*pi)                
                sw =-aj*tou*pi*cdexp(-aj*w*shift)
     &        /(cdexp(w*pi*tou/2)-cdexp(-w*pi*tou/2))
        return
        end if
!c   CASE 6. A pulse
        if ( Type_STF.eq.'Green' ) then
           sw = cmplx(1., 0.)
        return
        end if

!c   CASE 7. Quadratic ramp function:
!c     ra(t)=0.0     for t<0 and t>.6;
!c     ra(t)=t/.2     for 0<t<.2;
!c     ra(t)=1       for .2<t<.4;
!c     ra(t)=-t/.2+.3  for .4<t<.6.
!c    Its spectral is:
!c    RA(w)=(exp(-ai*4.*w)-1.)(1.-exp(-ai*2.*w))/(2.*w*w).
        if ( Type_STF.eq.'Trapezoid' ) then
          sw = aj*(1.-cdexp(-aj*0.4*w))*(1.-cdexp(-aj*0.2*w))
     &           /(0.2*w*w*w)
        return
        end if

!c   CASE 8. Heaviside step function:
!c     ra(t)=0.0     for t<0 ;
!c     ra(t)=1.0     for t>0 .
!c    Its spectral is:
!c    RA(w)=1./(ai*w).
        if ( Type_STF.eq.'Heaviside' ) then
          sw = 1./(aj*w)
        return
        end if

!---------------------------------------------------------------------
!c   CASE 9. Triangle ramp function:
!c     ra(t)=0.0     for t<0 and t>tou;
!c     ra(t)=2at/tou     for 0<=t<tou/2;
!c     ra(t)=-2at/tou+2a  for tou/2<=t<=tou.
!c     a=2/tou
!c    Its spectral is:
!c    RA(w)=-2*(exp(-ai*w*tou/2.)-1.)**2/(tou*w*w).
        if ( Type_STF.eq.'Triangle' ) then
          sw = -4.*(cdexp(-aj*w*tou/2.)-1.)*(cdexp(-aj*w*tou/2.)-1.)
     &           /(tou*tou*w*w)
          return
        end if

!---------------------------------------------------------------------
!c   CASE 10. Half cos function:
!c     ra(t)=0.0     for t<0 and t>tou;
!c     ra(t)=(1-cos(2*pi)/tou)/tou     for 0<=t<tou;
!c   Its spectral is:
!c   RA(w)=-((2*pi)^2/(aj*(w*tou)^3))*((1-exp(-aj*w*tou))/(1-(2*pi)^2/(w*tou)^2))
	  if ( Type_STF.eq.'COS' ) then
		sw = -(pi2*pi2)/(aj*(w*tou)**3)*(1-cdexp(-aj*w*tou))
     &		/(1-pi2*pi2/(w*tou)**2)
	    return
	  end if

!---------------------------------------------------------------------

!c   If smod dose not equal to 1, 2, 3, 4 or 5, it means the input 
!c   of smod is wrong.
        print *, 'Warrning: info@sw -- Wrong "smod" input,
     &            please re-input smod (1, 2, 3, or 4). '
        stop

        end     

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	real*8 function window(f)

	character*80 WinSwitch,WinType
	real*8       f,a,pi,f1,f2,f3,f4
      common      /filter/ WinSwitch, WinType
	common      /f1234/ f1,f2,f3,f4

      pi=4d0*atan(1d0)

!c   CASE 1 : All_pass filter :
      if (WinSwitch.eq.'OFF') then
	   window=1.
	   return
	end if

!c   CASE 2 : Hamming window :
	if(WinSwitch.eq.'ON'.and.WinType.eq.'Hamming') then

		if (f.lt.f1.or.f.gt.f4) then
			window=0d0
	        return
	    end if

		if (f.le.f3.and.f.ge.f2) then
			window=1d0
	        return
	    end	if

	    if (f.ge.f1.and.f.lt.f2.and.f1.ne.f2) then
			a=(f-f1)/(f2-f1)*pi
			window=(0.46d0-0.46d0*cos(a))/0.92d0
	        return
	    end	 if

	    if (f.gt.f3.and.f.le.f4.and.f3.ne.f4) then
			a=(f-f3)/(f4-f3)*pi
			window=(0.46d0+0.46d0*cos(a))/0.92d0
	        return
		end if

	end if

	end
