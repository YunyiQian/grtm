
        subroutine green_input

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Note: The purpose of this subroutine is to input parameters c
!c         from the existed data files. So, before run the main  c
!c         program, the input data files must be prepared in     c
!c         appropriate forms.                                    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        include    "green_com.inc"
        character  list*80
        integer    layernumb
        real*8     mus,x0,y0,xs,ys
!        common     /B_shift/ shift

!c   Reading media parameters:

 101    format ( I10   )
! 102    format ( f10.8 )
 103    format ( f10.1 )
 104    format ( f10.4 )
 105    format ( a10   )
 106    format ( a90   )

        open(10,file='grt.dat',form='formatted',status='old')
        read(10, 106) list
        read(10, 106) list
! for receiver
        read(10,  * ) r0
        read(10,  * ) z0
!	x0 = x0*1e-3
!	y0 = y0*1e-3
!	z0 = z0*1e-3
! for source
        read(10, 106) list
        read(10, 106) list
        read(10, 106) list
        read(10,  * ) M0
!	  read(10,  * ) mxx, myy, mzz, mxy, mxz, myz
        read(10,  * ) zs
!	xs = xs*1e-3
!	ys = ys*1e-3
!	zs = zs*1e-3
!        write(*,*)M0,dip,rake,str,xs,ys,zs
! for basic parameters used in computation
        read(10, 106) list    
        read(10, 106) list    
        read(10, 106) list
        read(10, 106) list
!c    (1). Basic control parameters:
        read(10, 101) m
        read(10,  * ) Twin
!c    (2). Source_type ('E':explosive_source, 'S':single_force,'D'double_couple)
        read(10, 106) list
        read(10, 106) list
        read(10, 105) SourceType
!c    (2). Source parameters:
        read(10, 106) list
        read(10, 106) list
        read(10, 105) Type_STF
        read(10, 106) list
        read(10,  * ) tou
        read(10,  * ) shift ! modified shift as the time windows start time
        read(10,  * ) fc
!c    (3). Window swith and parameters:
        read(10, 106) list
        read(10, 106) list
        read(10, 105) WinSwitch
        read(10, 105) WinType
        read(10,  * ) f1
        read(10,  * ) f2
        read(10,  * ) f3
        read(10,  * ) f4
!c    (4). Output Format:
        read(10, 106) list
        read(10, 106) list
        read(10, 105) Outname
! multilayers 
        read(10, 106) list
        read(10, 106) list
        read(10, 106) list
        read(10, 101) nly
        read(10, 106) list
        do i=1,nly
           read(10, *) j,z(i-1),rho(i),vs0(i),vp0(i),qs(i),qp(i)
        end do
        close (10)

!c       .....
!c       double-couple source:          
!c       .....
!c       DS0  -- slip amount;
!c       A0S  -- size of fault plane;
!c       zs   -- depth of source;
!c       dip  -- dip angle of fault plane;
!c       str  -- strike angle of fault;
!c       rak  -- rake angle;
!c       ....................
!c        print*,'open file in_s_d'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        if (M0.le.12) then
! when M0 is less than 12, M0 should be the magnitude of earthquake other than moment
! so, should be transform into scalor moment(dyne cm) according to Kanamori's equation
! (JGR, 1977,82,2981-2987)
!            M0 = 10**(1.5*M0+16.1)*1d-22
            M0 = 10**(1.5*M0-5.9)
        else
            M0 = M0*1d-22
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        


!c  CHECKING INPUTS:
!c  ---------------------

        ls   = layernumb(zs,nly,z)
        lo   = layernumb(z0,nly,z)
       
        mus = (vs0(ls)**2)*rho(ls)
!c        M0  = D0s*A0s*mus 

        if (nly.lt.ls.or.ls.lt.1) then
             print*,'Error: info@input-checking_1'
             stop
        end if
          
        if (nly.eq.ls) then                   
!c# introduce a fictitious interface to prevent that source is within 
!c  the last layer.
            if (zs.gt.z(nly-1)) then          
               print *, ' Add one fictious interface '
               z(nly)     = zs + 5.0  !!!!!!!!!                
               vs0(nly+1) = vs0(nly)   
               vp0(nly+1) = vp0(nly)
               qs(nly+1)  = qs(nly)
               qp(nly+1)  = qp(nly)
               rho(nly+1) = rho(nly)
               nly=nly+1
               print*,'nly=',nly,'     z(nly-1)=',z(nly-1)
             else
               print*,'Error: info@input-checking_2'
               stop
             end if
          end if

          if (nly .gt. NM) then
             print*,'nly>NM !!!'
             print*,'Current NM = ', NM, '  Current layer number=', nly
             stop
          endif

        END






