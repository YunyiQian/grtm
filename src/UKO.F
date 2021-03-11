      subroutine uko_D
!c     .......................

      include     "green_com.inc"
      complex*16  a0(2), a1(2), a2(2), b1, b2
      common      /ab012/a0, a1, a2, b1, b2

!cccccccccccc Case: j=s and z0<zs or j<s  cccccccccccccc
      if (lo.eq.ls.and.z0.le.zs.or.lo.lt.ls) then
           b1=yu0*su0_1
!           b2=yu0*su0_2
           do i=1,2
              a0(i)=yu(i,1)*su_0(1)+yu(i,2)*su_0(2)
              a1(i)=yu(i,1)*su_1(1)+yu(i,2)*su_1(2)
!              a2(i)=yu(i,1)*su_2(1)+yu(i,2)*su_2(2)
           end do
!cccccccccccc Case: j=s and z0>zs or j>s  cccccccccccccc
        else if(lo.eq.ls.and.z0.gt.zs.or.lo.gt.ls) then
           b1=yd0*sd0_1
!           b2=yd0*sd0_2
           do i=1,2
              a0(i)=yd(i,1)*sd_0(1)+yd(i,2)*sd_0(2)
              a1(i)=yd(i,1)*sd_1(1)+yd(i,2)*sd_1(2)
!              a2(i)=yd(i,1)*sd_2(1)+yd(i,2)*sd_2(2)
           end do
      end if
!c     ................
      return
      end 
