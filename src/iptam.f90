      subroutine ptam(no,nf,k,inte,counter,o)
      
      use comvar      
      include   "green_com.inc"
      
      integer :: no, nf, counter
      real*8  :: k, kk(3)
      real*8  :: xm, fm, igd
      complex*16 :: inte
! inte : one componet of displacement      
      
      if (counter <= 3) then
! initiate first three elements for peak and trough search      
            rp(no,counter)=real(inte)
            ip(no,counter)=imag(inte)
      else
! go ahead and update elements for peak and trough search
            rp(no,1)=rp(no,2)
            rp(no,2)=rp(no,3)
            rp(no,3)=real(inte)
            ip(no,1)=ip(no,2)
            ip(no,2)=ip(no,3)
            ip(no,3)=imag(inte)
            kk(1)=k-2d0*dk
            kk(2)=k-dk
            kk(3)=k
            
!-------------- real part -------------------            
            if(rp(no,2) == rp(no,1) .and. rp(no,2) == rp(no,3)) then
! if all three elements are equal, average is unneccesary
                  sumreal(no)=rp(no,2)
                  sumimag(no)=ip(no,2)
                  flagrp(no)=0
                  flagip(no)=0
! for real peak                  
            else if(rp(no,2) >= rp(no,1) .and. rp(no,2) > rp(no,3) .and. flag(no,1) == 1) then
                  num(no,1)=num(no,1)+1     ! counter the number of peak or trough
                  if(num(no,1) == MAX) flag(no,1)=0
                  ! fit three point and get coordinate of true peak
                  call extremum(rp)
                  pre(no,num(no,1),1)=xm
                  pre(no,num(no,1),2)=fm
            else if(rp(no,2) > rp(no,1) .and. rp(no,2) == rp(no,3)) then
                  call igdcalc(1)
                  if(rp(no,2) == rp(no,3) .and. igd == 0d0) then
                        sumreal(no)=rp(no,2)
                        sumimag(no)=ip(no,2)
                        flagrp(no)=0
                        flagip(no)=0
                  else if(flag(no,1) ==1 ) then
                        num(no,1)=num(no,1)+1
                        if(num(no,1) == MAX) flag(no,1)=0
                        call extremum(rp)
                        pre(no,num(no,1),1)=xm
                        pre(no,num(no,1),2)=fm
                  end if
! for real trough                  
            else if(rp(no,2) < rp(no,1) .and. rp(no,2) == rp(no,3)) then
                  call igdcalc(1)
                  if(rp(no,2) ==rp(no,3) .and. igd == 0d0) then
                        sumreal(no)=rp(no,2)
                        sumimag(no)=ip(no,2)
                        flagrp(no)=0
                        flagip(no)=0
                  else if(flag(no,1) ==1 ) then
                        num(no,1)=num(no,1)+1
                        if(num(no,1) == MAX) flag(no,1)=0
                        call extremum(rp)
                        pre(no,num(no,1),1)=xm
                        pre(no,num(no,1),2)=fm
                  end if
            else if(rp(no,2) < rp(no,1) .and. rp(no,2) < rp(no,3) .and. flag(no,2) ==1 ) then
                  num(no,2)=num(no,2)+1
                  if(num(no,2) == MAX) flag(no,2)=0
                  call extremum(rp)
                  tre(no,num(no,2),1)=xm
                  tre(no,num(no,2),2)=fm
            else if(rp(no,2) < rp(no,1) .and. rp(no,2) == rp(no,3)) then
                  call igdcalc(1)
                  if(rp(no,2) == rp(no,3) .and. igd == 0d0) then
                        sumreal(no)=rp(no,2)
                        sumimag(no)=ip(no,2)
                        flagrp(no)=0
                        flagip(no)=0
                  else if(flag(no,2) == 1) then
                        num(no,2)=num(no,2)+1
                        if(num(no,2) == MAX) flag(no,2)=0
                        call extremum(rp)
                        tre(no,num(no,2),1)=xm
                        tre(no,num(no,2),2)=fm
                  end if
            end if

!-------------- imaginary part -------------------
! for peak
            if(ip(no,2) == ip(no,1) .and. ip(no,2) == ip(no,3) .and. flagip(no) ==1) then
                  sumimag(no)=ip(no,2)
                  flagip(no)=0
            else if(ip(no,2) >= ip(no,1) .and. ip(no,2) > ip(no,3) .and. flagip(no) ==1 .and. flag(no,3) ==1) then
                  num(no,3)=num(no,3)+1
                  if(num(no,3) == MAX) flag(no,3)=0
                  call extremum(ip)
                  pim(no,num(no,3),1)=xm
                  pim(no,num(no,3),2)=fm
            else if(ip(no,2) > ip(no,1) .and. ip(no,2) == ip(no,3) .and. flagip(no) ==1) then
                  call igdcalc(2)
                  if(ip(no,2) == ip(no,3) .and. igd == 0d0) then
                        sumimag(no)=ip(no,2)
                        flagip(no)=0
                  else if(flag(no,3) == 1) then
                        num(no,3)=num(no,3)+1
                        if(num(no,3) == MAX) flag(no,3)=0
                        call extremum(ip)
                        pim(no,num(no,3),1)=xm
                        pim(no,num(no,3),2)=fm
                  end if
! for trough                  
            else if(ip(no,2) < ip(no,1) .and. ip(no,2) < ip(no,3) .and. flagip(no) == 1 .and. flag(no,4) == 1) then
                  num(no,4)=num(no,4)+1
                  if(num(no,4) == MAX) flag(no,4)=0
                  call extremum(ip)
                  tim(no,num(no,4),1)=xm
                  tim(no,num(no,4),2)=fm
            else if(ip(no,2) < ip(no,1) .and. ip(no,2) == ip(no,3) .and. flagip(no) ==1) then
                  call igdcalc(2)
                  if(ip(no,2) == ip(no,3) .and. igd == 0d0) then
                        sumimag(no)=ip(no,2)
                        flagip(no)=0
                  else if(flag(no,4) ==1) then
                        num(no,4)=num(no,4)+1
                        if(num(no,4) == MAX) flag(no,4)=0
                        call extremum(ip)
                        tim(no,num(no,4),1)=xm
                        tim(no,num(no,4),2)=fm
                  end if      
            end if
      end if
      
      if(flagrp(no) == 1 .and. flag(no,1)+flag(no,2) == 0) then
            call averagevalue(pre,tre,sumreal(no))
            flagrp(no)=0
      end if
      
      if(flagip(no) ==1 .and. flag(no,3)+flag(no,4) == 0) then
            call averagevalue(pim,tim,sumimag(no))
            flagip(no)=0
      end if
!      write(*,*)sumreal(no)
      contains
            subroutine extremum(f)
            real*8 :: a1, a2
! 10 integ            
            real*8 :: f(10,3)
            
            a1=2d0*f(no,3)-4d0*f(no,2)+2d0*f(no,1)
            a2=4d0*f(no,2)-f(no,3)-3d0*f(no,1)
            xm=kk(1)-a2*(kk(3)-kk(1))/(2d0*a1)
            fm=f(no,1)-a2*a2/(4d0*a1)
            if (a1 == 0d0) then
                xm = kk(2)
                fm = f(no,2)
            end if
!                write(*,*)a1,a2,fm
            end subroutine extremum
            
            subroutine averagevalue(peak,trough,aver)
	    
            integer :: i, count
! 10 integ	      
            real*8  :: peak(10, MAX, 2), trough(10, MAX, 2)
            real*8  :: s(MAX*2)
            real*8  :: aver
            
            if(peak(no,1,1) < trough(no,1,1)) then
                  do i=1,MAX
                        s(2*i-1) = peak(no,i,2)
                        s(2*i  ) = trough(no,i,2)
                  end do
            else
                  do i=1,MAX
                        s(2*i-1) = trough(no,i,2)
                        s(2*i  ) = peak(no,i,2)
                  end do
            end if
            do count=10,2,-1
                  do i=1,count-1
                        s(i)=(s(i)+s(i+1))/2d0
                  end do
            end do
            aver=s(1)
!            write(*,*)aver
            end subroutine averagevalue
            
            subroutine igdcalc(judge)
            
            integer :: judge, i
! 10 integ            
            complex*16 :: igd0(10)
            
            call integrand_calc(k, o, igd0)
                        
            do i=1,3
                  if(i == no .and. judge == 1) then
                        igd=real(igd0(i))
                  end if
                  if(i == no .and. judge == 2) then
                        igd=imag(igd0(i))
                  end if
            end do
            
            end subroutine igdcalc
            
      end subroutine ptam             
