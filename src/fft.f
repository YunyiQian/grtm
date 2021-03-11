!c########################################################

       subroutine fft(a,m,flag)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      按时间抽取的FFT子程序:
!c      a: 复数数组，长度为2**m，存放要进行fft的序列；
!c      m: 序列长度为 2**m.
!c      flag: FFT或IFFT的选择参数：
!c            flag=1  -------  FFT;
!c            flag=-1 ------- IFFT.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
	integer    l,m,n,n2,i,j,k,l1,l2,ip,flag
	real*8       pi
      complex*16   a(2**m),u,w,t

	pi=4d0*atan(1d0)
	n=2**m                 !   抽样点数
	n2=n/2            

	if(flag.eq.-1) then    ! 
	  do i=1,n             !
	   a(i)=conjg(a(i))    !   做ifft时对数组a取共轭
        end do               !
	end if                 !

!ccccccccccccccccccc  倒位序程序段  ccccccccccccccccccccc

	j=1
	do i=1,n-1
	   if(i.lt.j) then     ! 
	       t=a(j)          !
	       a(j)=a(i)       !  如果i<j，互换a(i)和a(i)
	       a(i)=t          !
	   end if              !

	   k=n2                !     
	   do while(k.lt.j)    !      
	       j=j-k           !  已知j,用反向加法
	       k=k/2           !  求下一个倒位序号
	   end do              !  
	   j=j+k               !
      end do
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do l=1,m              ! 第一层循环，l为级数，从1到m
	 l1=2**l              ! l1 : 第l级中的乘数个数 
	 l2=l1/2              ! l2 : 第l级的蝶形间隔
	 w=cmplx(cos(pi/float(l2)),-sin(pi/float(l2)))
	                      ! w  : 第l级的乘数底数
	 u=cmplx(1.0,0.0)
	 do j=1,l2            ! 第二层循环，乘数控制
	   do i=j,n,l1        ! 第三层循环，群控制
	      ip=i+l2         ! i:蝶形左上角序号，ip:左下脚序号
	      t=a(ip)*u
	      a(ip)=a(i)-t    ! 蝶形运算
	      a(i)=a(i)+t
	   end do
	   u=u*w              ! 第二层循环中调整乘数
	 end do 
	end do

	if(flag.eq.-1) then   ! 做ifft时对数组再a取共轭，并除以n
	  do i=1,n
	   a(i)=conjg(a(i))/n
        end do
	end if     

	return
	end 
	     
