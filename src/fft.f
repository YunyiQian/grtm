!c########################################################

       subroutine fft(a,m,flag)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      ��ʱ���ȡ��FFT�ӳ���:
!c      a: �������飬����Ϊ2**m�����Ҫ����fft�����У�
!c      m: ���г���Ϊ 2**m.
!c      flag: FFT��IFFT��ѡ�������
!c            flag=1  -------  FFT;
!c            flag=-1 ------- IFFT.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
	integer    l,m,n,n2,i,j,k,l1,l2,ip,flag
	real*8       pi
      complex*16   a(2**m),u,w,t

	pi=4d0*atan(1d0)
	n=2**m                 !   ��������
	n2=n/2            

	if(flag.eq.-1) then    ! 
	  do i=1,n             !
	   a(i)=conjg(a(i))    !   ��ifftʱ������aȡ����
        end do               !
	end if                 !

!ccccccccccccccccccc  ��λ������  ccccccccccccccccccccc

	j=1
	do i=1,n-1
	   if(i.lt.j) then     ! 
	       t=a(j)          !
	       a(j)=a(i)       !  ���i<j������a(i)��a(i)
	       a(i)=t          !
	   end if              !

	   k=n2                !     
	   do while(k.lt.j)    !      
	       j=j-k           !  ��֪j,�÷���ӷ�
	       k=k/2           !  ����һ����λ���
	   end do              !  
	   j=j+k               !
      end do
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do l=1,m              ! ��һ��ѭ����lΪ��������1��m
	 l1=2**l              ! l1 : ��l���еĳ������� 
	 l2=l1/2              ! l2 : ��l���ĵ��μ��
	 w=cmplx(cos(pi/float(l2)),-sin(pi/float(l2)))
	                      ! w  : ��l���ĳ�������
	 u=cmplx(1.0,0.0)
	 do j=1,l2            ! �ڶ���ѭ������������
	   do i=j,n,l1        ! ������ѭ����Ⱥ����
	      ip=i+l2         ! i:�������Ͻ���ţ�ip:���½����
	      t=a(ip)*u
	      a(ip)=a(i)-t    ! ��������
	      a(i)=a(i)+t
	   end do
	   u=u*w              ! �ڶ���ѭ���е�������
	 end do 
	end do

	if(flag.eq.-1) then   ! ��ifftʱ��������aȡ���������n
	  do i=1,n
	   a(i)=conjg(a(i))/n
        end do
	end if     

	return
	end 
	     
