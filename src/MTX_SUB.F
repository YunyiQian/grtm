
      subroutine mtxmtp(n, a, b, c)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Note: The purpose of this subroutine is to obtain matrix   c
!c         c, which is the result of matrix a multiply matix b. c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         
      complex*16 zero, a(4,4), b(4,4), c(4,4)
      integer    n, i, j, k

      zero = 0.0
      do j = 1 ,n
         do i = 1, n
            c(i, j) = zero
         end do
         do k = 1, n
            do i = 1, n
               c(i, j) = c(i,j) + a(i, k)*b(k, j)
            end do
         end do
      end do

      return
      end 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mtxinv_2(a0)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Note: The purpose of this subroutine is to obtain the inverse      c
!c         matrix of a given 4*4 matrix a0, note that only the elements c
!c         of first 2 lines and first 2 columns are nonzero. The result c
!c         is restored in the original matrix a0.                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         
      real*8         epsi
      parameter      ( epsi=1.0E-20 )
      complex*16     d, a(4, 4)
      complex*16     a0(4, 4)
      
      a(1,1) = a0(1,1)
      a(2,1) = a0(2,1)
      a(1,2) = a0(1,2)
      a(2,2) = a0(2,2)

      d = a(1,1)*a(2,2)-a(1,2)*a(2,1)
      if ( cdabs(d).lt.epsi ) then
        print*,' Singular Matrix_2'
        d = epsi
      endif

      a0(1,1) = a(2,2)/d
      a0(2,2) = a(1,1)/d
      a0(1,2) =-a(1,2)/d
      a0(2,1) =-a(2,1)/d
        
      return
      end 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine mtxinv_4(a)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Note: The purpose of this subroutine is to obtain the inverse      c
!c         matrix of a given 4*4 matrix a0. The result is restored in   c
!c               the original matrix a0.                                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        implicit none
        real*8        epsi,d
        parameter     (epsi=1.0e-10)
        integer       i,j,k,l(4),m(4)
        complex*16    a(4,4),t

        do k=1,4
           d=0.0
           do i=k,4
            do j=k,4
              if(cdabs(a(i,j)).gt.d) then
                d=cdabs(a(i,j))
                l(k)=i
                m(k)=j
              end if
                end do
           end do
           if(d.lt.epsi) then
             print*,d
             print*,'The matrix is singlur!'
             stop
           end if

          if(l(k).ne.k) then
           do j=1,4
             t=a(k,j)
              a(k,j)=a(l(k),j)
              a(l(k),j)=t
           end do
          end if
          if(m(k).ne.k) then
           do i=1,4
             t=a(i,k)
             a(i,k)=a(i,m(k))
             a(i,m(k))=t
           end do
          end if
        
           a(k,k)=1./a(k,k)
        
           do i=1,4
             if(i.ne.k) then
                  a(i,k)=-a(k,k)*a(i,k)
             end if
           end do

           do i=1,4
             do j=1,4
                if(i.ne.k.and.j.ne.k) then
                   a(i,j)=a(i,j)+a(i,k)*a(k,j)
                end if
             end do
           end do
           do j=1,4
             if(j.ne.k) then
                 a(k,j)=a(k,j)*a(k,k)
              end if
           end do

        end do

        do k=4,1,-1
           do j=1,4
             t=a(k,j)
             a(k,j)=a(m(k),j)
             a(m(k),j)=t
        end do
            do i=1,4
             t=a(i,k)
             a(i,k)=a(i,l(k))
             a(i,l(k))=t
            end do
        end do

        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE MTXE(LAY,KN,O,E)
        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Note: The purpose of this subroutine is to calculate the E matrix  c
!c         of the lay_th layer for P-SV case, cf. eq.(3-18).            c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INCLUDE     "green_com.inc"
        complex*16  E(4,4), C

        C=KN*KN+CSN(LAY)*CSN(LAY)
C       .......

        E(1,1)= VP(LAY)*KN/O
        E(1,2)= VS(LAY)*CSN(LAY)/O
        E(1,3)= E(1,1)
        E(1,4)= E(1,2)
        E(2,1)= VP(LAY)*CPN(LAY)/O
        E(2,2)= VS(LAY)*KN/O
        E(2,3)=-E(2,1)
        E(2,4)=-E(2,2)
        E(3,1)=-2*VP(LAY)*MU(LAY)*KN*CPN(LAY)/O
        E(3,2)=-VS(LAY)*MU(LAY)*C/O
        E(3,3)=-E(3,1)
        E(3,4)=-E(3,2)
        E(4,1)=-VP(LAY)*MU(LAY)*C/O
        E(4,2)=-2*VS(LAY)*MU(LAY)*KN*CSN(LAY)/O
        E(4,3)= E(4,1)
        E(4,4)= E(4,2)          
C       ...........
        RETURN
        END     
